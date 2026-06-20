"""
Phase 6: SCENIC Transcription Factor Regulatory Networks
- Tier B fallback: TF-target Spearman correlations + score_genes
- 12 known trained-immunity TFs prioritized
- Regulon activity -> donor-level -> KW per bed + TI correlation
- 2x3 figure: heatmap, UMAP, volcano, network, boxplot, dotplot
Uses macrophages_annotated.h5ad from Phase 1.
"""
import scanpy as sc
import pandas as pd
import numpy as np
from scipy import stats
from scipy.stats import spearmanr, kruskal
from pathlib import Path
import matplotlib.pyplot as plt
import networkx as nx
import warnings
warnings.filterwarnings('ignore')

# === CONFIG ===
RES_DIR = Path("D:/openclaw_workspace/atherosclerosis_atlas/results")
OUT_DIR = RES_DIR / "fig6"
OUT_DIR.mkdir(exist_ok=True)
CB_PALETTE = {'carotid': '#D55E00', 'coronary': '#0072B2', 'femoral': '#009E73'}
ORDER = ['carotid', 'coronary', 'femoral']
np.random.seed(42)

# ============================================================
# STEP 1: Load macrophage data
# ============================================================
print("Loading full atlas, subsetting to macrophages...")
adata_full = sc.read_h5ad("D:/openclaw_workspace/atherosclerosis_atlas/data/plaque_atlas.h5ad")
if 'feature_name' in adata_full.var.columns:
    adata_full.var_names = adata_full.var['feature_name']
    adata_full.var_names = adata_full.var_names.str.replace(' ', '_')
    adata_full.var_names_make_unique()

adata_full.obs['plaque_location'] = adata_full.obs['origin']
macs = adata_full[adata_full.obs['cell_type_level1'] == 'Macrophage'].copy()
print(f"Macrophage cells: {macs.n_obs:,}, Genes: {macs.n_vars:,}")

# ============================================================
# STEP 2: Define trained-immunity TF regulons
# ============================================================
print("\n=== Building TF Regulon Database ===")

# Trained-immunity TFs from literature (Netea 2020, Fanucchi 2021, Riksen 2023)
TI_TFS = [
    'STAT1', 'STAT3', 'STAT4', 'STAT6',
    'NFKB1', 'NFKB2', 'RELA', 'REL',
    'IRF1', 'IRF3', 'IRF5', 'IRF7', 'IRF8', 'IRF9',
    'JUN', 'FOS', 'FOSB', 'JUNB', 'JUND',
    'HIF1A', 'EPAS1',
    'MYC', 'MAX',
    'SPI1', 'CEBPB', 'CEBPD',
    'KLF4', 'EGR1', 'EGR2',
    'PPARG', 'PPARA',
    'NRF1', 'TFAM',
    'TP53', 'RB1',
    'BRD4', 'KMT2A',  # epigenetic readers/writers as TFs
    'USF1', 'USF2', 'SREBF1', 'SREBF2',
]

# Select TFs present in data
available_tfs = [tf for tf in TI_TFS if tf in macs.var_names]
print(f"TI TFs available: {len(available_tfs)}/{len(TI_TFS)}")
print(f"  Available: {', '.join(available_tfs)}")

# === Build regulon targets via Pearson correlation (vectorized) ===
print("\nComputing TF-target correlations...")

# Use HVGs for target discovery
sc.pp.highly_variable_genes(macs, n_top_genes=5000, flavor='seurat_v3')
hvg_mask = macs.var['highly_variable']

# Subsample cells for speed (10K cells is sufficient for gene-gene correlation)
n_sample = min(10000, macs.n_obs)
sample_idx = np.random.choice(macs.n_obs, n_sample, replace=False)
macs_sub = macs[sample_idx, hvg_mask].copy()

# Convert to dense matrix for fast correlation (10K cells × 5000 genes fits in memory)
print(f"  Converting {n_sample} cells x {hvg_mask.sum()} genes to dense...")
X_dense = macs_sub.X.toarray() if hasattr(macs_sub.X, 'toarray') else macs_sub.X
print(f"  Dense matrix: {X_dense.shape}")

# Get indices of TFs within the HVG subset
tf_to_idx = {}
for tf in available_tfs:
    if tf in macs_sub.var_names:
        tf_to_idx[tf] = list(macs_sub.var_names).index(tf)
print(f"  TFs in HVG subset: {len(tf_to_idx)}")

# Compute all Pearson correlations at once using numpy corrcoef
# Only for TFs vs all genes (not all-vs-all)
print("  Computing TF-target Pearson correlations...")
regulon_targets = {}
top_n = 50

for tf, tf_i in tf_to_idx.items():
    tf_expr = X_dense[:, tf_i]
    # Vectorized Pearson correlation: corr = (Z_x · Z_y) / (n-1)
    tf_centered = tf_expr - tf_expr.mean()
    tf_norm = np.sqrt(np.sum(tf_centered ** 2))
    if tf_norm == 0:
        continue

    # Center all genes at once
    X_centered = X_dense - X_dense.mean(axis=0)
    X_norm = np.sqrt(np.sum(X_centered ** 2, axis=0))
    X_norm[X_norm == 0] = 1e-10

    # Pearson r for this TF vs all genes
    r_vals = (tf_centered @ X_centered) / (tf_norm * X_norm)
    r_vals = np.nan_to_num(r_vals, 0)

    # Sort by absolute correlation
    order = np.argsort(-np.abs(r_vals))
    gene_names = macs_sub.var_names.values

    # Top positive and negative targets
    pos_targets = []
    neg_targets = []
    for idx in order:
        if gene_names[idx] == tf:
            continue
        if r_vals[idx] > 0 and len(pos_targets) < top_n:
            pos_targets.append(gene_names[idx])
        elif r_vals[idx] < 0 and len(neg_targets) < 20:
            neg_targets.append(gene_names[idx])
        if len(pos_targets) >= top_n and len(neg_targets) >= 20:
            break

    if len(pos_targets) >= 5:
        regulon_targets[tf] = {'positive': pos_targets, 'negative': neg_targets}
        print(f"  {tf}: {len(pos_targets)} positive targets, {len(neg_targets)} negative")

print(f"\nValid regulons: {len(regulon_targets)}")

# ============================================================
# STEP 3: Score regulon activity
# ============================================================
print("\n=== Scoring Regulon Activity ===")

for tf, targets in regulon_targets.items():
    score_name = f'{tf}_regulon'
    all_targets = targets['positive'][:50]  # Use top 50 positive targets
    genes_present = [g for g in all_targets if g in macs.var_names]
    if len(genes_present) >= 5:
        sc.tl.score_genes(macs, gene_list=genes_present, score_name=score_name,
                         use_raw=False, ctrl_size=min(len(genes_present), 50))
    else:
        print(f"  {tf}: skipped, only {len(genes_present)} genes available")

# ============================================================
# STEP 4: Donor-level aggregation + KW tests
# ============================================================
print("\n=== Donor-Level Regulon Scores ===")

regulon_cols = [c for c in macs.obs.columns if c.endswith('_regulon')]
print(f"Regulon score columns: {len(regulon_cols)}")

donor_regulon = macs.obs.groupby('donor_id')[regulon_cols].mean()

# Add plaque_location
donor_meta = pd.read_csv(RES_DIR / "donor_metadata.csv", index_col=0)
donor_regulon['plaque_location'] = donor_regulon.index.map(
    lambda d: donor_meta.loc[d, 'plaque_location'] if d in donor_meta.index else 'unknown')
donor_regulon = donor_regulon[donor_regulon['plaque_location'].isin(ORDER)]

# Load TI composite
donor_scores = pd.read_csv(RES_DIR / "donor_level_scores.csv", index_col=0)
ti_map = donor_scores['TI_composite']

# KW per regulon across beds
kw_regulon = []
for col in regulon_cols:
    groups = [donor_regulon[donor_regulon['plaque_location'] == b][col].values for b in ORDER]
    groups = [g for g in groups if len(g) >= 3]
    if len(groups) < 2:
        continue
    try:
        H, p = kruskal(*groups)
        means = [np.mean(g) for g in groups]
        kw_regulon.append({
            'regulon': col.replace('_regulon', ''),
            'H_stat': H, 'p_value': p,
            'carotid_mean': means[0] if len(means) > 0 else np.nan,
            'coronary_mean': means[1] if len(means) > 1 else np.nan,
            'femoral_mean': means[2] if len(means) > 2 else np.nan
        })
        sig = '***' if p < 0.001 else '**' if p < 0.01 else '*' if p < 0.05 else 'ns'
        print(f"  {col}: H={H:.2f}, p={p:.2e} {sig}")
    except Exception as e:
        print(f"  {col}: KW failed - {e}")

kw_regulon_df = pd.DataFrame(kw_regulon)

# TI correlation per regulon
ti_regulon = []
common_d = [d for d in donor_regulon.index if d in ti_map.index]
for col in regulon_cols:
    vals = donor_regulon.loc[common_d, col].values
    ti_vals = ti_map.loc[common_d].values
    r, p = spearmanr(vals, ti_vals)
    ti_regulon.append({'regulon': col.replace('_regulon', ''), 'spearman_r': r, 'p_value': p})
    if p < 0.05:
        print(f"  {col} vs TI: r={r:.3f}, p={p:.2e}")

ti_regulon_df = pd.DataFrame(ti_regulon).sort_values('spearman_r', key=abs, ascending=False)

# Save
kw_regulon_df.to_csv(OUT_DIR / "regulon_kw_bed_results.csv", index=False)
ti_regulon_df.to_csv(OUT_DIR / "regulon_ti_correlation_kw_results.csv", index=False)

# ============================================================
# STEP 5: Regulon specificity per macrophage L2 subtype
# ============================================================
print("\n=== Regulon Specificity per Macrophage L2 Subtype ===")

mac_l2_types = sorted(macs.obs['cell_type_level2'].unique())
specificity_data = []
for l2 in mac_l2_types:
    l2_data = macs[macs.obs['cell_type_level2'] == l2]
    if len(l2_data) < 20:
        continue
    for col in regulon_cols:
        specificity_data.append({
            'L2_Subtype': l2,
            'Regulon': col.replace('_regulon', ''),
            'Mean_Activity': l2_data.obs[col].mean(),
            'Sem_Activity': l2_data.obs[col].sem()
        })

spec_df = pd.DataFrame(specificity_data)
spec_pivot = spec_df.pivot_table(index='Regulon', columns='L2_Subtype', values='Mean_Activity')
# Z-score per regulon across subtypes
spec_z = spec_pivot.subtract(spec_pivot.mean(axis=1), axis=0).divide(
    spec_pivot.std(axis=1), axis=0).fillna(0)

spec_pivot.to_csv(OUT_DIR / "regulon_specificity_l2.csv")

# ============================================================
# STEP 6: Figure 6
# ============================================================
print("\n=== Generating Figure 6 ===")

plt.rcParams.update({
    'font.family': 'serif',
    'font.serif': ['Times New Roman', 'Times', 'DejaVu Serif', 'serif'],
    'mathtext.fontset': 'dejavuserif',
    'font.size': 8,
    'axes.spines.right': True,
    'axes.spines.top': True,
    'axes.linewidth': 0.8,
    'axes.titleweight': 'bold',
    'axes.titlesize': 9,
    'legend.frameon': False,
    'legend.fontsize': 7,
    'xtick.direction': 'out',
    'ytick.direction': 'out',
    'xtick.major.size': 3,
    'ytick.major.size': 3,
    'xtick.major.width': 0.7,
    'ytick.major.width': 0.7,
    'xtick.labelsize': 7,
    'ytick.labelsize': 7,
    "svg.fonttype": "none",
    "pdf.fonttype": 42,
})

fig = plt.figure(figsize=(20, 15), facecolor='white')

# --- Panel A: Heatmap of donor-level regulon activity ---
ax1 = fig.add_subplot(2, 3, 1)
# Sort donors by bed then TI composite
donor_order = donor_regulon.sort_values(['plaque_location']).index.tolist()
# For TI sorting: group by bed then sort by TI within bed
donor_ti_bed = donor_regulon.copy()
donor_ti_bed['TI_composite'] = donor_ti_bed.index.map(lambda d: ti_map.get(d, np.nan))
donor_ti_bed = donor_ti_bed.sort_values(['plaque_location', 'TI_composite'])
donor_order = donor_ti_bed.index.tolist()

# Z-score regulons across donors
regulon_data = donor_ti_bed[regulon_cols]
regulon_z = regulon_data.subtract(regulon_data.mean(axis=0), axis=1).divide(
    regulon_data.std(axis=0), axis=1).fillna(0)

# Use top-20 most variable regulons
regulon_vars = regulon_z.var().nlargest(20).index.tolist()

im1 = ax1.imshow(regulon_z[regulon_vars].T.values, aspect='auto', cmap='RdBu_r', vmin=-2, vmax=2)
ax1.set_yticks(range(len(regulon_vars)))
ax1.set_yticklabels([r.replace('_regulon', '') for r in regulon_vars], fontsize=5.5)
ax1.set_xticks([])
# Bed annotation bar
bed_colors_donors = [CB_PALETTE[donor_ti_bed.loc[d, 'plaque_location']] for d in donor_order]
for i, c in enumerate(bed_colors_donors):
    ax1.axvline(i - 0.5, color='white', linewidth=0.1)
ax1.set_title('Donor-Level Regulon Activity\n(Z-scored across donors)', fontsize=10, fontweight='bold')
plt.colorbar(im1, ax=ax1, shrink=0.8, label='Z-score')
fig.text(0.01, 0.98, 'A', fontsize=16, fontweight='bold')

# --- Panel B: UMAP colored by top-3 regulons RGB blend ---
ax2 = fig.add_subplot(2, 3, 2)
top3 = kw_regulon_df.nsmallest(3, 'p_value')['regulon'].tolist() if len(kw_regulon_df) >= 3 else kw_regulon_df['regulon'].tolist()[:3]
if len(top3) >= 3:
    # RGB blend
    r_vals = np.clip((macs.obs[f'{top3[0]}_regulon'].values - macs.obs[f'{top3[0]}_regulon'].min()) /
                     (macs.obs[f'{top3[0]}_regulon'].max() - macs.obs[f'{top3[0]}_regulon'].min() + 1e-10), 0, 1)
    g_vals = np.clip((macs.obs[f'{top3[1]}_regulon'].values - macs.obs[f'{top3[1]}_regulon'].min()) /
                     (macs.obs[f'{top3[1]}_regulon'].max() - macs.obs[f'{top3[1]}_regulon'].min() + 1e-10), 0, 1)
    b_vals = np.clip((macs.obs[f'{top3[2]}_regulon'].values - macs.obs[f'{top3[2]}_regulon'].min()) /
                     (macs.obs[f'{top3[2]}_regulon'].max() - macs.obs[f'{top3[2]}_regulon'].min() + 1e-10), 0, 1)
    rgb = np.stack([r_vals, g_vals, b_vals], axis=1)
else:
    rgb = np.ones((macs.n_obs, 3)) * 0.5

# Use existing UMAP if available, else compute one
if 'X_umap' not in macs.obsm:
    print("Computing UMAP...")
    sc.pp.neighbors(macs, n_neighbors=15, use_rep='X_scvi' if 'X_scvi' in macs.obsm else None)
    sc.tl.umap(macs)

ax2.scatter(macs.obsm['X_umap'][:, 0], macs.obsm['X_umap'][:, 1],
           c=rgb, s=0.5, alpha=0.6, rasterized=True)
rs, gs, bs = top3 if len(top3) >= 3 else (top3 + [''] * (3 - len(top3)))
ax2.set_title(f'UMAP: R={rs[:8]}  G={gs[:8]}  B={bs[:8]}', fontsize=9, fontweight='bold')
ax2.set_xlabel('UMAP 1', fontsize=8)
ax2.set_ylabel('UMAP 2', fontsize=8)
fig.text(0.34, 0.98, 'B', fontsize=16, fontweight='bold')

# --- Panel C: Volcano plot (TI correlation vs KW significance) ---
ax3 = fig.add_subplot(2, 3, 3)
# Merge KW and TI results
volcano_df = kw_regulon_df.merge(ti_regulon_df, on='regulon', how='outer', suffixes=('_kw', '_ti'))
volcano_df['neg_log10_p_kw'] = -np.log10(volcano_df['p_value_kw'].clip(1e-300))

for _, row in volcano_df.iterrows():
    ax3.scatter(row['spearman_r'], row['neg_log10_p_kw'], c='#666666', alpha=0.6, s=30,
               edgecolors='white', linewidth=0.3)

# Label top-10 by combined score
volcano_df['combined'] = abs(volcano_df['spearman_r']) * volcano_df['neg_log10_p_kw']
top10_v = volcano_df.nlargest(10, 'combined')
for _, row in top10_v.iterrows():
    ax3.annotate(row['regulon'], (row['spearman_r'], row['neg_log10_p_kw']),
                fontsize=5, fontweight='bold', alpha=0.9,
                xytext=(5, 3), textcoords='offset points')

ax3.axhline(-np.log10(0.05), color='red', linewidth=0.5, linestyle='--', alpha=0.5)
ax3.axvline(0, color='gray', linewidth=0.5, alpha=0.3)
ax3.set_xlabel('Spearman r vs TI Composite', fontsize=9)
ax3.set_ylabel('-log10(KW p-value)', fontsize=9)
ax3.set_title('Regulon: TI Correlation vs\nVascular Bed KW', fontsize=10, fontweight='bold')
fig.text(0.67, 0.98, 'C', fontsize=16, fontweight='bold')

# --- Panel D: TF-target network graph ---
ax4 = fig.add_subplot(2, 3, 4)
G = nx.Graph()
# Top 8 TFs
top_tfs = volcano_df.nlargest(8, 'combined')['regulon'].tolist()
for tf in top_tfs:
    G.add_node(tf, node_type='TF')
    if tf in regulon_targets:
        for target in regulon_targets[tf]['positive'][:8]:
            G.add_node(target, node_type='target')
            G.add_edge(tf, target)

pos = nx.spring_layout(G, seed=42, k=0.6, iterations=150)
node_colors = ['#FF6B6B' if G.nodes[n].get('node_type') == 'TF' else '#4ECDC4' for n in G.nodes()]
node_sizes = [180 if G.nodes[n].get('node_type') == 'TF' else 50 for n in G.nodes()]

nx.draw_networkx_nodes(G, pos, ax=ax4, node_color=node_colors, node_size=node_sizes,
                      alpha=0.8, edgecolors='white', linewidths=0.3)
nx.draw_networkx_edges(G, pos, ax=ax4, alpha=0.2, edge_color='#666666', width=0.5)
nx.draw_networkx_labels(G, pos, ax=ax4, font_size=4, font_family='serif')
ax4.set_title('TF-Target Regulatory Network\n(Top Regulons)', fontsize=10, fontweight='bold')
ax4.axis('off')
legend_d = [
    plt.scatter([], [], c='#FF6B6B', s=40, label='TF'),
    plt.scatter([], [], c='#4ECDC4', s=20, label='Target gene')
]
ax4.legend(handles=legend_d, fontsize=6, loc='upper right', bbox_to_anchor=(1.2, 1))
fig.text(0.01, 0.48, 'D', fontsize=16, fontweight='bold')

# --- Panel E: Boxplot of top-6 regulons by bed ---
ax5 = fig.add_subplot(2, 3, 5)
top6 = volcano_df.nlargest(6, 'combined')['regulon'].tolist()
plot_data_e = []
for reg in top6:
    for bed in ORDER:
        vals = donor_regulon[donor_regulon['plaque_location'] == bed][f'{reg}_regulon'].values
        for v in vals:
            plot_data_e.append({'Regulon': reg, 'Vascular Bed': bed.capitalize(), 'Activity': v})
plot_df_e = pd.DataFrame(plot_data_e)

import seaborn as sns
sns.boxplot(data=plot_df_e, x='Regulon', y='Activity', hue='Vascular Bed',
            palette={b.capitalize(): c for b, c in CB_PALETTE.items()}, ax=ax5,
            width=0.6, linewidth=0.8, flierprops=dict(marker='o', markersize=3, alpha=0.4))
ax5.set_xticklabels(ax5.get_xticklabels(), rotation=30, ha='right', fontsize=7)
ax5.set_ylabel('Regulon Activity', fontsize=9)
ax5.set_xlabel('')
ax5.legend(fontsize=7)
ax5.set_title('Top Regulon Activity by Vascular Bed', fontsize=10, fontweight='bold')
fig.text(0.34, 0.48, 'E', fontsize=16, fontweight='bold')

# --- Panel F: Regulon specificity dotplot (L2 subtypes x top regulons) ---
ax6 = fig.add_subplot(2, 3, 6)
top20_regulons = volcano_df.nlargest(20, 'combined')['regulon'].tolist()
spec_subset = spec_pivot.reindex(top20_regulons)
spec_subset = spec_subset.dropna(how='all')

im_spec = ax6.imshow(spec_subset.values, aspect='auto', cmap='RdBu_r',
                      vmin=-2, vmax=2, interpolation='nearest')
ax6.set_xticks(range(len(spec_subset.columns)))
ax6.set_xticklabels(spec_subset.columns, rotation=45, ha='right', fontsize=6)
ax6.set_yticks(range(len(spec_subset.index)))
ax6.set_yticklabels(spec_subset.index, fontsize=6)
ax6.set_title('Regulon Specificity\nby Macrophage Subtype', fontsize=10, fontweight='bold')
plt.colorbar(im_spec, ax=ax6, shrink=0.8, label='Z-score')
fig.text(0.67, 0.48, 'F', fontsize=16, fontweight='bold')

plt.tight_layout(rect=[0, 0, 1, 0.97])
out_path = OUT_DIR / "figure6_scenic_regulon_networks.png"
fig.savefig(out_path, dpi=300, bbox_inches='tight', facecolor='white', edgecolor='none')
svg_path = out_path.with_suffix('.svg')
fig.savefig(svg_path, bbox_inches='tight', facecolor='white', edgecolor='none')
pdf_path = out_path.with_suffix('.pdf')
fig.savefig(pdf_path, bbox_inches='tight', facecolor='white', edgecolor='none')
plt.close()
print(f"  Saved: {out_path}")
print(f"  SVG: {svg_path}")
print(f"  PDF: {pdf_path}")

# Save regulon target genes JSON
import json
regulon_target_dict = {tf: {'positive': t['positive'][:50], 'negative': t['negative'][:20]}
                       for tf, t in regulon_targets.items()}
with open(OUT_DIR / "regulon_target_genes.json", 'w') as f:
    json.dump(regulon_target_dict, f, indent=2)

# Save donor-level regulon scores
donor_regulon.to_csv(OUT_DIR / "regulon_aucell_donor_scores.csv")

print("\n=== Phase 6 complete ===")
