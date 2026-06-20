"""
Phase 4: Full-Atlas Cell Type Proportions & Vascular Bed Differences
— Donor-level L1/L2 composition across all 13 cell types + 23 subtypes
— Vascular bed KW tests with permutation validation for femoral (n=7)
— Links global cell type shifts to TI composite (Phase 2 output)
Uses output from phase1_macrophage_analysis.py + full plaque_atlas.h5ad
"""

import scanpy as sc
import pandas as pd
import numpy as np
from scipy import stats
from scipy.stats import spearmanr, kruskal
from pathlib import Path
import matplotlib.pyplot as plt
import seaborn as sns
import warnings
warnings.filterwarnings('ignore')

# === CONFIG ===
RES_DIR = Path("D:/openclaw_workspace/atherosclerosis_atlas/results")
OUT_DIR = RES_DIR / "fig4"
OUT_DIR.mkdir(exist_ok=True)

sc.set_figure_params(dpi=100, facecolor='white', frameon=True)
CB_PALETTE = {'carotid': '#D55E00', 'coronary': '#0072B2', 'femoral': '#009E73'}
ORDER = ['carotid', 'coronary', 'femoral']
np.random.seed(42)

# ============================================================
# STEP 1: Load full atlas
# ============================================================
print("Loading full atlas...")
adata = sc.read_h5ad("D:/openclaw_workspace/atherosclerosis_atlas/data/plaque_atlas.h5ad")
print(f"Cells: {adata.n_obs:,}, Genes: {adata.n_vars:,}")

# Convert var_names to gene symbols (same as Phase 1)
if 'feature_name' in adata.var.columns:
    adata.var_names = adata.var['feature_name']
    adata.var_names = adata.var_names.str.replace(' ', '_')
    adata.var_names_make_unique()
    print("var_names converted to gene symbols")

# Create canonical columns
adata.obs['plaque_location'] = adata.obs['origin']
donor_meta = pd.read_csv(RES_DIR / "donor_metadata.csv", index_col=0)
print(f"Donors: {adata.obs['donor_id'].nunique()}")
print(f"L1 types: {adata.obs['cell_type_level1'].nunique()}")
print(f"L2 types: {adata.obs['cell_type_level2'].nunique()}")

# ============================================================
# STEP 2: Donor-level cell type composition
# ============================================================
print("\n=== Donor-Level Cell Type Composition ===")

# L1 proportions per donor
l1_counts = adata.obs.groupby(['donor_id', 'cell_type_level1'], observed=False).size().unstack(fill_value=0)
l1_props = l1_counts.div(l1_counts.sum(axis=1), axis=0)

# L2 proportions per donor
l2_counts = adata.obs.groupby(['donor_id', 'cell_type_level2'], observed=False).size().unstack(fill_value=0)
l2_props = l2_counts.div(l2_counts.sum(axis=1), axis=0)

# Add metadata (plaque_location from donor_meta)
donor_meta['plaque_location'] = donor_meta['plaque_location'].astype(str)
l1_props['plaque_location'] = l1_props.index.map(
    lambda d: donor_meta.loc[d, 'plaque_location'] if d in donor_meta.index else 'unknown')
l2_props['plaque_location'] = l2_props.index.map(
    lambda d: donor_meta.loc[d, 'plaque_location'] if d in donor_meta.index else 'unknown')

# Filter to known beds
l1_props = l1_props[l1_props['plaque_location'].isin(ORDER)]
l2_props = l2_props[l2_props['plaque_location'].isin(ORDER)]

# Remove 'Undefined' L2 if present
if 'Undefined' in l2_props.columns:
    l2_props = l2_props.drop(columns=['Undefined'])

cell_types_l1 = [c for c in l1_props.columns if c != 'plaque_location']
cell_types_l2 = [c for c in l2_props.columns if c != 'plaque_location']

print(f"L1 cell types: {len(cell_types_l1)}")
print(f"L2 cell types: {len(cell_types_l2)}")

# ============================================================
# STEP 3: KW tests + eta-squared for cell type proportions
# ============================================================
print("\n=== KW Tests: Cell Type Proportions by Vascular Bed ===")

def kruskal_eta2(groups):
    """Eta-squared effect size for Kruskal-Wallis."""
    k = len(groups)
    N = sum(len(g) for g in groups)
    all_vals = np.concatenate(groups)
    grand_median = np.median(all_vals)
    ssb = sum(len(g) * (np.median(g) - grand_median)**2 for g in groups)
    sst = sum((all_vals - grand_median)**2)
    return ssb / sst if sst > 0 else 0.0

def permutation_kw(groups, n_perm=10000):
    """Permutation-validated KW p-value (robust for small n)."""
    all_vals = np.concatenate(groups)
    labels = np.concatenate([np.full(len(g), i) for i, g in enumerate(groups)])
    obs_H, _ = kruskal(*groups)
    perm_H = []
    for _ in range(n_perm):
        np.random.shuffle(labels)
        perm_groups = [all_vals[labels == i] for i in range(len(groups))]
        if all(len(g) >= 2 for g in perm_groups):
            try:
                H, _ = kruskal(*perm_groups)
                perm_H.append(H)
            except ValueError:
                pass
    perm_H = np.array(perm_H)
    p_perm = np.mean(perm_H >= obs_H) if len(perm_H) > 0 else 1.0
    return obs_H, p_perm

kw_results = []
for ct in cell_types_l1:
    groups = [l1_props[l1_props['plaque_location'] == b][ct].values for b in ORDER]
    groups = [g for g in groups if len(g) >= 3]
    if len(groups) < 2:
        continue
    try:
        H, p = kruskal(*groups)
        eta2 = kruskal_eta2(groups)
        _, p_perm = permutation_kw(groups, n_perm=5000)
        means = [np.mean(g) for g in groups]
        kw_results.append({
            'cell_type': ct, 'level': 'L1',
            'H_stat': H, 'p_value': p, 'p_perm': p_perm,
            'eta_squared': eta2,
            'carotid_mean': means[0] if len(means) > 0 else np.nan,
            'coronary_mean': means[1] if len(means) > 1 else np.nan,
            'femoral_mean': means[2] if len(means) > 2 else np.nan
        })
        sig = '***' if p < 0.001 else '**' if p < 0.01 else '*' if p < 0.05 else 'ns'
        print(f"  {ct}: H={H:.2f}, p={p:.2e} {sig}, eta2={eta2:.4f}")
    except Exception as e:
        print(f"  {ct}: KW failed - {e}")

kw_df = pd.DataFrame(kw_results)
kw_df.to_csv(OUT_DIR / "celltype_proportions_bed_kw.csv", index=False)

# Top discriminating cell types
top_disc = kw_df.nlargest(5, 'eta_squared')['cell_type'].tolist()
print(f"\nTop-5 discriminating L1 types: {top_disc}")

# ============================================================
# STEP 4: Link composition to TI composite
# ============================================================
print("\n=== TI Composite vs Cell Type Proportions ===")

donor_scores = pd.read_csv(RES_DIR / "donor_level_scores.csv", index_col=0)
ti_map = donor_scores['TI_composite']

# Align proportions with TI scores
common_donors = [d for d in l1_props.index if d in ti_map.index]
l1_props_aligned = l1_props.loc[common_donors]
ti_aligned = ti_map.loc[common_donors]

ti_corr_results = []
for ct in cell_types_l1:
    prop = l1_props_aligned[ct].values
    if np.std(prop) == 0:
        continue
    r, p = spearmanr(prop, ti_aligned.values)
    ti_corr_results.append({'cell_type': ct, 'spearman_r': r, 'p_value': p})
    if p < 0.05:
        print(f"  {ct}: r={r:.3f}, p={p:.2e}")

ti_corr_df = pd.DataFrame(ti_corr_results).sort_values('spearman_r', key=abs, ascending=False)
ti_corr_df.to_csv(OUT_DIR / "celltype_ti_correlation.csv", index=False)

# ============================================================
# STEP 5: Module scores for all L1 cell types
# ============================================================
print("\n=== Module Scores per L1 Cell Type ===")

# Define cell-type-specific marker modules
l1_markers = {
    'T_cell_score': ['CD3D', 'CD3E', 'CD3G', 'CD2', 'CD28', 'LCK', 'ZAP70'],
    'CD8_T_score': ['CD8A', 'CD8B', 'GZMK', 'GZMB', 'PRF1', 'NKG7'],
    'CD4_T_score': ['CD4', 'IL7R', 'CCR7', 'SELL', 'TCF7', 'LEF1'],
    'NK_score': ['NKG7', 'GNLY', 'KLRD1', 'KLRF1', 'PRF1', 'GZMB', 'NCAM1'],
    'B_cell_score': ['CD19', 'CD79A', 'CD79B', 'MS4A1', 'PAX5', 'CD22', 'BLNK'],
    'Plasma_score': ['SDC1', 'MZB1', 'SLAMF7', 'XBP1', 'JCHAIN', 'IGHG1'],
    'Myeloid_score': ['CD14', 'CD68', 'FCGR3A', 'CSF1R', 'ITGAM', 'CD33'],
    'DC_score': ['CLEC10A', 'CLEC9A', 'CD1C', 'FCER1A', 'XCR1', 'BATF3', 'IRF8'],
    'Mast_score': ['KIT', 'CPA3', 'TPSAB1', 'FCER1A', 'HDC', 'CMA1'],
    'Neutrophil_score': ['FCGR3B', 'CXCR2', 'CSF3R', 'ELANE', 'MPO', 'MMP8'],
    'EC_score': ['PECAM1', 'CDH5', 'VWF', 'ENG', 'KDR', 'TEK', 'CLDN5'],
    'SMC_score': ['ACTA2', 'MYH11', 'CNN1', 'TAGLN', 'MYOCD', 'LMOD1'],
    'Fibroblast_score': ['COL1A1', 'COL1A2', 'DCN', 'LUM', 'FAP', 'PDGFRA', 'THY1'],
    'Fibromyocyte_score': ['ACTA2', 'COL1A1', 'MYH11', 'TAGLN', 'CNN1', 'DCN'],
}

for score_name, gene_list in l1_markers.items():
    genes_present = [g for g in gene_list if g in adata.var_names]
    if len(genes_present) >= 3:
        sc.tl.score_genes(adata, gene_list=genes_present, score_name=score_name,
                         use_raw=False, ctrl_size=min(len(genes_present), 50))
        print(f"  {score_name}: {len(genes_present)}/{len(gene_list)} genes scored")

# Aggregate to donor-level per L1 cell type
print("\nAggregating to donor-level per cell type...")
score_cols = [s for s in l1_markers.keys() if s in adata.obs.columns]
donor_celltype_scores = adata.obs.groupby(['donor_id', 'cell_type_level1'], observed=False)[score_cols].mean()
donor_celltype_scores = donor_celltype_scores.reset_index()

# Add plaque_location
donor_celltype_scores['plaque_location'] = donor_celltype_scores['donor_id'].map(
    lambda d: donor_meta.loc[d, 'plaque_location'] if d in donor_meta.index else 'unknown')
donor_celltype_scores = donor_celltype_scores[donor_celltype_scores['plaque_location'].isin(ORDER)]

# Save
donor_celltype_scores.to_csv(OUT_DIR / "l1_module_scores_by_bed.csv", index=False)

# Compute mean Z-score per L1 type per bed for heatmap
z_data = []
for ct in cell_types_l1:
    ct_data = donor_celltype_scores[donor_celltype_scores['cell_type_level1'] == ct]
    for bed in ORDER:
        bed_ct = ct_data[ct_data['plaque_location'] == bed]
        if len(bed_ct) >= 3:
            row = {'cell_type': ct, 'vascular_bed': bed}
            for sc_name in score_cols:
                row[sc_name] = bed_ct[sc_name].mean()
            z_data.append(row)

z_df = pd.DataFrame(z_data)
# Z-score per row (gene) for heatmap
z_pivot = z_df.pivot_table(index='cell_type', columns='vascular_bed', values=score_cols, aggfunc='mean')
# For each score, Z-score across cell types
z_heatmap_data = {}
for sc_name in score_cols:
    if sc_name in z_df.columns:
        vals = z_df.groupby(['cell_type', 'vascular_bed'])[sc_name].mean().unstack()
        z_heatmap_data[sc_name] = vals

# ============================================================
# STEP 6: Figure 4
# ============================================================
print("\n=== Generating Figure 4 ===")

plt.rcParams.update({
    'font.family': 'serif',
    'font.serif': ['Times New Roman', 'Times', 'DejaVu Serif', 'serif'],
    'mathtext.fontset': 'custom',
    'mathtext.rm': 'Times New Roman',
    'mathtext.it': 'Times New Roman:italic',
    'mathtext.bf': 'Times New Roman:bold',
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

# L1 color palette (13 types, distinguishable)
L1_COLORS = {
    'T cell': '#1F77B4', 'Macrophage': '#D55E00', 'Smooth Muscle Cell': '#2CA02C',
    'Monocyte': '#FF7F0E', 'NK cell': '#9467BD', 'EC': '#8C564B',
    'B cell': '#E377C2', 'Fibroblast': '#7F7F7F', 'Fibromyocyte': '#BCBD22',
    'Dendritic cell': '#17BECF', 'Neutrophil': '#D62728', 'Mast cell': '#FFBB78',
    'Plasma cell': '#98DF8A'
}

# --- Panel A: Stacked bar of L1 proportions per bed ---
ax1 = fig.add_subplot(2, 3, 1)
bed_means = l1_props.groupby('plaque_location')[cell_types_l1].mean()
bed_sems = l1_props.groupby('plaque_location')[cell_types_l1].sem()

# Sort cell types by overall mean abundance
type_order = bed_means.mean().sort_values(ascending=False).index.tolist()
bed_means = bed_means[type_order]
bed_sems = bed_sems[type_order]

x = np.arange(len(ORDER))
w = 0.6
bottom = np.zeros(len(ORDER))
for ct in type_order:
    vals = bed_means[ct].values
    ax1.bar(x, vals, w, bottom=bottom, label=ct, color=L1_COLORS.get(ct, '#999999'), edgecolor='white', linewidth=0.3)
    bottom += vals
ax1.set_xticks(x)
ax1.set_xticklabels([b.capitalize() for b in ORDER])
ax1.set_ylabel('Proportion', fontsize=10)
ax1.set_title('L1 Cell Type Composition by Vascular Bed', fontsize=10, fontweight='bold')
ax1.legend(fontsize=5.5, loc='upper left', bbox_to_anchor=(1.01, 1), ncol=1)
fig.text(0.01, 0.98, 'A', fontsize=16, fontweight='bold')

# --- Panel B: L2 proportion heatmap ---
ax2 = fig.add_subplot(2, 3, 2)
l2_bed_means = l2_props.groupby('plaque_location')[cell_types_l2].mean()
# Z-score each L2 type across beds
l2_z = l2_bed_means.subtract(l2_bed_means.mean(axis=1), axis=0).divide(
    l2_bed_means.std(axis=1), axis=0).fillna(0)
# Add eta-squared annotations from L1 (parent) as row colors
sns.heatmap(l2_z.T, cmap='RdBu_r', center=0, ax=ax2,
            annot=l2_bed_means.round(3).T, fmt='.3f',
            annot_kws={'fontsize': 5.5},
            cbar_kws={'label': 'Z-score', 'shrink': 0.8},
            linewidths=0.3, linecolor='#EEEEEE')
ax2.set_title('L2 Subtype Proportions\n(Z-score by subtype)', fontsize=10, fontweight='bold')
ax2.tick_params(labelsize=6.5)
fig.text(0.34, 0.98, 'B', fontsize=16, fontweight='bold')

# --- Panel C: Boxplot of top discriminating L1 types ---
ax3 = fig.add_subplot(2, 3, 3)
n_top = min(5, len(top_disc))
plot_data = []
for ct in top_disc[:n_top]:
    for bed in ORDER:
        vals = l1_props[l1_props['plaque_location'] == bed][ct].values
        for v in vals:
            plot_data.append({'Cell Type': ct, 'Vascular Bed': bed.capitalize(), 'Proportion': v})
plot_df = pd.DataFrame(plot_data)
sns.boxplot(data=plot_df, x='Cell Type', y='Proportion', hue='Vascular Bed',
            palette={b.capitalize(): c for b, c in CB_PALETTE.items()}, ax=ax3,
            width=0.6, linewidth=0.8,
            flierprops=dict(marker='o', markersize=3, alpha=0.4))
# Add KW annotations
for i, ct in enumerate(top_disc[:n_top]):
    row = kw_df[kw_df['cell_type'] == ct]
    if len(row) > 0:
        p = row.iloc[0]['p_value']
        sig = '***' if p < 0.001 else '**' if p < 0.01 else '*' if p < 0.05 else 'ns'
        ax3.text(i, plot_df[plot_df['Cell Type'] == ct]['Proportion'].max() * 1.05,
                sig, ha='center', fontsize=9, fontweight='bold')
ax3.set_xticklabels(ax3.get_xticklabels(), rotation=30, ha='right', fontsize=7)
ax3.set_ylabel('Proportion', fontsize=10)
ax3.set_xlabel('')
ax3.legend(fontsize=7, title='', loc='upper right')
ax3.set_title('Top Discriminating L1 Types', fontsize=10, fontweight='bold')
fig.text(0.67, 0.98, 'C', fontsize=16, fontweight='bold')

# --- Panel D: TI composite vs key cell type proportions ---
ax4 = fig.add_subplot(2, 3, 4)
# Use top TI-correlated cell type, or Inflammatory Macrophage proportion if available
ti_target = 'Macrophage'
if ti_target in l1_props_aligned.columns:
    x_vals = l1_props_aligned[ti_target].values
    y_vals = ti_aligned.values
    for bed in ORDER:
        mask = l1_props_aligned['plaque_location'] == bed
        ax4.scatter(x_vals[mask], y_vals[mask], c=CB_PALETTE[bed], alpha=0.6, s=40,
                   edgecolors='white', linewidth=0.5, label=bed)
        if mask.sum() >= 5:
            xb = x_vals[mask]; yb = y_vals[mask]
            lr = stats.linregress(xb, yb)
            x_line = np.linspace(xb.min(), xb.max(), 50)
            ax4.plot(x_line, lr.slope * x_line + lr.intercept, color=CB_PALETTE[bed], linewidth=1.5)
        r, p = spearmanr(x_vals[mask], y_vals[mask])
        sig = '***' if p < 0.001 else '**' if p < 0.01 else '*' if p < 0.05 else ''
        ax4.annotate(f'{bed}: rs={r:.2f}{sig}', xy=(0.05, 0.92 - 0.08 * ORDER.index(bed)),
                    xycoords='axes fraction', fontsize=8, color=CB_PALETTE[bed], fontweight='bold')
    ax4.set_xlabel(f'{ti_target} Proportion', fontsize=10)
    ax4.set_ylabel('TI Composite Score', fontsize=10)
    ax4.set_title('TI Composite vs Macrophage Proportion', fontsize=10, fontweight='bold')
    ax4.legend(fontsize=7)
fig.text(0.01, 0.48, 'D', fontsize=16, fontweight='bold')

# --- Panel E: 3 donut charts per bed ---
def plot_donut(ax, values, labels, colors, title, total_cells):
    """Single donut chart with truncated labels."""
    wedges, texts = ax.pie(values, labels=None, colors=colors, startangle=90,
                          wedgeprops=dict(width=0.35, edgecolor='white', linewidth=0.5))
    ax.text(0, 0, f'n={total_cells:,}', ha='center', va='center', fontsize=8, fontweight='bold')
    ax.set_title(title, fontsize=9, fontweight='bold')

# Use gridspec sub-grid within panel E position (row 2, col 2 = subplot 5)
from matplotlib.gridspec import GridSpecFromSubplotSpec
gs_e = GridSpecFromSubplotSpec(1, 3, subplot_spec=fig.add_subplot(2, 3, 5).get_subplotspec(),
                                wspace=0.3)
for i, bed in enumerate(ORDER):
    ax_e = fig.add_subplot(gs_e[0, i])
    bed_cells = adata[adata.obs['plaque_location'] == bed]
    bed_counts = bed_cells.obs['cell_type_level1'].value_counts()
    threshold = 0.03 * bed_counts.sum()
    main_types = bed_counts[bed_counts >= threshold]
    other_count = bed_counts[bed_counts < threshold].sum()
    if other_count > 0:
        main_types['Other'] = other_count
    main_types = main_types.sort_values(ascending=False)
    colors_e = [L1_COLORS.get(ct, '#CCCCCC') for ct in main_types.index]
    plot_donut(ax_e, main_types.values, main_types.index, colors_e,
              bed.capitalize(), bed_cells.n_obs)
fig.text(0.34, 0.48, 'E', fontsize=16, fontweight='bold')

# --- Panel F: Module score heatmap (L1 types x scores) ---
ax6 = fig.add_subplot(2, 3, 6)
# Compute mean of each score per L1 type (across all beds), Z-score across cell types
score_summary = donor_celltype_scores.groupby('cell_type_level1')[score_cols].mean()
score_z = score_summary.subtract(score_summary.mean(axis=0), axis=1).divide(
    score_summary.std(axis=0), axis=1).fillna(0)
# Clean up labels
score_z.index = [idx.replace('_score', '').replace('_', ' ') for idx in score_z.index]
score_z.columns = [c.replace('_score', '').replace('_', ' ') for c in score_z.columns]
sns.heatmap(score_z, cmap='RdBu_r', center=0, ax=ax6,
            annot=score_summary.round(2).values, fmt='.2f',
            annot_kws={'fontsize': 5.5},
            cbar_kws={'label': 'Z-score', 'shrink': 0.8},
            linewidths=0.3, linecolor='#EEEEEE')
ax6.set_title('Module Scores by Cell Type\n(Z-score across types)', fontsize=10, fontweight='bold')
ax6.tick_params(labelsize=7)
fig.text(0.67, 0.48, 'F', fontsize=16, fontweight='bold')

plt.tight_layout(rect=[0, 0, 1, 0.97])
out_path = OUT_DIR / "figure4_celltype_proportions.png"
fig.savefig(out_path, dpi=300, bbox_inches='tight', facecolor='white', edgecolor='none')
svg_path = out_path.with_suffix('.svg')
fig.savefig(svg_path, bbox_inches='tight', facecolor='white', edgecolor='none')
pdf_path = out_path.with_suffix('.pdf')
fig.savefig(pdf_path, bbox_inches='tight', facecolor='white', edgecolor='none')
plt.close()
print(f"  SVG: {svg_path}")
print(f"  PDF: {pdf_path}")

# Save full donor composition
l1_props.to_csv(OUT_DIR / "donor_full_composition.csv")

# Save top TI-correlated cell types table
ti_corr_df.to_csv(OUT_DIR / "celltype_ti_correlation.csv", index=False)

print("\n=== Phase 4 complete ===")
