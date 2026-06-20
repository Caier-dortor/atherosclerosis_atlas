"""
Phase 7: Pseudotime Trajectory (Macrophage-Myeloid Lineage)
- Myeloid subset: Macrophage L2 + Monocyte + Dendritic cell
- Diffusion Pseudotime (DPT) with CD14 root
- PAGA graph topology
- GAM gene trends along pseudotime
- K-S test by vascular bed
- Regulon-in-pseudotime heatmap (Phase 6 regulons if available)
Uses full atlas for myeloid cells.
"""
import scanpy as sc
import pandas as pd
import numpy as np
from scipy import stats
from scipy.stats import spearmanr, ks_2samp
from pathlib import Path
import matplotlib.pyplot as plt
import json
import warnings
warnings.filterwarnings('ignore')

# === CONFIG ===
RES_DIR = Path("D:/openclaw_workspace/atherosclerosis_atlas/results")
OUT_DIR = RES_DIR / "fig7"
OUT_DIR.mkdir(exist_ok=True)
CB_PALETTE = {'carotid': '#D55E00', 'coronary': '#0072B2', 'femoral': '#009E73'}
ORDER = ['carotid', 'coronary', 'femoral']
np.random.seed(42)

# ============================================================
# STEP 1: Load and subset myeloid cells
# ============================================================
print("Loading full atlas for myeloid subset...")
adata = sc.read_h5ad("D:/openclaw_workspace/atherosclerosis_atlas/data/plaque_atlas.h5ad")
if 'feature_name' in adata.var.columns:
    adata.var_names = adata.var['feature_name']
    adata.var_names = adata.var_names.str.replace(' ', '_')
    adata.var_names_make_unique()

adata.obs['plaque_location'] = adata.obs['origin']

# Subset to myeloid lineage
myeloid_types = ['Macrophage', 'Monocyte', 'Dendritic cell']
myeloid_mask = adata.obs['cell_type_level1'].isin(myeloid_types)
mye = adata[myeloid_mask].copy()
print(f"Myeloid cells: {mye.n_obs:,} (Mac={mye.obs['cell_type_level1'].eq('Macrophage').sum():,}, "
      f"Mono={mye.obs['cell_type_level1'].eq('Monocyte').sum():,}, "
      f"DC={mye.obs['cell_type_level1'].eq('Dendritic cell').sum():,})")

# Create combined L2 labels for PAGA
mye.obs['myeloid_type'] = mye.obs['cell_type_level1'].astype(str)
mac_mask = mye.obs['cell_type_level1'] == 'Macrophage'
mye.obs.loc[mac_mask, 'myeloid_type'] = mye.obs.loc[mac_mask, 'cell_type_level2']
print(f"Myeloid groups: {mye.obs['myeloid_type'].nunique()}")
print(mye.obs['myeloid_type'].value_counts().to_string())

# ============================================================
# STEP 2: Preprocess for trajectory
# ============================================================
print("\n=== Preprocessing for Trajectory ===")

# HVG selection
sc.pp.highly_variable_genes(mye, n_top_genes=3000, flavor='seurat_v3')
mye = mye[:, mye.var['highly_variable']].copy()
print(f"HVG subset: {mye.n_vars} genes")

# PCA
sc.tl.pca(mye, n_comps=50)
print("PCA done")

# Neighbors
sc.pp.neighbors(mye, n_neighbors=15, n_pcs=30)
print("Neighbors done")

# ============================================================
# STEP 3: Diffusion map + DPT
# ============================================================
print("\n=== Diffusion Pseudotime ===")

# Diffusion map
sc.tl.diffmap(mye, n_comps=10)

# Set root: highest CD14 expression
if 'CD14' in mye.var_names:
    cd14_idx = mye.var_names.get_loc('CD14')
    cd14_expr = mye.X[:, cd14_idx].toarray().flatten() if hasattr(mye.X, 'toarray') else mye.X[:, cd14_idx]
    root_idx = np.argmax(cd14_expr)
else:
    # Fallback: use Monocyte with highest Myeloid_score
    root_idx = np.where(mye.obs['cell_type_level1'] == 'Monocyte')[0]
    if len(root_idx) > 0:
        root_idx = root_idx[0]
    else:
        root_idx = 0
print(f"Root cell index: {root_idx}")

# DPT
mye.uns['iroot'] = root_idx
sc.tl.dpt(mye, n_dcs=10)
print(f"DPT range: [{mye.obs['dpt_pseudotime'].min():.2f}, {mye.obs['dpt_pseudotime'].max():.2f}]")

# ============================================================
# STEP 4: PAGA graph
# ============================================================
print("\n=== PAGA ===")
sc.tl.paga(mye, groups='myeloid_type')
print(f"PAGA connectivity: {mye.uns['paga']['connectivities'].shape}")

# ============================================================
# STEP 5: Gene trends along pseudotime
# ============================================================
print("\n=== Gene Trends ===")

# Key TI and macrophage genes
key_genes = [
    'CD14', 'FCGR3A', 'CSF1R', 'ITGAM',  # Myeloid
    'IL1B', 'TNF', 'CCL2', 'CXCL8',       # Inflammatory
    'TREM2', 'APOE', 'FABP5', 'LPL',      # Foamy/Lipid
    'FOLR2', 'LYVE1', 'MRC1',             # Resident
    'HMOX1', 'S100A8', 'S100A9',          # Stress/Acute
    'SPP1', 'FN1',                         # Fibrosis
]

available_key = [g for g in key_genes if g in mye.var_names]
print(f"Key genes available: {len(available_key)}/{len(key_genes)}")

# Bin pseudotime into 50 bins
n_bins = 50
mye.obs['dpt_bin'] = pd.cut(mye.obs['dpt_pseudotime'], bins=n_bins, labels=False)

gene_trends = {}
for gene in available_key:
    gene_idx = mye.var_names.get_loc(gene)
    expr = mye.X[:, gene_idx]
    if hasattr(expr, 'toarray'):
        expr_vals = expr.toarray().flatten()
    else:
        expr_vals = expr.A.flatten()
    # Group by dpt_bin using positional indexing
    bin_means = {}
    for bin_id in range(n_bins):
        bin_mask = mye.obs['dpt_bin'] == bin_id
        if bin_mask.sum() >= 5:
            bin_means[bin_id] = np.mean(expr_vals[bin_mask.values])
        else:
            bin_means[bin_id] = np.nan
    gene_trends[gene] = np.array([bin_means.get(i, np.nan) for i in range(n_bins)])

trend_df = pd.DataFrame(gene_trends, index=range(n_bins))
trend_df.to_csv(OUT_DIR / "gene_pseudotime_trends.csv")

# Per-bed smooth curves
bed_gene_trends = {}
gene_indices = [mye.var_names.get_loc(g) for g in available_key]
for bed in ORDER:
    bed_mye = mye[mye.obs['plaque_location'] == bed]
    if len(bed_mye) < 50:
        continue
    bed_bin_means_arr = []
    for bin_id in range(n_bins):
        bin_mask = bed_mye.obs['dpt_bin'] == bin_id
        if bin_mask.sum() >= 5:
            bin_expr = bed_mye.X[bin_mask.values, :][:, gene_indices]
            if hasattr(bin_expr, 'toarray'):
                bed_bin_means_arr.append(bin_expr.toarray().mean(axis=0))
            else:
                bed_bin_means_arr.append(bin_expr.mean(axis=0))
        else:
            bed_bin_means_arr.append(np.full(len(available_key), np.nan))
    bed_gene_trends[bed] = np.array([v for v in bed_bin_means_arr if not np.all(np.isnan(v))])

# ============================================================
# STEP 6: K-S test by vascular bed
# ============================================================
print("\n=== K-S Test: Pseudotime Distribution by Bed ===")

pt_beds = [mye.obs[mye.obs['plaque_location'] == b]['dpt_pseudotime'].values for b in ORDER]
pt_beds = [p for p in pt_beds if len(p) >= 10]

ks_results = []
for i, b1 in enumerate(ORDER):
    for j, b2 in enumerate(ORDER):
        if i >= j:
            continue
        p1 = mye.obs[mye.obs['plaque_location'] == b1]['dpt_pseudotime'].values
        p2 = mye.obs[mye.obs['plaque_location'] == b2]['dpt_pseudotime'].values
        if len(p1) < 10 or len(p2) < 10:
            continue
        ks_stat, ks_p = ks_2samp(p1, p2)
        mean1, mean2 = np.mean(p1), np.mean(p2)
        ks_results.append({
            'Comparison': f'{b1} vs {b2}',
            'KS_stat': ks_stat, 'p_value': ks_p,
            f'{b1}_mean_pt': mean1, f'{b2}_mean_pt': mean2
        })
        sig = '***' if ks_p < 0.001 else '**' if ks_p < 0.01 else '*' if ks_p < 0.05 else 'ns'
        print(f"  {b1} vs {b2}: KS={ks_stat:.3f}, p={ks_p:.2e} {sig}")

ks_df = pd.DataFrame(ks_results)
ks_df.to_csv(OUT_DIR / "pseudotime_ks_bed_results.csv", index=False)

# ============================================================
# STEP 7: Donor-level pseudotime summary
# ============================================================
print("\n=== Donor-Level Pseudotime Summary ===")

donor_pt = mye.obs.groupby('donor_id')['dpt_pseudotime'].agg(['mean', 'std', 'median']).rename(
    columns={'mean': 'pt_mean', 'std': 'pt_std', 'median': 'pt_median'})
donor_pt['n_cells'] = mye.obs.groupby('donor_id').size()

donor_meta = pd.read_csv(RES_DIR / "donor_metadata.csv", index_col=0)
donor_pt['plaque_location'] = donor_pt.index.map(
    lambda d: donor_meta.loc[d, 'plaque_location'] if d in donor_meta.index else 'unknown')
donor_pt = donor_pt[donor_pt['plaque_location'].isin(ORDER)]
donor_pt.to_csv(OUT_DIR / "pseudotime_donor_summary.csv")

# ============================================================
# STEP 8: Figure 7
# ============================================================
print("\n=== Generating Figure 7 ===")

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

# Myeloid type colors
myeloid_cmap = {
    'Monocyte': '#FF7F0E', 'Dendritic cell': '#17BECF',
    'HMOX1+ Macrophage': '#E41A1C', 'Inflammatory Macrophage': '#FF7F00',
    'Other Macrophage': '#4DAF4A', 'PLIN2+/TREM1+ Macrophage': '#377EB8',
    'TREM2+/Foamy Macrophage': '#984EA3'
}

# UMAP for panels A and later use
if 'X_umap' not in mye.obsm:
    sc.tl.umap(mye)

# --- Panel A: UMAP colored by pseudotime ---
ax1 = fig.add_subplot(2, 3, 1)
pts = mye.obs['dpt_pseudotime'].values
# Subsample for scatter speed
n_dots = min(30000, mye.n_obs)
idx_a = np.random.choice(mye.n_obs, n_dots, replace=False)
sc1 = ax1.scatter(mye.obsm['X_umap'][idx_a, 0], mye.obsm['X_umap'][idx_a, 1],
                 c=pts[idx_a], cmap='viridis', s=0.8, alpha=0.7, rasterized=True)
# Mark root
ax1.scatter(mye.obsm['X_umap'][root_idx, 0], mye.obsm['X_umap'][root_idx, 1],
           c='red', s=80, marker='*', edgecolors='white', linewidth=0.5, label='Root (CD14$^{hi}$)')
ax1.legend(fontsize=6, loc='upper right')
ax1.set_title('Myeloid UMAP colored by\nDiffusion Pseudotime', fontsize=10, fontweight='bold')
ax1.set_xlabel('UMAP 1', fontsize=8)
ax1.set_ylabel('UMAP 2', fontsize=8)
plt.colorbar(sc1, ax=ax1, shrink=0.8, label='DPT')
fig.text(0.01, 0.98, 'A', fontsize=16, fontweight='bold')

# --- Panel B: PAGA graph ---
ax2 = fig.add_subplot(2, 3, 2)
sc.pl.paga(mye, ax=ax2, show=False, node_size_scale=1.8, edge_width_scale=0.8,
           title='PAGA: Myeloid Lineage Topology', fontsize=8,
           node_label_size=6)
ax2.set_title('PAGA: Myeloid Lineage\nTopology', fontsize=10, fontweight='bold')
ax2.text(0.5, -0.12, 'PAGA connectivity suggests monocyte-to-macrophage\nconversion potential (directionality not inferred)',
         transform=ax2.transAxes, fontsize=7, fontstyle='italic', ha='center', va='top', color='grey')
fig.text(0.34, 0.98, 'B', fontsize=16, fontweight='bold')

# --- Panel C: Gene trend plots ---
ax3 = fig.add_subplot(2, 3, 3)
bin_centers = np.linspace(0, 1, n_bins)
genes_to_plot = [g for g in ['TREM1', 'PLIN2', 'CD14', 'CD68', 'HLA-DRA', 'PPARG',
                   'IL1B', 'TNF', 'TREM2', 'APOE', 'HMOX1', 'SPP1', 'FCGR3A'] if g in gene_trends]
if len(genes_to_plot) == 0:
    print('WARNING: No requested genes in gene_trends, using all available')
    genes_to_plot = list(gene_trends.keys())[:8]
n_genes = len(genes_to_plot)
print(f'Panel C: plotting {n_genes} genes along pseudotime')
colors_c = plt.cm.tab10(np.linspace(0, 1, n_genes))

# Smooth trends
bin_valid = np.arange(n_bins)
for k, gene in enumerate(genes_to_plot):
    vals = gene_trends[gene]
    valid = ~np.isnan(vals)
    if valid.sum() < 5:
        continue
    # Lowess smooth
    from scipy.interpolate import make_interp_spline
    x_valid = bin_valid[valid] / n_bins
    y_valid = vals[valid]
    # Z-score normalize
    y_z = (y_valid - np.nanmean(y_valid)) / (np.nanstd(y_valid) + 1e-10)
    ax3.plot(x_valid, y_z, color=colors_c[k], alpha=0.4, linewidth=0.5)
    # Smooth with BSpline
    try:
        x_smooth = np.linspace(0, 1, 200)
        if len(x_valid) >= 4:
            spl = make_interp_spline(x_valid, y_z, k=min(3, len(x_valid) - 1))
            ax3.plot(x_smooth, spl(x_smooth), color=colors_c[k], linewidth=1.5, label=gene)
    except Exception:
        pass

ax3.legend(fontsize=5.5, ncol=2, loc='upper right')
ax3.set_xlabel('Pseudotime (normalized)', fontsize=9)
ax3.set_ylabel('Z-scored expression', fontsize=9)
ax3.set_title('Key Gene Expression\nalong Pseudotime', fontsize=10, fontweight='bold')
fig.text(0.67, 0.98, 'C', fontsize=16, fontweight='bold')

# --- Panel D: Pseudotime distribution by bed ---
ax4 = fig.add_subplot(2, 3, 4)
bed_pt_data = [mye.obs[mye.obs['plaque_location'] == b]['dpt_pseudotime'].values for b in ORDER]

# Violin-like KDE
for i, (bed, pts_bed) in enumerate(zip(ORDER, bed_pt_data)):
    if len(pts_bed) < 10:
        continue
    from scipy.stats import gaussian_kde
    try:
        kde = gaussian_kde(pts_bed)
        x_kde = np.linspace(pts_bed.min(), pts_bed.max(), 200)
        y_kde = kde(x_kde)
        y_kde = y_kde / y_kde.max() * 0.4
        offset = i + 1
        ax4.fill_betweenx(x_kde, offset - y_kde, offset + y_kde, alpha=0.4, color=CB_PALETTE[bed])
        ax4.fill_betweenx(x_kde, offset - y_kde, offset + y_kde, alpha=0.7, color=CB_PALETTE[bed], linewidth=0.8)
        # Mean line
        ax4.axhline(np.median(pts_bed), xmin=(offset - 0.35) / (len(ORDER) + 1),
                   xmax=(offset + 0.35) / (len(ORDER) + 1), color=CB_PALETTE[bed], linewidth=2)
    except Exception:
        pass

# K-S annotations
for i, b1 in enumerate(ORDER):
    for j, b2 in enumerate(ORDER):
        if i >= j:
            continue
        row_match = [r for r in ks_results if r['Comparison'].startswith(b1) and b2 in r['Comparison']]
        if row_match:
            p = row_match[0]['p_value']
            sig = '***' if p < 0.001 else '**' if p < 0.01 else '*' if p < 0.05 else 'ns'
            y_pos = mye.obs['dpt_pseudotime'].max() * 1.05 + (ORDER.index(b2)) * 0.03
            ax4.text((i + j + 2) / 2, y_pos, sig, ha='center', fontsize=8, fontweight='bold')

ax4.set_xlim(0.5, 3.5)
ax4.set_xticks([1, 2, 3])
ax4.set_xticklabels([b.capitalize() for b in ORDER])
ax4.set_ylabel('Diffusion Pseudotime', fontsize=9)
ax4.set_title('Pseudotime Distribution\nby Vascular Bed', fontsize=10, fontweight='bold')
fig.text(0.01, 0.48, 'D', fontsize=16, fontweight='bold')

# --- Panel E: UMAP by bed (3 sub-panels within Panel E position) ---
from matplotlib.gridspec import GridSpecFromSubplotSpec
gs_e = GridSpecFromSubplotSpec(1, 3, subplot_spec=fig.add_subplot(2, 3, 5).get_subplotspec(),
                                wspace=0.2)
for i, bed in enumerate(ORDER):
    ax_e = fig.add_subplot(gs_e[0, i])
    bed_idx = mye.obs['plaque_location'] == bed
    bed_cells = mye[bed_idx]
    if len(bed_cells) < 10:
        ax_e.text(0.5, 0.5, f'{bed}\nn<10', transform=ax_e.transAxes, ha='center', fontsize=6)
        ax_e.axis('off')
        continue
    n_d = min(5000, bed_cells.n_obs)
    idx_d = np.random.choice(bed_cells.n_obs, n_d, replace=False)
    ax_e.scatter(bed_cells.obsm['X_umap'][idx_d, 0], bed_cells.obsm['X_umap'][idx_d, 1],
                c=bed_cells.obs['dpt_pseudotime'].values[idx_d], cmap='viridis', s=0.8, alpha=0.6, rasterized=True)
    ax_e.set_title(bed.capitalize(), fontsize=8, fontweight='bold', color=CB_PALETTE[bed])
    ax_e.set_xticks([])
    ax_e.set_yticks([])

# Panel E title
fig.text(0.34, 0.48, 'E', fontsize=16, fontweight='bold')

# --- Panel F: Regulon-in-pseudotime heatmap ---
ax6 = fig.add_subplot(2, 3, 6)
# Score regulons in myeloid cells using Phase 6 targets
import json
regulon_json_path = RES_DIR / "fig6" / "regulon_target_genes.json"
if regulon_json_path.exists():
    print("Scoring Phase 6 regulons in myeloid cells...")
    with open(regulon_json_path, 'r') as f:
        regulon_target_dict = json.load(f)

    regulon_cols_6 = []
    for tf, targets in regulon_target_dict.items():
        score_name = f'{tf}_regulon'
        if score_name in mye.obs.columns:
            regulon_cols_6.append(score_name)
            continue
        genes_present = [g for g in targets['positive'][:50] if g in mye.var_names]
        if len(genes_present) >= 5:
            sc.tl.score_genes(mye, gene_list=genes_present, score_name=score_name,
                            use_raw=False, ctrl_size=min(len(genes_present), 50))
            regulon_cols_6.append(score_name)
    print(f"  Scored {len(regulon_cols_6)} regulons in myeloid cells")
else:
    # Fallback: score TF expression directly
    print("No Phase 6 regulons found, using TF expression directly...")
    fallback_tfs = ['SPI1', 'CEBPB', 'IRF8', 'STAT1', 'NFKB1', 'JUN', 'FOS', 'HIF1A', 'KLF4',
                    'STAT4', 'MYC', 'PPARG', 'EGR1', 'FOSB', 'EPAS1']
    regulon_cols_6 = []
    for tf in fallback_tfs:
        if tf not in mye.var_names:
            continue
        tf_idx = mye.var_names.get_loc(tf)
        score_name = f'{tf}_regulon'
        expr = mye.X[:, tf_idx]
        mye.obs[score_name] = expr.toarray().flatten() if hasattr(expr, 'toarray') else expr
        regulon_cols_6.append(score_name)
    print(f"  Fallback regulons: {len(regulon_cols_6)}")

# Bin regulon activity by pseudotime
if regulon_cols_6:
    regulon_bins = mye.obs.groupby('dpt_bin')[regulon_cols_6].mean()
    regulon_bins_z = regulon_bins.subtract(regulon_bins.mean(axis=0), axis=1).divide(
        regulon_bins.std(axis=0), axis=1).fillna(0)

    # Keep top-20 with highest variance
    top_reg_vars = regulon_bins_z.var().nlargest(min(20, len(regulon_cols_6))).index.tolist()
    im6 = ax6.imshow(regulon_bins_z[top_reg_vars].T.values, aspect='auto', cmap='RdBu_r',
                      vmin=-2, vmax=2, interpolation='nearest')
    ax6.set_yticks(range(len(top_reg_vars)))
    ax6.set_yticklabels([r.replace('_regulon', '') for r in top_reg_vars], fontsize=5.5)
    ax6.set_xticks(np.linspace(0, n_bins - 1, 10).astype(int))
    ax6.set_xticklabels([f'{i/100:.2f}' for i in np.linspace(0, 100, 10).astype(int)])
    ax6.set_xlabel('Pseudotime (normalized)', fontsize=9)
    ax6.set_title('Regulon Activity\nalong Pseudotime', fontsize=10, fontweight='bold')
    plt.colorbar(im6, ax=ax6, shrink=0.8, label='Z-score')
else:
    ax6.text(0.5, 0.5, 'No regulon data available', ha='center', va='center', transform=ax6.transAxes)
    ax6.set_title('Regulon Activity\n(awaiting Phase 6)', fontsize=10, fontweight='bold')

fig.text(0.67, 0.48, 'F', fontsize=16, fontweight='bold')

plt.tight_layout(rect=[0, 0, 1, 0.97])
out_path = OUT_DIR / "figure7_pseudotime_trajectory.png"
fig.savefig(out_path, dpi=300, bbox_inches='tight', facecolor='white', edgecolor='none')
svg_path = out_path.with_suffix('.svg')
fig.savefig(svg_path, bbox_inches='tight', facecolor='white', edgecolor='none')
pdf_path = out_path.with_suffix('.pdf')
fig.savefig(pdf_path, bbox_inches='tight', facecolor='white', edgecolor='none')
plt.close()
print(f"  Saved: {out_path}")
print(f"  SVG: {svg_path}")
print(f"  PDF: {pdf_path}")

print("\n=== Phase 7 complete ===")
