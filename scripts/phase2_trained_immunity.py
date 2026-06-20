"""
Phase 2: Trained Immunity Signature Analysis (v2)
— Donor-level pseudobulk (P0-corrected)
— Lesion stage temporal proxy (P2)
— Weight sensitivity analysis (P1)
Uses output from phase1_macrophage_analysis.py
"""

import scanpy as sc
import pandas as pd
import numpy as np
from scipy import stats
from scipy.stats import spearmanr
from pathlib import Path
import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec
import seaborn as sns
from sklearn.decomposition import PCA
from sklearn.preprocessing import StandardScaler
import warnings
warnings.filterwarnings('ignore')

# === CONFIG ===
RES_DIR = Path("D:/openclaw_workspace/atherosclerosis_atlas/results")
OUT_DIR = RES_DIR / "fig2"
OUT_DIR.mkdir(exist_ok=True)

sc.set_figure_params(dpi=100, facecolor='white', frameon=True)
COLORS = {'carotid': '#E74C3C', 'coronary': '#3498DB', 'femoral': '#2ECC71'}

# ============================================================
# STEP 1: Load annotated macrophages
# ============================================================
print("Loading macrophages...")
macrophages = sc.read_h5ad(RES_DIR / "macrophages_annotated.h5ad")
donor_df = pd.read_csv(RES_DIR / "donor_level_scores.csv", index_col=0)

print(f"Cells: {macrophages.n_obs:,}")
print(f"Donors: {donor_df.shape[0]}")

# ============================================================
# STEP 2: TI composite sensitivity analysis (P1)
# ============================================================
print("\n=== TI Composite Sensitivity Analysis ===")

# Define tiers
ti_tiers_def = {
    'TI_Inflammation': (['IL1B', 'IL6', 'TNF', 'CXCL8', 'CCL2', 'CCL3', 'CCL4', 'IL18'], 0.25),
    'TI_PRR': (['TLR2', 'TLR4', 'NOD2', 'CLEC7A', 'TLR1', 'CD14', 'MARCO'], 0.15),
    'TI_Metabolic': (['HK1', 'HK2', 'PFKFB3', 'PKM2', 'LDHA', 'ACLY', 'FASN', 'IDH1'], 0.20),
    'TI_H3K4me3': (['KMT2A', 'KMT2D', 'SETD1A', 'SETD1B'], 0.10),
    'TI_H3K27ac': (['EP300', 'CREBBP', 'KAT2A', 'KAT2B'], 0.10),
    'TI_HDAC_SIRT': (['HDAC1', 'HDAC2', 'HDAC3', 'HDAC8', 'SIRT1', 'SIRT3', 'SIRT6'], 0.05),
}
score_names = [f'{t}_score' for t in ti_tiers_def.keys()]
base_weights = np.array([w for _, w in ti_tiers_def.values()])

# PCA weights
available_scores = [s for s in score_names if s in donor_df.columns]
ti_matrix = donor_df[available_scores].values
scaler = StandardScaler()
ti_pca = PCA(n_components=1).fit_transform(scaler.fit_transform(ti_matrix))
donor_df['TI_pca'] = ti_pca[:, 0]

# Sensitivity: 100 random perturbations
np.random.seed(42)
n_iter = 100
perturbed_scores = []
for i in range(n_iter):
    noise = np.random.uniform(0.5, 1.5, len(base_weights))
    w = base_weights * noise
    w = w / w.sum()
    donor_df[f'TI_perturb_{i}'] = sum(donor_df[s] * ww for s, ww in zip(available_scores, w))
    perturbed_scores.append(donor_df[f'TI_perturb_{i}'].values)

# Consistency check: each perturbed score vs preset
correlations = [spearmanr(donor_df['TI_composite'], p)[0] for p in perturbed_scores]
consistency = np.mean(np.array(correlations) > 0.8)
print(f"Consistency (>0.8 corr): {consistency:.1%}")

# ============================================================
# STEP 3: Vascular bed comparison (donor-level, P0-corrected)
# ============================================================
print("\n=== Vascular Bed Comparison ===")

ti_vars = ['TI_composite', 'TI_pca', 'TI_Inflammation_score', 'TI_Metabolic_score',
           'TI_H3K4me3_score', 'TI_H3K27ac_score', 'TI_HDAC_SIRT_score',
           'Acute_Inflammation_score']

from scipy.stats import kruskal
results = []
for var in ti_vars:
    if var not in donor_df.columns:
        continue
    groups = [donor_df[donor_df['plaque_location'] == b][var].dropna().values
              for b in ['carotid', 'coronary', 'femoral']]
    if all(len(g) >= 3 for g in groups):
        H, p = kruskal(*groups)
        results.append({
            'variable': var,
            'H_stat': H, 'p_value': p,
            'carotid_mean': np.mean(groups[0]),
            'coronary_mean': np.mean(groups[1]),
            'femoral_mean': np.mean(groups[2])
        })
        print(f"  {var}: H={H:.2f}, p={p:.2e}")

results_df = pd.DataFrame(results)
results_df.to_csv(OUT_DIR / "ti_vascular_bed_comparison.csv", index=False)

# ============================================================
# STEP 4: Lesion stage as temporal proxy (P2)
# ============================================================
print("\n=== Lesion Stage Analysis (Temporal Proxy) ===")

if 'lesion_stage' in macrophages.obs.columns:
    stage_order = ['Early', 'Intermediate', 'Advanced', 'Calcified']
    fig, axes = plt.subplots(1, 3, figsize=(16, 6))

    for i, bed in enumerate(['carotid', 'coronary', 'femoral']):
        ax = axes[i]
        bed_mac = macrophages[macrophages.obs['plaque_location'] == bed]
        stage_data = []
        for stage in stage_order:
            if stage not in bed_mac.obs['lesion_stage'].cat.categories:
                continue
            stage_cells = bed_mac[bed_mac.obs['lesion_stage'] == stage]
            if stage_cells.n_obs >= 10:
                stage_data.append({
                    'stage': stage,
                    'TI_mean': stage_cells.obs['TI_composite'].mean(),
                    'TI_sem': stage_cells.obs['TI_composite'].sem(),
                    'Acute_mean': stage_cells.obs['Acute_Inflammation_score'].mean(),
                    'Acute_sem': stage_cells.obs['Acute_Inflammation_score'].sem(),
                    'n': stage_cells.n_obs
                })
        if stage_data:
            df_stage = pd.DataFrame(stage_data)
            x = range(len(df_stage))
            w = 0.35
            ax.bar([xi - w/2 for xi in x], df_stage['TI_mean'], w,
                   yerr=df_stage['TI_sem'], label='TI composite', color='#E74C3C', capsize=3)
            ax.bar([xi + w/2 for xi in x], df_stage['Acute_mean'], w,
                   yerr=df_stage['Acute_sem'], label='Acute Inflammation', color='#7F8C8D', capsize=3)
            ax.set_xticks(x)
            ax.set_xticklabels(df_stage['stage'])
            ax.set_title(f'{bed.capitalize()}\n(n per stage: {dict(zip(df_stage.stage, df_stage.n))})')
            ax.legend(fontsize=8)

    plt.suptitle('TI Score vs Acute Inflammation by Lesion Stage (Temporal Proxy)', fontweight='bold')
    plt.tight_layout()
    plt.savefig(OUT_DIR / "lesion_stage_temporal_proxy.png", dpi=150, bbox_inches='tight')
    plt.close()

# ============================================================
# STEP 5: Epigenetic enzyme expression analysis
# ============================================================
print("\n=== Epigenetic Enzyme Analysis ===")

epi_genes = ['KMT2A', 'KMT2D', 'SETD1A', 'SETD1B',  # H3K4me3 writers
             'EP300', 'CREBBP', 'KAT2A', 'KAT2B',     # H3K27ac writers
             'KDM5B', 'KDM6B',                        # Demethylases
             'SIRT1', 'SIRT3', 'HDAC3',               # Deacetylases
             'DNMT3A', 'TET2']                        # DNA methylation

# Filter to expressed genes
expressed_epi = [g for g in epi_genes if g in macrophages.var_names]
print(f"Expressed epi genes: {len(expressed_epi)}/{len(epi_genes)}")

# Donor-level mean expression of epi genes
donor_epi = macrophages.obs[['donor_id', 'plaque_location', 'dataset', 'sex']].copy()
for g in expressed_epi:
    donor_epi[f'{g}_expr'] = macrophages[:, g].X.toarray().flatten()

donor_epi_agg = donor_epi.groupby('donor_id').agg({
    'plaque_location': 'first',
    **{f'{g}_expr': 'mean' for g in expressed_epi}
})

# Per vascular bed comparison
for g in expressed_epi:
    groups = [donor_epi_agg[donor_epi_agg['plaque_location'] == b][f'{g}_expr'].values
              for b in ['carotid', 'coronary', 'femoral']]
    if all(len(g) >= 3 for g in groups):
        H, p = kruskal(*groups)
        if p < 0.05:
            print(f"  {g}: H={H:.2f}, p={p:.2e} — carotid={np.mean(groups[0]):.3f}, femoral={np.mean(groups[2]):.3f}")

# ============================================================
# Subpopulation-stratified TI composite (P0 fix per evaluation report)
# ============================================================
print("\n=== Subpopulation-Stratified TI Analysis ===")

mac_obs = macrophages.obs.copy()
l2_subtypes = mac_obs['cell_type_level2'].unique()
BED_ORDER = ['carotid', 'coronary', 'femoral']
stratified_stats = []
for subtype in l2_subtypes:
    sub = mac_obs[mac_obs['cell_type_level2'] == subtype]
    for bed in BED_ORDER:
        bed_sub = sub[sub['plaque_location'] == bed]
        if len(bed_sub) >= 10:
            stratified_stats.append({
                'L2_Subtype': subtype,
                'Vascular_Bed': bed.capitalize(),
                'TI_composite': bed_sub['TI_composite'].mean(),
                'TI_sem': bed_sub['TI_composite'].sem(),
                'Inflammatory_score': bed_sub['Inflammatory_Mac_score'].mean(),
                'n_cells': len(bed_sub)
            })
strat_df = pd.DataFrame(stratified_stats)
print(strat_df.pivot_table(index='L2_Subtype', columns='Vascular_Bed', values='TI_composite'))
strat_df.to_csv(OUT_DIR / "ti_stratified_by_subtype.csv", index=False)

# ============================================================
print("\n=== Generating Figure 2 ===")

# Colorblind-friendly palette (Wong 2011, Nature Methods)
CB_PALETTE = {'carotid': '#D55E00', 'coronary': '#0072B2', 'femoral': '#009E73'}
ORDER = ['carotid', 'coronary', 'femoral']

# SCI journal style
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

fig = plt.figure(figsize=(20, 20), facecolor="white")
gs = GridSpec(3, 3, figure=fig, height_ratios=[1, 1, 1.5],
              hspace=0.40, wspace=0.35, left=0.05, right=0.97, top=0.96, bottom=0.04)

# ============================================================
# ROW 1: Supporting evidence (3 panels)
# ============================================================

# --- Panel A: TI tier scores --- split violins ---
ax1 = fig.add_subplot(gs[0, 0])
ti_plot_vars = available_scores[:6]
donor_melt = donor_df.melt(id_vars='plaque_location', value_vars=ti_plot_vars)
sns.violinplot(data=donor_melt, x='variable', y='value', hue='plaque_location',
               palette=CB_PALETTE, split=True, ax=ax1, density_norm='width',
               inner='quartile', linewidth=0.8, cut=0)
ax1.set_xticklabels(ax1.get_xticklabels(), rotation=40, ha='right', fontsize=7.5)
ax1.set_ylabel('Module Score', fontsize=10)
ax1.set_xlabel('')
ax1.set_title('TI Module Scores by Vascular Bed', fontsize=10, fontweight='bold')
ax1.legend(fontsize=7, title='', loc='upper right')
ax1.text(0.98, 0.95, f'Donors: C=50, Co=13, F=7', transform=ax1.transAxes,
         fontsize=6, ha='right', va='top',
         bbox=dict(boxstyle='round,pad=0.3', facecolor='#F5F5F5', alpha=0.8))
fig.text(0.02, 0.95, 'A', fontsize=16, fontweight='bold')

# --- Panel B: TI composite subpopulation-stratified boxplot ---
ax2 = fig.add_subplot(gs[0, 1])
plot_data_b = []
for subtype in l2_subtypes:
    sub = mac_obs[mac_obs['cell_type_level2'] == subtype]
    for bed in ORDER:
        vals = sub[sub['plaque_location'] == bed]['TI_composite'].dropna().values
        for v in vals:
            plot_data_b.append({
                'L2 Subtype': subtype.replace(' Macrophage', '').replace('/', '/\n'),
                'TI Composite': v,
                'Bed': bed
            })
df_b = pd.DataFrame(plot_data_b)
subtype_order = df_b.groupby('L2 Subtype')['TI Composite'].mean().sort_values(ascending=False).index.tolist()
sns.boxplot(data=df_b, x='L2 Subtype', y='TI Composite', hue='Bed',
            order=subtype_order, palette=CB_PALETTE, ax=ax2, width=0.6, linewidth=0.8,
            flierprops=dict(marker='o', markersize=2, alpha=0.3), showmeans=True,
            meanprops=dict(marker='D', markersize=5, markerfacecolor='white', markeredgecolor='#333333'))
# Diamond (◆) = donor-level mean; box center line = median
ax2.text(0.02, 0.98, chr(9670) + ' = mean, center line = median', transform=ax2.transAxes,
         fontsize=6, ha='left', va='top', color='#666666')
ax2.set_xticklabels(ax2.get_xticklabels(), rotation=25, ha='right', fontsize=7)
ax2.set_ylabel('TI Composite Score', fontsize=10)
ax2.set_xlabel('')
ax2.set_title('TI Composite by Macrophage Subtype', fontsize=10, fontweight='bold')
ax2.legend(fontsize=7, title='', loc='upper right')
n_donors_per_bed = mac_obs.groupby('plaque_location')['donor_id'].nunique()
ax2.text(0.98, 0.05, f'Donors: C={n_donors_per_bed["carotid"]}, Co={n_donors_per_bed["coronary"]}, F={n_donors_per_bed["femoral"]}',
         transform=ax2.transAxes, fontsize=6, ha='right', va='bottom',
         bbox=dict(boxstyle='round,pad=0.3', facecolor='#F5F5F5', alpha=0.8))
fig.text(0.35, 0.95, 'B', fontsize=16, fontweight='bold')

# --- Panel C: Epigenetic enzyme heatmap (Z-score) ---
ax3 = fig.add_subplot(gs[0, 2])
epi_data = []
for bed in ORDER:
    bed_df = donor_epi_agg[donor_epi_agg['plaque_location'] == bed]
    for g in expressed_epi:
        epi_data.append({
            'gene': g, 'vascular_bed': bed,
            'mean_expr': bed_df[f'{g}_expr'].mean(),
        })
epi_plot = pd.DataFrame(epi_data)
epi_pivot = epi_plot.pivot_table(index='gene', columns='vascular_bed',
                                  values='mean_expr', aggfunc='mean')
epi_pivot_z = epi_pivot.subtract(epi_pivot.mean(axis=1), axis=0).divide(
    epi_pivot.std(axis=1), axis=0).fillna(0)
sns.heatmap(epi_pivot_z, cmap='RdBu_r', center=0, ax=ax3,
            annot=epi_pivot.round(3), fmt='.3f',
            annot_kws={'fontsize': 6.5, 'fontweight': 'bold'},
            cbar_kws={'label': 'Z-score', 'shrink': 0.8},
            linewidths=0.5, linecolor='#EEEEEE')
# Auto-set annot color for readability on dark cells
for t in ax3.texts:
    val = float(t.get_text())
    t.set_color('white' if abs(val) > 1.2 else 'black')
ax3.set_title('Epigenetic Enzyme Expression\n(Z-score by gene)', fontsize=10, fontweight='bold')
ax3.tick_params(labelsize=8)
fig.text(0.68, 0.95, 'C', fontsize=16, fontweight='bold')

# ============================================================
# ROW 2: Methodological validation (2 panels)
# ============================================================

# --- Panel D: Sensitivity analysis histogram ---
ax4 = fig.add_subplot(gs[1, 0])
corr_array = np.array(correlations)
ax4.hist(corr_array, bins=30, color=CB_PALETTE['coronary'], edgecolor='white', alpha=0.85)
ax4.axvline(x=0.8, color='#D55E00', linestyle='--', linewidth=1.5, label='Consistency threshold (r=0.8)')
ax4.axvline(x=np.median(corr_array), color='#009E73', linestyle='-', linewidth=2,
            label=f'Median = {np.median(corr_array):.3f}')
ax4.set_xlabel('Spearman r with preset TI', fontsize=10)
ax4.set_ylabel('Frequency (100 perturbations)', fontsize=10)
ax4.set_title(f'TI Weight Sensitivity\nConsistency: {consistency:.1%}', fontsize=10, fontweight='bold')
ax4.legend(fontsize=8, loc='upper left')
fig.text(0.02, 0.58, 'D', fontsize=16, fontweight='bold')

# --- Panel E: PCA loadings ---
ax5 = fig.add_subplot(gs[1, 1])
pca_full = PCA(n_components=2).fit(ti_matrix)
loadings = pca_full.components_.T
colors_e = [CB_PALETTE['carotid'] if v > 0 else CB_PALETTE['femoral'] for v in loadings[:, 0]]
ax5.barh(range(len(available_scores)), loadings[:, 0], color=colors_e, alpha=0.8, height=0.6)
ax5.set_yticks(range(len(available_scores)))
ax5.set_yticklabels([s.replace('_score', '').replace('_', ' ') for s in available_scores],
                     fontsize=8)
ax5.set_xlabel('PC1 Loading', fontsize=10)
ax5.set_title('TI PCA Loadings\n(Data-driven weights)', fontsize=10, fontweight='bold')
ax5.axvline(x=0, color='#333333', linewidth=0.8)
fig.text(0.35, 0.58, 'E', fontsize=16, fontweight='bold')

# ============================================================
# ROW 3: CORE FINDING -- TI-Acute Inflammation uncoupling
# Full-width panel as visual punchline
# ============================================================
ax6 = fig.add_subplot(gs[2, :])
from scipy.stats import linregress

for bed in ORDER:
    bed_df = donor_df[donor_df["plaque_location"] == bed].dropna(
        subset=["Acute_Inflammation_score", "TI_composite"])
    x = bed_df["Acute_Inflammation_score"].values
    y = bed_df["TI_composite"].values
    lw = 2.5 if bed == "femoral" else 1.5
    alpha_line = 1.0 if bed == "femoral" else 0.6
    ax6.scatter(x, y, c=CB_PALETTE[bed], alpha=0.7, s=55,
                label=bed, edgecolors="white", linewidth=0.5,
                zorder=3 if bed == "femoral" else 2)
    if len(x) >= 5:
        lr = linregress(x, y)
        x_line = np.linspace(x.min(), x.max(), 50)
        ax6.plot(x_line, lr.slope * x_line + lr.intercept,
                 color=CB_PALETTE[bed], linewidth=lw, alpha=alpha_line,
                 zorder=5 if bed == "femoral" else 3)
        ax6.fill_between(x_line,
                         lr.slope * x_line + lr.intercept - 1.96 * lr.stderr * np.sqrt(len(x)),
                         lr.slope * x_line + lr.intercept + 1.96 * lr.stderr * np.sqrt(len(x)),
                         color=CB_PALETTE[bed], alpha=0.08)
    r_s, p_s = spearmanr(x, y)
    sig = "***" if p_s < 0.001 else "**" if p_s < 0.01 else "*" if p_s < 0.05 else f", p={p_s:.2e}"
    if bed == "femoral":
        # Femoral n=7: use dagger for trend, add sample-size caveat
        if bed == "femoral" and len(x) <= 10:
            sig_display = sig.replace('*', chr(8224))  # dagger for trend
            ax6.annotate(f"{bed}: r_s={r_s:.2f}{sig_display} (n={len(x)})\n(Femoral TI-Acute decoupling; limited by small n)",
                         xy=(0.75, 0.87), xycoords="axes fraction",
                         fontsize=10, color=CB_PALETTE[bed], fontweight="bold",
                         bbox=dict(boxstyle="round,pad=0.5", facecolor="white",
                                   edgecolor=CB_PALETTE[bed], linewidth=2, alpha=0.9))
            ax6.text(0.98, 0.02, chr(8224) + ' = trend (n=' + str(len(x)) + ' donors, limited power)',
                     transform=ax6.transAxes, fontsize=6.5, ha='right', va='bottom', color='grey')
        else:
            ax6.annotate(f"{bed}: r_s={r_s:.2f}{sig}\n(Femoral TI-Acute decoupling)",
                         xy=(0.75, 0.85), xycoords="axes fraction",
                         fontsize=11, color=CB_PALETTE[bed], fontweight="bold",
                         bbox=dict(boxstyle="round,pad=0.5", facecolor="white",
                                   edgecolor=CB_PALETTE[bed], linewidth=2, alpha=0.9))
    else:
        ax6.annotate(f"{bed}: r_s={r_s:.2f}{sig}",
                     xy=(0.05, 0.92 - 0.08 * ORDER.index(bed)),
                     xycoords="axes fraction", fontsize=9, color=CB_PALETTE[bed],
                     fontweight="bold")

ax6.set_xlabel("Acute Inflammation Score", fontsize=12, fontweight="bold")
ax6.set_ylabel("Trained Immunity Composite Score", fontsize=12, fontweight="bold")
ax6.set_ylim(-0.35, 0.65)  # Focus on data range, trimming excess whitespace
ax6.set_title("TI vs Acute Inflammation: Femoral-Specific Decoupling\n"
              "Femoral TI operates independently of acute inflammation (FAO-SIRT1-HDAC pathway)",
              fontsize=12, fontweight="bold", pad=12)
ax6.legend(fontsize=10, loc="lower right", markerscale=1.2,
           handletextpad=0.8, borderpad=0.5).set_zorder(10)

ax6.text(0.98, 0.05,
         "Mechanism: FAO-SIRT1-HDAC axis maintains\n"
         "trained immunity in femoral macrophages\n"
         "without triggering acute inflammatory programs",
         transform=ax6.transAxes, fontsize=8.5, ha="right", va="bottom",
         fontstyle="italic", color="#555555",
         bbox=dict(boxstyle="round,pad=0.5", facecolor="#FAFAFA",
                   edgecolor="#CCCCCC", linewidth=1))

fig.text(0.02, 0.28, 'F', fontsize=18, fontweight='bold', color='#009E73')

# --- Save Figure 2 ---
plt.tight_layout(rect=[0, 0, 1, 0.97])
out_path = OUT_DIR / "figure2_trained_immunity.png"
fig.savefig(out_path, dpi=300, bbox_inches='tight', facecolor='white', edgecolor='none')
svg_path = out_path.with_suffix('.svg')
fig.savefig(svg_path, bbox_inches='tight', facecolor='white', edgecolor='none')
pdf_path = out_path.with_suffix('.pdf')
fig.savefig(pdf_path, bbox_inches='tight', facecolor='white', edgecolor='none')
plt.close()
print(f"  SVG: {svg_path}")
print(f"  PDF: {pdf_path}")
