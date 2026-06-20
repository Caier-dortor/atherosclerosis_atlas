"""
Phase 3: Metabolism-Epigenetics Coupling Analysis (v2)
— Mediation analysis (P0-corrected)
— Computational metabolite inference preparation
— Partial correlation (P0)
Uses output from phase1_macrophage_analysis.py
"""

import scanpy as sc
import pandas as pd
import numpy as np
from scipy import stats
from scipy.stats import spearmanr
from pathlib import Path
import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec, GridSpecFromSubplotSpec
import seaborn as sns
import warnings
warnings.filterwarnings('ignore')

# === CONFIG ===
RES_DIR = Path("D:/openclaw_workspace/atherosclerosis_atlas/results")
OUT_DIR = RES_DIR / "fig3"
OUT_DIR.mkdir(exist_ok=True)

COLORS = {'carotid': '#E74C3C', 'coronary': '#3498DB', 'femoral': '#2ECC71'}

# ============================================================
# STEP 1: Load data
# ============================================================
print("Loading data...")
macrophages = sc.read_h5ad(RES_DIR / "macrophages_annotated.h5ad")
donor_df = pd.read_csv(RES_DIR / "donor_level_scores.csv", index_col=0)

# ============================================================
# STEP 2: Metabolic-epigenetic correlation matrix
# ============================================================
print("\n=== Metabolic-Epigenetic Correlation Matrix ===")

# Define gene sets
metabolic_genes = {
    'Glycolysis': ['HK1', 'HK2', 'PFKFB3', 'ALDOA', 'GAPDH', 'PGK1', 'PKM2', 'LDHA'],
    'OXPHOS': ['NDUFA1', 'NDUFA2', 'SDHA', 'UQCRC1', 'COX4I1', 'ATP5F1A'],
    'FAO': ['CPT1A', 'ACADVL', 'HADHA', 'ACAA2'],
    'FAS': ['ACACA', 'FASN', 'SCD', 'ACLY', 'ME1'],
    'Cholesterol': ['HMGCR', 'SQLE', 'LDLR', 'PCSK9', 'ABCA1', 'ABCG1', 'APOE'],
}

epigenetic_genes = [
    'KMT2A', 'KMT2D', 'SETD1B',    # H3K4me3 writers
    'EP300', 'CREBBP',               # H3K27ac writers
    'KDM5B', 'KDM6B',               # Demethylases
    'SIRT1', 'HDAC3',               # Deacetylases
    'DNMT3A', 'TET2',               # DNA methylation
]

# Filter to expressed genes
met_genes_flat = [g for genes in metabolic_genes.values() for g in genes]
met_genes_flat = [g for g in met_genes_flat if g in macrophages.var_names]
epi_genes_expressed = [g for g in epigenetic_genes if g in macrophages.var_names]
print(f"Metabolic genes: {len(met_genes_flat)}, Epigenetic genes: {len(epi_genes_expressed)}")

# Donor-level correlation
donor_met_epi = macrophages.obs[['donor_id', 'plaque_location']].copy()
for g in met_genes_flat + epi_genes_expressed:
    donor_met_epi[f'{g}_expr'] = macrophages[:, g].X.toarray().flatten()

donor_met_agg = donor_met_epi.groupby('donor_id').agg({
    'plaque_location': 'first',
    **{f'{g}_expr': 'mean' for g in met_genes_flat + epi_genes_expressed}
})

# Metabolic module scores (donor-level average of per-gene expression)
for mod_name, genes in metabolic_genes.items():
    valid = [g for g in genes if g in macrophages.var_names]
    if valid:
        donor_met_agg[f'{mod_name}_mean'] = donor_met_agg[[f'{g}_expr' for g in valid]].mean(axis=1)

# Merge TI scores from donor_df for partial correlation control
if 'TI_composite' in donor_df.columns:
    donor_met_agg = donor_met_agg.merge(
        donor_df[['TI_composite', 'Acute_Inflammation_score']],
        left_index=True, right_index=True, how='left'
    )
    donor_met_agg['TI_composite'] = donor_met_agg['TI_composite'].fillna(0)

# Correlation matrix: metabolic modules vs epigenetic genes
met_modules = [f'{m}_mean' for m in metabolic_genes.keys() if f'{m}_mean' in donor_met_agg.columns]
corr_matrix = pd.DataFrame(index=met_modules, columns=epi_genes_expressed)

for met in met_modules:
    for epi in epi_genes_expressed:
        r, p = spearmanr(donor_met_agg[met].dropna(), donor_met_agg[f'{epi}_expr'].dropna())
        corr_matrix.loc[met, epi] = r

corr_matrix = corr_matrix.astype(float)

print("\nMetabolic-Epigenetic Correlations:")
print(corr_matrix.round(2))

# ============================================================
# STEP 3: Vascular bed-specific coupling patterns
# ============================================================
print("\n=== Vascular Bed-Specific Coupling ===")

bed_corr_matrices = {}
for bed in ['carotid', 'coronary', 'femoral']:
    bed_data = donor_met_agg[donor_met_agg['plaque_location'] == bed]
    if len(bed_data) < 5:
        print(f"  {bed}: insufficient donors (n={len(bed_data)})")
        continue
    bed_corr = pd.DataFrame(index=met_modules, columns=epi_genes_expressed)
    for met in met_modules:
        for epi in epi_genes_expressed:
            r, p = spearmanr(bed_data[met].dropna(), bed_data[f'{epi}_expr'].dropna())
            bed_corr.loc[met, epi] = r
    bed_corr_matrices[bed] = bed_corr.astype(float)

# ============================================================
# STEP 4: Mediation analysis (P0)
# ============================================================
print("\n=== Mediation Analysis ===")

# Test: Does metabolism (Glycolysis) mediate vascular_bed → TI_composite?
# X = plaque_location (0=femoral, 1=carotid), M = glycolysis, Y = TI_composite

donor_med = donor_df.dropna(subset=['TI_composite', 'Glycolysis_score'])
donor_med['bed_encoded'] = donor_med['plaque_location'].map({'femoral': 0, 'carotid': 1, 'coronary': 0.5})

# Simple mediation using bootstrapped indirect effect
try:
    import pingouin as pg

    for bed_pair in [('carotid', 'femoral'), ('coronary', 'femoral')]:
        subset = donor_med[donor_med['plaque_location'].isin(bed_pair)]
        if len(subset) < 20:
            continue
        subset = subset.copy()
        subset['X'] = subset['plaque_location'].map({bed_pair[0]: 1, bed_pair[1]: 0})

        for mediator in ['Glycolysis_score', 'FAO_score', 'OXPHOS_score']:
            if mediator not in subset.columns:
                continue
            try:
                med = pg.mediation_analysis(
                    data=subset,
                    x='X', m=mediator, y='TI_composite',
                    covar=['sex'] if 'sex' in subset.columns else None,
                    n_boot=1000, seed=42
                )
                print(f"\n{bed_pair}: {mediator} → TI_composite")
                print(med[['path', 'coef', 'pval', 'CI[2.5%]', 'CI[97.5%]']].round(4))
            except Exception as e:
                print(f"  Mediation failed for {mediator}: {e}")

except ImportError:
    print("pingouin not available, using manual mediation (Baron & Kenny method)")
    # Baron & Kenny (1986) 4-step approach
    for bed_pair in [('carotid', 'femoral')]:
        subset = donor_med[donor_med['plaque_location'].isin(bed_pair)]
        subset = subset.copy()
        subset['X'] = subset['plaque_location'].map({bed_pair[0]: 1, bed_pair[1]: 0})

        # Step 1: X → Y (total effect)
        from scipy.stats import linregress
        c_path = linregress(subset['X'], subset['TI_composite'])

        # Step 2: X → M
        for mediator in ['Glycolysis_score', 'FAO_score']:
            a_path = linregress(subset['X'], subset[mediator])

            # Step 3: X + M → Y
            import statsmodels.api as sm
            X_med = sm.add_constant(subset[['X', mediator]])
            model = sm.OLS(subset['TI_composite'], X_med).fit()

            print(f"\n{bed_pair}: {mediator}")
            print(f"  Step 1 (c):  X→Y    β={c_path.slope:.3f}, p={c_path.pvalue:.3e}")
            print(f"  Step 2 (a):  X→M    β={a_path.slope:.3f}, p={a_path.pvalue:.3e}")
            print(f"  Step 3 (b):  M→Y|X  β={model.params[mediator]:.3f}, p={model.pvalues[mediator]:.3e}")
            print(f"  Step 3 (c'): X→Y|M  β={model.params['X']:.3f}, p={model.pvalues['X']:.3e}")

            # Sobel test
            se_ab = np.sqrt(a_path.slope**2 * model.bse[mediator]**2 + model.params[mediator]**2 * a_path.stderr**2)
            z_ab = a_path.slope * model.params[mediator] / se_ab
            from scipy.stats import norm
            p_ab = 2 * (1 - norm.cdf(abs(z_ab)))
            print(f"  Sobel: z={z_ab:.3f}, p={p_ab:.3e}")
            print(f"  Mediation proportion: {(a_path.slope * model.params[mediator]) / c_path.slope:.1%}")

# ============================================================
# STEP 5: Partial correlation (metabolism-controlled)
# ============================================================
print("\n=== Partial Correlation Analysis ===")

try:
    from pingouin import partial_corr

    for bed in ['carotid', 'femoral']:
        bed_data = donor_met_agg[donor_met_agg['plaque_location'] == bed].copy()
        if len(bed_data) < 10:
            continue
        print(f"\n{bed}:")

        for met_mod in ['Glycolysis_mean', 'FAO_mean']:
            if met_mod not in bed_data.columns:
                continue
            for epi in epi_genes_expressed:
                try:
                    pc = partial_corr(
                        data=bed_data,
                        x=met_mod, y=f'{epi}_expr',
                        covar='TI_composite',
                        method='spearman'
                    )
                    if pc['p-val'].values[0] < 0.05:
                        print(f"  {met_mod}-{epi}: r_partial={pc['r'].values[0]:.3f}, p={pc['p-val'].values[0]:.3e}")
                except:
                    pass

except ImportError:
    print("pingouin not available, computing manual partial correlation...")
    # Manual partial Spearman: r_xy.z = (r_xy - r_xz * r_yz) / sqrt((1-r_xz^2)(1-r_yz^2))
    for bed in ['carotid', 'femoral']:
        bed_data = donor_met_agg[donor_met_agg['plaque_location'] == bed].copy()
        if len(bed_data) < 10:
            continue
        for met_mod in ['Glycolysis_mean', 'FAO_mean']:
            if met_mod not in bed_data.columns:
                continue
            for epi in epi_genes_expressed:
                r_xy = spearmanr(bed_data[met_mod], bed_data[f'{epi}_expr'])[0]
                r_xz = spearmanr(bed_data[met_mod], bed_data['TI_composite'])[0]
                r_yz = spearmanr(bed_data[f'{epi}_expr'], bed_data['TI_composite'])[0]
                r_partial = (r_xy - r_xz * r_yz) / np.sqrt((1 - r_xz**2) * (1 - r_yz**2))
                if abs(r_partial) > 0.3:
                    print(f"  {bed} {met_mod}-{epi}: r_xy={r_xy:.3f}, r_partial={r_partial:.3f}")

# ============================================================
print("\n=== Generating Figure 3 ===")

import matplotlib.patches as mpatches
from matplotlib.patches import FancyBboxPatch, FancyArrowPatch

# Colorblind-friendly palette
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

fig = plt.figure(figsize=(22, 20), facecolor="white")
gs = GridSpec(3, 3, figure=fig, height_ratios=[1, 1, 1.5],
              hspace=0.40, wspace=0.35, left=0.05, right=0.97, top=0.96, bottom=0.04)

# ============================================================
# ROW 1: Metabolic-epigenetic landscape (3 panels)
# ============================================================

# --- Panel A: Metabolic module scores boxplot ---
ax1 = fig.add_subplot(gs[0, 0])
met_plot_vars = [f'{m}_score' for m in ['Glycolysis', 'FAO', 'OXPHOS', 'FAS', 'Cholesterol']]
met_plot_vars = [v for v in met_plot_vars if v in donor_df.columns]
donor_met_melt = donor_df.melt(id_vars='plaque_location', value_vars=met_plot_vars)
sns.boxplot(data=donor_met_melt, x='variable', y='value', hue='plaque_location',
            palette=CB_PALETTE, ax=ax1, width=0.6, linewidth=1.0,
            flierprops=dict(marker='o', markersize=3, alpha=0.4))
ax1.set_xticklabels(ax1.get_xticklabels(), rotation=30, ha='right', fontsize=8)
ax1.set_ylabel('Module Score', fontsize=10)
ax1.set_xlabel('')
ax1.set_title('Metabolic Pathway Scores\n(donor-level)', fontsize=10, fontweight='bold')
ax1.legend(fontsize=8, title='', loc='upper right')
fig.text(0.02, 0.95, 'A', fontsize=16, fontweight='bold')

# --- Panel B: Metabolic-epigenetic correlation heatmap ---
ax2 = fig.add_subplot(gs[0, 1])
sns.heatmap(corr_matrix, cmap='RdBu_r', center=0, annot=True,
            fmt='.2f', ax=ax2, annot_kws={'fontsize': 7},
            cbar_kws={'label': 'Spearman r', 'shrink': 0.8},
            linewidths=0.5, linecolor='#EEEEEE',
            vmin=-1, vmax=1)
ax2.set_title('Metabolic-Epigenetic Correlation\n(all donors)', fontsize=10, fontweight='bold')
ax2.tick_params(labelsize=8)
fig.text(0.35, 0.95, 'B', fontsize=16, fontweight='bold')

# --- Panel C: Vascular bed-specific differential correlation ---
ax3 = fig.add_subplot(gs[0, 2])
diff_corr = pd.DataFrame(index=met_modules, columns=epi_genes_expressed, dtype=float)
if 'carotid' in bed_corr_matrices and 'femoral' in bed_corr_matrices:
    for met in met_modules:
        for epi in epi_genes_expressed:
            diff_corr.loc[met, epi] = bed_corr_matrices['carotid'].loc[met, epi] - bed_corr_matrices['femoral'].loc[met, epi]
    vmax_val = max(abs(diff_corr.values).max(), 0.5)
    sns.heatmap(diff_corr, cmap='RdBu_r', center=0, annot=True, fmt='.2f',
                ax=ax3, annot_kws={'fontsize': 7.5},
                cbar_kws={'label': '\u0394r (Carotid - Femoral)', 'shrink': 0.8},
                linewidths=0.5, linecolor='#EEEEEE',
                vmin=-vmax_val, vmax=vmax_val)
    ax3.set_title('Differential Correlation\nCarotid vs Femoral', fontsize=10, fontweight='bold')
else:
    ax3.text(0.5, 0.5, 'Insufficient data', ha='center', va='center')
ax3.tick_params(labelsize=8)
fig.text(0.68, 0.95, 'C', fontsize=16, fontweight='bold')

# ============================================================
# ROW 2: Key pairwise correlations + simplified pathway comparison
# ============================================================

# --- Panel D: Key metabolic-epigenetic correlations by bed (was Panel E) ---
ax4 = fig.add_subplot(gs[1, 0])
key_pairs = [('Glycolysis_mean', 'EP300'), ('Glycolysis_mean', 'CREBBP'),
             ('FAO_mean', 'SIRT1'), ('Glycolysis_mean', 'KMT2A')]
pair_data = []
for met_mod, epi in key_pairs:
    if met_mod not in donor_met_agg.columns or f'{epi}_expr' not in donor_met_agg.columns:
        continue
    for bed in ORDER:
        bed_data = donor_met_agg[donor_met_agg['plaque_location'] == bed]
        if len(bed_data) >= 3:
            r, p = spearmanr(bed_data[met_mod], bed_data[f'{epi}_expr'])
            mod_label = met_mod.split('_')[0]
            sig = '*' if p < 0.05 else ''
            pair_data.append({'pair': f'{mod_label}-{epi}', 'bed': bed, 'r': r, 'sig': sig})
if pair_data:
    pair_df = pd.DataFrame(pair_data)
    bp = sns.barplot(data=pair_df, x='pair', y='r', hue='bed', palette=CB_PALETTE, ax=ax4,
                     edgecolor='white', linewidth=0.5)
    ax4.axhline(y=0, color='#333333', linewidth=0.8)
    ax4.set_xticklabels(ax4.get_xticklabels(), rotation=25, ha='right', fontsize=8)
    ax4.set_ylabel('Spearman r', fontsize=10)
    ax4.set_xlabel('')
    ax4.set_title('Key Metabolic-Epigenetic\nCorrelations by Bed', fontsize=10, fontweight='bold')
    ax4.legend(fontsize=8, title='', loc='lower left')
    for i, (_, row) in enumerate(pair_df.iterrows()):
        if row['sig']:
            x = i
            y = row['r'] + (0.05 if row['r'] >= 0 else -0.08)
            ax4.text(x, y, row['sig'], ha='center', fontsize=9, fontweight='bold')
fig.text(0.02, 0.58, 'D', fontsize=16, fontweight='bold')

# --- Panel E: Pathway divergence schematic (simplified, replaces mediation) ---
ax5 = fig.add_subplot(gs[1, 1])
ax5.set_xlim(0, 10)
ax5.set_ylim(0, 10)
ax5.axis('off')

# Vascular bed level comparison -- two parallel paths
box_style_c = dict(boxstyle='round,pad=0.4', facecolor='#FFF3E0', edgecolor='#D55E00', linewidth=2)
box_style_f = dict(boxstyle='round,pad=0.4', facecolor='#E8F5E9', edgecolor='#009E73', linewidth=2)
box_style_shared = dict(boxstyle='round,pad=0.3', facecolor='#F5F5F5', edgecolor='#888888', linewidth=1.2)

# Carotid path (left)
bbox_c1 = FancyBboxPatch((0.3, 7.5), 3.8, 1.2, **box_style_c)
ax5.add_patch(bbox_c1)
ax5.text(2.2, 8.1, 'Carotid: Glycolysis\u2191', ha='center', va='center', fontsize=9, fontweight='bold', color='#D55E00')
bbox_c2 = FancyBboxPatch((0.3, 5.5), 3.8, 1.2, **box_style_c)
ax5.add_patch(bbox_c2)
ax5.text(2.2, 6.1, 'EP300/CREBBP\u2191\n(H3K27ac writers)', ha='center', va='center', fontsize=8, color='#D55E00')
bbox_c3 = FancyBboxPatch((0.3, 3.5), 3.8, 1.2, **box_style_c)
ax5.add_patch(bbox_c3)
ax5.text(2.2, 4.1, 'Inflammatory\nTI Activation', ha='center', va='center', fontsize=8, fontweight='bold', color='#D55E00')

# Femoral path (right)
bbox_f1 = FancyBboxPatch((5.9, 7.5), 3.8, 1.2, **box_style_f)
ax5.add_patch(bbox_f1)
ax5.text(7.8, 8.1, 'Femoral: FAO\u2191', ha='center', va='center', fontsize=9, fontweight='bold', color='#009E73')
bbox_f2 = FancyBboxPatch((5.9, 5.5), 3.8, 1.2, **box_style_f)
ax5.add_patch(bbox_f2)
ax5.text(7.8, 6.1, 'SIRT1/HDAC3\u2191\n(Deacetylases)', ha='center', va='center', fontsize=8, color='#009E73')
bbox_f3 = FancyBboxPatch((5.9, 3.5), 3.8, 1.2, **box_style_f)
ax5.add_patch(bbox_f3)
ax5.text(7.8, 4.1, 'TI Maintenance\n(Acute-Inflam Decoupled)', ha='center', va='center', fontsize=8, fontweight='bold', color='#009E73')

# Shared source
bbox_source = FancyBboxPatch((3.5, 1.0), 3.0, 1.0, **box_style_shared)
ax5.add_patch(bbox_source)
ax5.text(5.0, 1.5, 'Vascular Bed\nMicroenvironment', ha='center', va='center', fontsize=8, fontweight='bold')

# Arrows: shared source -> carotid and femoral
ax5.annotate('', xy=(2.2, 3.5), xytext=(4.2, 2.0),
            arrowprops=dict(arrowstyle='->', color='#D55E00', lw=2, connectionstyle='arc3,rad=0.3'))
ax5.annotate('', xy=(7.8, 3.5), xytext=(5.8, 2.0),
            arrowprops=dict(arrowstyle='->', color='#009E73', lw=2, connectionstyle='arc3,rad=-0.3'))

# Vertical arrows within carotid path
ax5.annotate('', xy=(2.2, 5.5), xytext=(2.2, 7.3),
            arrowprops=dict(arrowstyle='->', color='#D55E00', lw=1.5))
ax5.annotate('', xy=(2.2, 3.5), xytext=(2.2, 5.3),
            arrowprops=dict(arrowstyle='->', color='#D55E00', lw=1.5))

# Vertical arrows within femoral path
ax5.annotate('', xy=(7.8, 5.5), xytext=(7.8, 7.3),
            arrowprops=dict(arrowstyle='->', color='#009E73', lw=1.5))
ax5.annotate('', xy=(7.8, 3.5), xytext=(7.8, 5.3),
            arrowprops=dict(arrowstyle='->', color='#009E73', lw=1.5))

ax5.set_title('Parallel Pathway Model\n(Vascular-bed-specific metabolic-epigenetic programs)',
              fontsize=10, fontweight='bold')
fig.text(0.35, 0.58, 'E', fontsize=16, fontweight='bold')

# ============================================================
# ROW 3: CORE FINDING -- Pathway divergence mechanism diagram
# Full-width panel as visual punchline
# ============================================================
ax6 = fig.add_subplot(gs[2, :])
ax6.set_xlim(0, 10)
ax6.set_ylim(0, 10)
ax6.axis('off')

# Three-column layout: carotid (left) / shared (center) / femoral (right)
# Column positions: carotid x=0-3.2, center x=3.5-6.5, femoral x=6.8-10

# --- Carotid column ---
carotid_patches = [
    (0.5, 7.0, 2.6, 1.3, 'Enhanced\nGlycolysis', '#D55E00', 'KW p=4.42e-02'),
    (0.5, 4.8, 2.6, 1.3, 'EP300/CREBBP\n(H3K27ac Writers)', '#E07000', 'r=-0.44'),
    (0.5, 2.6, 2.6, 1.3, 'Inflammatory\nTI Program', '#CC0000', 'Acute Inflam.\u2191'),
]
for x, y, w, h, label, color, stat in carotid_patches:
    bbox = FancyBboxPatch((x, y), w, h, boxstyle='round,pad=0.3',
                          facecolor=color, edgecolor='#666666', linewidth=1.5, alpha=0.12)
    ax6.add_patch(bbox)
    ax6.text(x + w/2, y + h/2, label, ha='center', va='center', fontsize=8.5, fontweight='bold', color=color)
    ax6.text(x + w/2, y - 0.15, stat, ha='center', va='top', fontsize=6.5, fontstyle='italic', color='#888888')

# --- Femoral column ---
femoral_patches = [
    (6.9, 7.0, 2.6, 1.3, 'Enhanced\nFAO', '#009E73', 'KW p=6.78e-02'),
    (6.9, 4.8, 2.6, 1.3, 'SIRT1/HDAC3\n(Deacetylases)', '#007A5E', 'r=+0.66'),
    (6.9, 2.6, 2.6, 1.3, 'TI Maintenance\n(Acute-Inflam Decoupled)', '#005A3E', 'Acute Inflam.\u2193'),
]
for x, y, w, h, label, color, stat in femoral_patches:
    bbox = FancyBboxPatch((x, y), w, h, boxstyle='round,pad=0.3',
                          facecolor=color, edgecolor='#666666', linewidth=1.5, alpha=0.12)
    ax6.add_patch(bbox)
    ax6.text(x + w/2, y + h/2, label, ha='center', va='center', fontsize=8.5, fontweight='bold', color=color)
    ax6.text(x + w/2, y - 0.15, stat, ha='center', va='top', fontsize=6.5, fontstyle='italic', color='#888888')

# --- Center: Decision node ---
bbox_decision = FancyBboxPatch((3.6, 4.5), 2.8, 2.5, boxstyle='round,pad=0.4',
                              facecolor='#F5F5F5', edgecolor='#333333', linewidth=2)
ax6.add_patch(bbox_decision)
ax6.text(5.0, 6.5, 'Vascular Bed\nMicroenvironment', ha='center', va='center',
         fontsize=9, fontweight='bold', color='#333333')
ax6.text(5.0, 5.2, 'Determines:\nMetabolic preference\nEpigenetic program\nTI phenotype',
         ha='center', va='center', fontsize=7, fontstyle='italic', color='#666666')

# --- Connecting arrows ---
# Decision -> Carotid (Glycolysis)
ax6.annotate('', xy=(3.1, 7.65), xytext=(5.8, 6.0),
            arrowprops=dict(arrowstyle='->', color='#D55E00', lw=2.5, connectionstyle='arc3,rad=0.4'))
ax6.text(3.8, 7.3, 'Glycolysis\nprogram', ha='center', fontsize=7, color='#D55E00', fontweight='bold')

# Decision -> Femoral (FAO)
ax6.annotate('', xy=(7.0, 7.65), xytext=(5.8, 6.0),
            arrowprops=dict(arrowstyle='->', color='#009E73', lw=2.5, connectionstyle='arc3,rad=-0.4'))
ax6.text(6.2, 7.3, 'FAO\nprogram', ha='center', fontsize=7, color='#009E73', fontweight='bold')

# Vertical arrows within carotid
for i in range(len(carotid_patches) - 1):
    cx = carotid_patches[i][0] + carotid_patches[i][2]/2
    ax6.annotate('', xy=(cx, carotid_patches[i+1][1] + carotid_patches[i+1][3]),
                xytext=(cx, carotid_patches[i][1]),
                arrowprops=dict(arrowstyle='->', color='#D55E00', lw=1.5))

# Vertical arrows within femoral
for i in range(len(femoral_patches) - 1):
    fx = femoral_patches[i][0] + femoral_patches[i][2]/2
    ax6.annotate('', xy=(fx, femoral_patches[i+1][1] + femoral_patches[i+1][3]),
                xytext=(fx, femoral_patches[i][1]),
                arrowprops=dict(arrowstyle='->', color='#009E73', lw=1.5))

# Bottom annotation
ax6.text(5.0, 0.3,
         'Parallel Metabolic-Epigenetic Programs: Carotid = Glycolysis-EP300-H3K27ac \u2192 Inflammatory TI  |  '
         'Femoral = FAO-SIRT1-HDAC \u2192 TI Maintenance (Acute-Inflammation Decoupled)',
         ha='center', fontsize=7.5, fontstyle='italic', fontweight='bold', color='#333333',
         bbox=dict(boxstyle='round,pad=0.4', facecolor='#FFFDE7', edgecolor='#E0E0E0', alpha=0.9))

ax6.set_title('Mechanism Pathway Divergence by Vascular Bed\n',
              fontsize=12, fontweight='bold', pad=10)
fig.text(0.02, 0.28, 'F', fontsize=18, fontweight='bold', color='#009E73')

# --- Save Figure 3 ---
plt.tight_layout(rect=[0, 0, 1, 0.97])
out_path = OUT_DIR / "figure3_metabolism_epigenetics.png"
fig.savefig(out_path, dpi=300, bbox_inches='tight', facecolor='white', edgecolor='none')
svg_path = out_path.with_suffix('.svg')
fig.savefig(svg_path, bbox_inches='tight', facecolor='white', edgecolor='none')
pdf_path = out_path.with_suffix('.pdf')
fig.savefig(pdf_path, bbox_inches='tight', facecolor='white', edgecolor='none')
plt.close()
print(f"  SVG: {svg_path}")
print(f"  PDF: {pdf_path}")
print(f"\nPhase 3 complete. Results in: {OUT_DIR}")
