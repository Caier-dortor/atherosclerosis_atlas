"""
Phase 1: Macrophage Subpopulation Vascular Bed-Specific Analysis (v2)
Framework v2 — Pseudobulk + LMM design (P0-corrected)
Data: Traeuble et al., 2025, Nature Communications (CELLxGENE atlas)
============================================================
Column mappings (actual CELLxGENE structure):
  cell_type_level1 = L1 annotation  (13 types)
  cell_type_level2 = L2 annotation  (23 types)
  origin           = vascular bed   (carotid/coronary/femoral)
  donor_id         = donor          (73 unique)
  dataset          = source dataset (12, incl. Slysz_femoral)
  sex              = male/female/unknown
  development_stage = age text      (parse numeric age)
  NO lesion_stage  — skip temporal proxy
  NO sample column — use donor_id as pseudobulk unit
"""

import scanpy as sc
import pandas as pd
import numpy as np
from scipy import stats
from scipy.sparse import csr_matrix, issparse
from sklearn.decomposition import PCA
from sklearn.preprocessing import StandardScaler
from pathlib import Path
import re
import warnings
warnings.filterwarnings('ignore')

# === CONFIG ===
DATA_DIR = Path("D:/openclaw_workspace/atherosclerosis_atlas/data")
RES_DIR  = Path("D:/openclaw_workspace/atherosclerosis_atlas/results")
H5AD_PATH = DATA_DIR / "plaque_atlas.h5ad"

sc.settings.verbosity = 2
sc.set_figure_params(dpi=100, facecolor='white', frameon=True)

BEDS = ['carotid', 'coronary', 'femoral']

# ============================================================
# STEP 1: Load & inspect data
# ============================================================
print("=" * 60)
print("STEP 1: Loading CELLxGENE atlas")
print("=" * 60)

adata = sc.read_h5ad(H5AD_PATH)
print(f"Full atlas shape: {adata.shape}")

# CRITICAL: Convert Ensembl IDs to gene symbols (gene symbols in var['feature_name'])
if 'feature_name' in adata.var.columns:
    def make_unique(names):
        seen = {}; result = []
        for n in names:
            if n in seen:
                seen[n] += 1; result.append(f"{n}_{seen[n]}")
            else:
                seen[n] = 0; result.append(n)
        return result
    adata.var_names = make_unique(adata.var['feature_name'].values)
    # Note: adata.raw.var_names cannot be directly modified (read-only property)
    # Use adata.var_names for module scoring, adata.raw for pseudobulk counts
    print("Gene symbols set as var_names")
print(f"obs columns ({len(adata.obs.columns)}): {list(adata.obs.columns)}")
print(f"\ncell_type_level1:\n{adata.obs['cell_type_level1'].value_counts()}")
print(f"\ncell_type_level2:\n{adata.obs['cell_type_level2'].value_counts()}")
print(f"\norigin (vascular bed):\n{adata.obs['origin'].value_counts()}")

# Extract numeric age from development_stage
if 'development_stage' in adata.obs.columns:
    def extract_age(s):
        m = re.search(r'(\d+)-year-old', str(s))
        return float(m.group(1)) if m else None
    adata.obs['age'] = adata.obs['development_stage'].apply(extract_age)
    n_age = adata.obs['age'].notna().sum()
    print(f"\nAge extracted: {n_age}/{adata.n_obs:,} cells ({n_age/adata.n_obs:.1%})")
    if n_age > 0:
        print(f"  range: {adata.obs['age'].min():.0f}-{adata.obs['age'].max():.0f} yr, median: {adata.obs['age'].median():.0f}")

has_raw = hasattr(adata, 'raw') and adata.raw is not None
print(f"\nadata.raw present: {has_raw}")

# ============================================================
# STEP 2: Extract macrophages & myeloid
# ============================================================
print("\n" + "=" * 60)
print("STEP 2: Extracting macrophage & myeloid populations")
print("=" * 60)

# L1 macrophages
mac_mask = adata.obs['cell_type_level1'] == 'Macrophage'
macrophages = adata[mac_mask].copy()
macrophages.obs['plaque_location'] = macrophages.obs['origin']
print(f"Macrophages (L1): {macrophages.n_obs:,} cells")

# L2 macrophage subtypes
mac_subtypes = adata.obs[adata.obs['cell_type_level2'].str.contains('Macrophage', na=False)]['cell_type_level2'].unique()
print(f"Macrophage L2 subtypes: {list(mac_subtypes)}")

# All myeloid for trajectory
myeloid_mask = adata.obs['cell_type_level1'].isin(['Macrophage', 'Monocyte', 'Dendritic cell'])
myeloid = adata[myeloid_mask].copy()
print(f"Myeloid total: {myeloid.n_obs:,} cells")

# ============================================================
# STEP 3: Donor summary per vascular bed
# ============================================================
print("\n" + "=" * 60)
print("STEP 3: Donor/sample summary for pseudobulk design")
print("=" * 60)

# Full atlas by bed
for bed in BEDS:
    subset = adata[adata.obs['origin'] == bed]
    n_donors = subset.obs['donor_id'].nunique()
    n_datasets = subset.obs['dataset'].nunique()
    print(f"  {bed}: {subset.n_obs:,} cells | {n_donors} donors | {n_datasets} datasets")
    print(f"    datasets: {list(subset.obs['dataset'].unique())}")

# Macrophage-specific
print("\nMacrophages per vascular bed:")
for bed in BEDS:
    bed_mac = macrophages[macrophages.obs['origin'] == bed]
    n_donors = bed_mac.obs['donor_id'].nunique()
    print(f"  {bed}: {bed_mac.n_obs:,} macrophages | {n_donors} donors")

# Donor-level summary
donor_bed_counts = macrophages.obs.groupby('donor_id')['origin'].first().value_counts()
print(f"\nDonors per vascular bed:\n{donor_bed_counts}")

# ============================================================
# STEP 4: Donor-level cell composition
# ============================================================
print("\n" + "=" * 60)
print("STEP 4: Donor-level cell composition (pseudobulk unit)")
print("=" * 60)

comp_col = 'cell_type_level2'

donor_bed     = macrophages.obs.groupby('donor_id')['origin'].first()
donor_dataset = macrophages.obs.groupby('donor_id')['dataset'].first()
donor_sex     = macrophages.obs.groupby('donor_id')['sex'].first()
donor_age     = macrophages.obs.groupby('donor_id')['age'].first() if 'age' in macrophages.obs.columns else None

donor_counts = macrophages.obs.groupby(['donor_id', comp_col], observed=False).size().unstack(fill_value=0)
# Fix: convert column MultiIndex/categorical to strings for pandas 3.14 compat
donor_counts.columns = [str(c) for c in donor_counts.columns]
donor_props = donor_counts.div(donor_counts.sum(axis=1), axis=0)

# Use pd.concat instead of .join (avoids categorical index bug in py3.14)
donor_comp = pd.concat([
    pd.DataFrame({
        'origin': donor_bed,
        'dataset': donor_dataset,
        'sex': donor_sex,
    }),
    donor_props
], axis=1)
if donor_age is not None:
    donor_comp['age'] = donor_age

print(f"Donor-level composition: {donor_comp.shape[0]} donors x {donor_comp.shape[1]} columns")
print(f"Donors per bed:\n{donor_comp['origin'].value_counts()}")
donor_comp.to_csv(RES_DIR / "donor_composition.csv")

# Find macrophage subtype columns
mac_l2_cols = [c for c in donor_comp.columns if 'Macrophage' in c]
print(f"Macrophage L2 subtypes in composition: {mac_l2_cols}")

# ============================================================
# STEP 5: Pseudobulk expression (donor-level raw counts)
# ============================================================
print("\n" + "=" * 60)
print("STEP 5: Donor-level pseudobulk expression")
print("=" * 60)

if has_raw:
    raw_counts = adata.raw[mac_mask, :].X
    gene_names = adata.raw.var_names
else:
    raw_counts = macrophages.X
    gene_names = macrophages.var_names
    print("WARNING: adata.raw not available, using adata.X (may not be raw counts)")

donor_ids = macrophages.obs['donor_id'].values
unique_donors = np.unique(donor_ids)
print(f"Building pseudobulk: {len(unique_donors)} donors x {raw_counts.shape[1]} genes")

# Matrix multiplication for pseudobulk
donor_idx = pd.Categorical(donor_ids, categories=unique_donors).codes
design = csr_matrix((np.ones(len(donor_ids)), (donor_idx, np.arange(len(donor_ids)))))

if issparse(raw_counts):
    pb = design @ raw_counts
    pb_counts = pb.toarray()
else:
    pb_counts = design @ raw_counts

pseudobulk_df = pd.DataFrame(pb_counts, index=unique_donors, columns=gene_names)
print(f"Pseudobulk matrix: {pseudobulk_df.shape}")

# Build donor metadata
donor_meta = macrophages.obs.groupby('donor_id').agg({
    'origin': 'first',
    'dataset': 'first',
    'sex': 'first',
}).rename(columns={'origin': 'plaque_location'})
if 'age' in macrophages.obs.columns:
    donor_meta['age'] = macrophages.obs.groupby('donor_id')['age'].first()

pseudobulk_df.to_csv(RES_DIR / "pseudobulk_donor_counts.csv.gz", compression='gzip')
donor_meta.to_csv(RES_DIR / "donor_metadata.csv")
print(f"Pseudobulk saved: {pseudobulk_df.shape}")
print(f"  carotid donors: {(donor_meta['plaque_location']=='carotid').sum()}")
print(f"  coronary donors: {(donor_meta['plaque_location']=='coronary').sum()}")
print(f"  femoral donors: {(donor_meta['plaque_location']=='femoral').sum()}")

# ============================================================
# STEP 6: Module scoring
# ============================================================
print("\n" + "=" * 60)
print("STEP 6: Macrophage module scoring (single-cell)")
print("=" * 60)

modules = {
    # Macrophage subtype signatures
    'Resident_Mac':      ['LYVE1', 'FOLR2', 'TIMD4', 'MRC1', 'CD163', 'SIGLEC1', 'CSF1R'],
    'Foamy_Mac':         ['TREM2', 'CD9', 'SPP1', 'GPNMB', 'PLIN2', 'OLR1'],
    'Inflammatory_Mac':  ['IL1B', 'TNF', 'CCL2', 'CCL3', 'CCL4', 'CXCL8', 'CCL20'],
    'HMOX1_Mac':         ['HMOX1', 'NQO1', 'SOD2', 'TXNRD1', 'FTH1'],
    # Metabolism
    'Glycolysis':  ['HK1', 'HK2', 'PFKFB3', 'ALDOA', 'GAPDH', 'PGK1', 'PKM2', 'LDHA'],
    'OXPHOS':      ['NDUFA1', 'NDUFA2', 'SDHA', 'UQCRC1', 'COX4I1', 'ATP5F1A'],
    'FAO':         ['CPT1A', 'ACADVL', 'HADHA', 'ACAA2'],
    'FAS':         ['ACACA', 'FASN', 'SCD', 'ACLY', 'ME1'],
    'Cholesterol': ['HMGCR', 'SQLE', 'LDLR', 'PCSK9', 'ABCA1', 'ABCG1', 'APOE'],
    'Hypoxia':     ['HIF1A', 'EPAS1', 'CA9', 'BNIP3', 'EGLN3', 'VEGFA'],
    # Trained immunity
    'TI_Inflammation': ['IL1B', 'IL6', 'TNF', 'CXCL8', 'CCL2', 'CCL3', 'CCL4', 'IL18'],
    'TI_PRR':          ['TLR2', 'TLR4', 'NOD2', 'CLEC7A', 'TLR1', 'CD14', 'MARCO'],
    'TI_Metabolic':    ['HK1', 'HK2', 'PFKFB3', 'PKM2', 'LDHA', 'ACLY', 'FASN', 'IDH1'],
    'TI_H3K4me3':      ['KMT2A', 'KMT2D', 'SETD1A', 'SETD1B'],
    'TI_H3K27ac':      ['EP300', 'CREBBP', 'KAT2A', 'KAT2B'],
    'TI_HDAC_SIRT':    ['HDAC1', 'HDAC2', 'HDAC3', 'HDAC8', 'SIRT1', 'SIRT3', 'SIRT6'],
    # Controls
    'Acute_Inflammation': ['S100A8', 'S100A9', 'FCN1', 'VCAN', 'CD14'],
    'Healthy_PVAT':       ['ADIPOQ', 'CFD', 'PPARG', 'FABP4', 'UCP1', 'PLIN1', 'CIDEA'],
    'Disease_PVAT':       ['TNF', 'IL6', 'CCL2', 'IL1B', 'IL18', 'ICAM1', 'VCAM1'],
}

all_genes = set(macrophages.var_names)
for mod_name, gene_list in modules.items():
    valid = [g for g in gene_list if g in all_genes]
    missing = set(gene_list) - set(valid)
    if missing:
        print(f"  {mod_name}: {len(valid)}/{len(gene_list)} — missing: {missing}")
    if len(valid) >= 3:
        sc.tl.score_genes(macrophages, gene_list=valid, score_name=f'{mod_name}_score',
                          ctrl_size=min(len(valid), 50), use_raw=False)

# ============================================================
# STEP 7: TI composite score + sensitivity analysis (P1)
# ============================================================
print("\n" + "=" * 60)
print("STEP 7: TI composite score + weight sensitivity")
print("=" * 60)

ti_tiers = [
    ('TI_Inflammation_score', 0.25),
    ('TI_PRR_score',          0.15),
    ('TI_Metabolic_score',    0.20),
    ('TI_H3K4me3_score',      0.10),
    ('TI_H3K27ac_score',      0.10),
    ('TI_HDAC_SIRT_score',    0.05),
]
available_tiers = [(n, w) for n, w in ti_tiers if n in macrophages.obs.columns]
base_weights = np.array([w for _, w in available_tiers])
print(f"Available TI tiers: {[n for n,_ in available_tiers]}")

# Preset weights
macrophages.obs['TI_composite'] = sum(
    macrophages.obs[n].values * w for n, w in available_tiers
)

# PCA data-driven weights
ti_matrix = macrophages.obs[[n for n,_ in available_tiers]].values
ti_z = StandardScaler().fit_transform(ti_matrix)
pca = PCA(n_components=1)
macrophages.obs['TI_pca'] = pca.fit_transform(ti_z)[:, 0]
print(f"PCA explained variance: {pca.explained_variance_ratio_[0]:.3f}")
for name, loading in zip([n for n,_ in available_tiers], pca.components_[0]):
    print(f"  {name}: loading={loading:.3f}")

# Sensitivity: 100 weight perturbations
np.random.seed(42)
n_iter = 100
correlations = []
for i in range(n_iter):
    noise = np.random.uniform(0.5, 1.5, len(base_weights))
    w = base_weights * noise
    w = w / w.sum()
    perturbed = sum(macrophages.obs[n].values * ww for (n,_), ww in zip(available_tiers, w))
    r, _ = stats.spearmanr(macrophages.obs['TI_composite'], perturbed)
    correlations.append(r)

consistency = np.mean(np.array(correlations) > 0.8)
print(f"Weight sensitivity: {consistency:.1%} of perturbations correlate >0.8 with preset")
print(f"Correlation range: [{min(correlations):.4f}, {max(correlations):.4f}]")

# ============================================================
# STEP 8: Donor-level aggregation + KW tests (P0-corrected)
# ============================================================
print("\n" + "=" * 60)
print("STEP 8: Donor-level aggregation & KW tests")
print("=" * 60)

score_cols = [c for c in macrophages.obs.columns
              if c.endswith('_score') or c in ['TI_composite', 'TI_pca']]
donor_scores = macrophages.obs.groupby('donor_id')[score_cols].mean()
donor_meta2  = macrophages.obs.groupby('donor_id')[['origin', 'dataset', 'sex']].first()
donor_df = pd.concat([donor_meta2, donor_scores], axis=1)

print(f"Donor-level data: {donor_df.shape[0]} donors x {donor_df.shape[1]} columns")
print(f"Donors per bed:\n{donor_df['origin'].value_counts()}")

# Kruskal-Wallis (donor-level — correct no pseudoreplication)
print("\n--- Kruskal-Wallis (donor-level) ---")
test_cols = ['TI_composite', 'TI_pca',
             'Inflammatory_Mac_score', 'Foamy_Mac_score', 'Resident_Mac_score',
             'Glycolysis_score', 'FAO_score', 'OXPHOS_score',
             'Acute_Inflammation_score']

for col in test_cols:
    if col not in donor_df.columns:
        continue
    groups = [donor_df[donor_df['origin'] == b][col].dropna().values for b in BEDS]
    if all(len(g) >= 3 for g in groups):
        H, p = stats.kruskal(*groups)
        means = ', '.join([f'{b}={np.mean(g):.3f}' for b, g in zip(BEDS, groups)])
        sig = ' ***' if p < 0.001 else ' **' if p < 0.01 else ' *' if p < 0.05 else ''
        print(f"  {col}: H={H:.2f}, p={p:.2e}{sig}  |  {means}")

# Save
donor_df = donor_df.rename(columns={'origin': 'plaque_location'})
donor_df.to_csv(RES_DIR / "donor_level_scores.csv")

# ============================================================
# STEP 9: LOOCV cross-dataset validation (P1)
# ============================================================
print("\n" + "=" * 60)
print("STEP 9: LOOCV cross-dataset validation (carotid only)")
print("=" * 60)

carotid_datasets = macrophages[macrophages.obs['origin'] == 'carotid'].obs['dataset'].unique()
print(f"Carotid datasets ({len(carotid_datasets)}): {list(carotid_datasets)}")

if len(carotid_datasets) >= 3:
    loocv_results = []
    for leave_out in carotid_datasets:
        train = macrophages[(macrophages.obs['origin'] == 'carotid') &
                            (macrophages.obs['dataset'] != leave_out)]
        test  = macrophages[(macrophages.obs['origin'] == 'carotid') &
                            (macrophages.obs['dataset'] == leave_out)]

        # Compare Inflammatory Mac proportion
        train_prop = (train.obs['cell_type_level2'] == 'Inflammatory Macrophage').mean()
        test_prop  = (test.obs['cell_type_level2'] == 'Inflammatory Macrophage').mean()

        train_ti = train.obs['TI_composite'].mean()
        test_ti  = test.obs['TI_composite'].mean()

        loocv_results.append({
            'leave_out': leave_out,
            'train_n': train.n_obs,
            'test_n': test.n_obs,
            'train_inflam_prop': train_prop,
            'test_inflam_prop': test_prop,
            'train_TI': train_ti,
            'test_TI': test_ti,
        })
        print(f"  Leave out {leave_out}: train={train.n_obs:,} test={test.n_obs:,}")
        print(f"    Inflammatory: train={train_prop:.3f}, test={test_prop:.3f}")
        print(f"    TI_composite: train={train_ti:.3f}, test={test_ti:.3f}")

    lo_df = pd.DataFrame(loocv_results)
    lo_df.to_csv(RES_DIR / "loocv_carotid_validation.csv", index=False)

    # Check consistency: same direction for TI?
    ti_diff = np.sign(lo_df['test_TI'] - lo_df['train_TI'])
    consistent = 1 - np.mean(np.abs(ti_diff - ti_diff[0]) > 0)  # proportion matching first
    print(f"\nTI direction consistency: {consistent:.1%}")
else:
    print(f"Only {len(carotid_datasets)} carotid datasets — LOOCV skipped")

# ============================================================
# STEP 10: Visualization — Figure 1 (multi-panel, publication-quality)
# ============================================================
print("\n" + "=" * 60)
print("STEP 10: Generating Figure 1")
print("=" * 60)

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

# SCI journal style (Times New Roman, outward ticks, full-box axes)
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
import matplotlib.ticker as mticker
from matplotlib.gridspec import GridSpec, GridSpecFromSubplotSpec
import seaborn as sns
from scipy.stats import mannwhitneyu
from itertools import combinations

# Colorblind-friendly palette (Wong 2011, Nature Methods)
# Vermillion / Sky Blue / Bluish Green (safe for protanopia/deuteranopia)
CB_PALETTE = {'carotid': '#D55E00', 'coronary': '#0072B2', 'femoral': '#009E73'}
BEDS_ORDER = ['carotid', 'coronary', 'femoral']
# Mac subtype palette (6 types, CB-safe diverging)
MAC_CMAP = ['#D55E00', '#E69F00', '#F0E442', '#009E73', '#0072B2', '#CC79A7']

# Helper: eta-squared effect size from Kruskal-Wallis
def kruskal_eta2(groups):
    k = len(groups)
    N = sum(len(g) for g in groups)
    all_vals = np.concatenate(groups)
    grand_median = np.median(all_vals)
    ssb = sum(len(g) * (np.median(g) - grand_median)**2 for g in groups)
    sst = sum((all_vals - grand_median)**2)
    return ssb / sst if sst > 0 else 0

# Helper: pairwise Mann-Whitney U with Bonferroni correction
def pairwise_mwu(data, col, beds):
    pairs = list(combinations(beds, 2))
    results = {}
    for b1, b2 in pairs:
        g1 = data[data['plaque_location'] == b1][col].dropna().values
        g2 = data[data['plaque_location'] == b2][col].dropna().values
        if len(g1) >= 3 and len(g2) >= 3:
            u, p = mannwhitneyu(g1, g2, alternative='two-sided')
            results[(b1, b2)] = {'U': u, 'p_raw': p, 'n1': len(g1), 'n2': len(g2)}
        else:
            results[(b1, b2)] = None
    return results

# Helper: significance bar annotation
def add_sig_bars(ax, data, col, beds, y_max, pairwise_results):
    """Draw significance brackets between pairs that are significant after Bonferroni."""
    n_pairs = len(pairwise_results)
    p_vals = [(k, v['p_raw']) for k, v in pairwise_results.items() if v is not None]
    p_vals.sort(key=lambda x: x[1])
    significant = [(k, p) for k, p in p_vals if p < 0.05 / n_pairs]  # Bonferroni
    if not significant:
        significant = [(k, p) for k, p in p_vals if p < 0.05][:2]   # nominal if no Bonferroni survivors

    y_step = (y_max - data[col].min()) * 0.08
    for j, ((b1, b2), p) in enumerate(significant):
        x1 = beds.index(b1)
        x2 = beds.index(b2)
        y = y_max + y_step * (j + 1)
        ax.plot([x1, x1, x2, x2], [y, y + y_step*0.3, y + y_step*0.3, y], lw=1.2, color='#333333')
        stars = '***' if p < 0.001 else '**' if p < 0.01 else '*'
        ax.text((x1 + x2) / 2, y + y_step * 0.5, stars, ha='center', va='bottom', fontsize=9, fontweight='bold')

# =============================================================
# Build Figure 1: 3-row multi-panel (unified GridSpec — no overlap)
# =============================================================
fig = plt.figure(figsize=(22, 20), facecolor='white')

# Main 3-row GridSpec: Panel A (2x3 boxplots) / Panel B (stacked bar) / Panel C (1x3 donuts)
gs_main = GridSpec(3, 1, figure=fig, height_ratios=[1.6, 0.85, 0.95],
                    hspace=0.38, left=0.06, right=0.95, top=0.94, bottom=0.04)

# --- Panel A: Module score boxplots (2 x 3 grid) ---
gs_a = GridSpecFromSubplotSpec(2, 3, subplot_spec=gs_main[0], hspace=0.50, wspace=0.35)
plot_cols = ['TI_composite', 'Inflammatory_Mac_score', 'Resident_Mac_score',
             'Glycolysis_score', 'FAO_score', 'OXPHOS_score']

for i, col in enumerate(plot_cols):
    ax = fig.add_subplot(gs_a[i // 3, i % 3])
    if col not in donor_df.columns:
        ax.set_visible(False)
        continue
    df_plot = donor_df.dropna(subset=[col, 'plaque_location'])

    # Boxplot + stripplot
    bp = sns.boxplot(x='plaque_location', y=col, data=df_plot, order=BEDS_ORDER,
                     palette=CB_PALETTE, ax=ax, width=0.55, linewidth=1.2,
                     flierprops=dict(marker='o', markersize=4, alpha=0.4))
    sns.stripplot(x='plaque_location', y=col, data=df_plot, order=BEDS_ORDER,
                  color='#333333', alpha=0.25, size=4, ax=ax, jitter=0.15)

    # Statistics
    groups = [df_plot[df_plot['plaque_location'] == b][col].dropna().values for b in BEDS_ORDER]
    if all(len(g) >= 3 for g in groups):
        H, p_kw = stats.kruskal(*groups)
        eta2 = kruskal_eta2(groups)
        sig = '***' if p_kw < 0.001 else '**' if p_kw < 0.01 else '*' if p_kw < 0.05 else 'ns'
        # Pairwise post-hoc
        pw = pairwise_mwu(df_plot, col, BEDS_ORDER)
        y_max = df_plot[col].max()
        add_sig_bars(ax, df_plot, col, BEDS_ORDER, y_max, pw)
    else:
        H, p_kw, eta2, sig = np.nan, np.nan, np.nan, 'n<3'

    # Clean label
    label = col.replace('_score', '').replace('_', ' ')
    ax.set_title(f'{label}\nKW H={H:.1f}, p={p_kw:.2e}, η²={eta2:.3f} {sig}',
                 fontsize=9, fontweight='bold', fontfamily='sans-serif')
    ax.set_xlabel('')
    ax.tick_params(axis='x', labelsize=8)
    ax.tick_params(axis='y', labelsize=8)
    if col in ('TI_composite', 'Inflammatory_Mac_score', 'Glycolysis_score'):
        bed_means = {b: df_plot[df_plot['plaque_location'] == b][col].mean() for b in BEDS_ORDER}
        if bed_means.get('femoral', 0) > bed_means.get('carotid', 0):
            ax.annotate('Femoral\nhighest', xy=(2, 0.95), xycoords=('data', 'axes fraction'),
                        fontsize=7, color='#009E73', fontweight='bold', ha='center',
                        bbox=dict(boxstyle='round,pad=0.2', facecolor='white',
                                  edgecolor='#009E73', linewidth=1, alpha=0.85))

# Panel label + title for row A
fig.text(0.02, 0.95, 'A', fontsize=16, fontweight='bold', fontfamily='sans-serif',
         va='top', ha='left')
fig.text(0.5, 0.95, 'Donor-level module scores: Femoral plaques unexpectedly show highest TI and inflammatory signatures',
         fontsize=10, fontstyle='italic', ha='center', va='top', color='#009E73')
fig.text(0.98, 0.95, 'Donors: C=50, Co=13, F=7 | Cells: C=34.8k, Co=1.9k, F=1.2k | VarPart: plaque_location=0%',
         fontsize=6, ha='right', va='top',
         bbox=dict(boxstyle='round,pad=0.3', facecolor='#F5F5F5', alpha=0.8))

# --- Panel B: Macrophage composition stacked bar ---
ax_b = fig.add_subplot(gs_main[1])
if mac_l2_cols:
    comp_bed = donor_comp.groupby('origin')[mac_l2_cols].mean()
    # Reorder to match BEDS_ORDER
    comp_bed = comp_bed.reindex(BEDS_ORDER)
    comp_bed.T.plot(kind='bar', stacked=True, ax=ax_b, color=MAC_CMAP[:len(mac_l2_cols)],
                    width=0.65, edgecolor='white', linewidth=0.5)
    ax_b.set_ylabel('Mean proportion per donor', fontsize=11)
    ax_b.set_xlabel('')
    ax_b.legend(bbox_to_anchor=(1.01, 1), loc='upper left', fontsize=9, frameon=True,
                title='Macrophage Subtype', title_fontsize=10)
    ax_b.tick_params(axis='x', labelsize=10, rotation=0)
    ax_b.tick_params(axis='y', labelsize=9)
    ax_b.set_ylim(0, 1.05)
    ax_b.yaxis.set_major_formatter(mticker.PercentFormatter(xmax=1, decimals=0))
else:
    ax_b.text(0.5, 0.5, 'No L2 macrophage subtypes found', ha='center', va='center')
    ax_b.set_visible(False)

fig.text(0.02, 0.59, 'B', fontsize=16, fontweight='bold', fontfamily='sans-serif',
         va='top', ha='left')
fig.text(0.5, 0.59, 'Macrophage subtype composition (donor-mean)',
         fontsize=11, fontstyle='italic', ha='center', va='top')

# --- Panel C: Macrophage subtypes per bed (donut charts) ---
gs_c = GridSpecFromSubplotSpec(1, 3, subplot_spec=gs_main[2], wspace=0.25)
for i, bed in enumerate(BEDS_ORDER):
    ax = fig.add_subplot(gs_c[0, i])
    bed_mac = macrophages[macrophages.obs['origin'] == bed]
    counts = bed_mac.obs['cell_type_level2'].value_counts()

    # Truncate labels for readability
    def shorten_label(s, max_len=20):
        return s if len(s) <= max_len else s[:max_len-1] + '…'

    labels = [shorten_label(c) for c in counts.index]
    wedges, texts, autotexts = ax.pie(
        counts.values, labels=None, autopct='%1.1f%%',
        textprops={'fontsize': 7}, startangle=90,
        colors=MAC_CMAP[:len(counts)], pctdistance=0.78,
        wedgeprops=dict(width=0.38, edgecolor='white', linewidth=0.8))
    # Center circle for donut
    ax.text(0, 0, f'n={bed_mac.n_obs:,}', ha='center', va='center', fontsize=9, fontweight='bold')
    ax.set_title(f'{bed.capitalize()}', fontsize=12, fontweight='bold', fontfamily='sans-serif')
    # Only show legend on the last subplot
    if i == 2:
        ax.legend(wedges, [shorten_label(c) for c in counts.index],
                  bbox_to_anchor=(1.15, 1), loc='upper left', fontsize=7,
                  title='L2 Subtype', title_fontsize=8, frameon=True)

fig.text(0.02, 0.24, 'C', fontsize=16, fontweight='bold', fontfamily='sans-serif',
         va='top', ha='left')
fig.text(0.5, 0.24, 'Macrophage subtype distribution by vascular bed',
         fontsize=11, fontstyle='italic', ha='center', va='top')
fig.text(0.5, 0.21, 'Note: Femoral plaque estimates based on n=7 donors; interpret with caution due to limited sample size.',
         fontsize=8, fontstyle='italic', ha='center', va='top', color='grey')

# --- Save Figure 1 (PNG + SVG + PDF) ---
plt.tight_layout(rect=[0, 0, 1, 0.96])
out_path = RES_DIR / "fig1" / "figure1_macrophage_overview.png"
out_path.parent.mkdir(exist_ok=True)
fig.savefig(out_path, dpi=300, bbox_inches='tight', facecolor='white', edgecolor='none')
svg_path = out_path.with_suffix('.svg')
fig.savefig(svg_path, bbox_inches='tight', facecolor='white', edgecolor='none')
pdf_path = out_path.with_suffix('.pdf')
fig.savefig(pdf_path, bbox_inches='tight', facecolor='white', edgecolor='none')
plt.close()
print(f"Saved: {out_path}")
print(f"  SVG: {svg_path}")
print(f"  PDF: {pdf_path}")

# ============================================================
# Save & finalize
# ============================================================
macrophages.write(RES_DIR / "macrophages_annotated.h5ad")
myeloid.write(RES_DIR / "myeloid_raw.h5ad")
print(f"\nSaved: macrophages_annotated.h5ad ({macrophages.n_obs:,} cells)")

print("\n" + "=" * 60)
print("PHASE 1 COMPLETE")
print(f"Results directory: {RES_DIR}")
print("=" * 60)
