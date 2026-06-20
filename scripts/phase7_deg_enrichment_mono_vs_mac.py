"""Monocyte vs Macrophage DEG + GO enrichment for Fig S3 conversion mechanism."""
import scanpy as sc, numpy as np, pandas as pd
from scipy.stats import hypergeom
from pathlib import Path
import warnings
warnings.filterwarnings('ignore')

RES_DIR = Path("D:/openclaw_workspace/atherosclerosis_atlas/data")
OUT_DIR = Path("D:/openclaw_workspace/atherosclerosis_atlas/results/figS3")
OUT_DIR.mkdir(exist_ok=True, parents=True)

# ── Load full atlas, subset to Monocyte + Macrophage ──
print("Loading plaque atlas...")
adata = sc.read_h5ad(RES_DIR / "plaque_atlas.h5ad")
print(f"  {adata.n_obs} cells, {adata.n_vars} genes")

adata = adata[adata.obs['cell_type_level1'].isin(['Monocyte', 'Macrophage'])].copy()
print(f"  Monocyte + Macrophage: {adata.n_obs} cells")

# Save myeloid subset for later use
myeloid_path = RES_DIR.parent / "results" / "myeloid_mono_mac.h5ad"
adata.write(myeloid_path)
print(f"  Saved to {myeloid_path}")

# ── DEG using scanpy rank_genes_groups ──
print("\n=== Monocyte vs Macrophage DEG (scanpy) ===")
mac_subtypes = sorted(adata.obs[adata.obs['cell_type_level1'] == 'Macrophage']['cell_type_level2'].unique())
print(f"Macrophage subtypes: {mac_subtypes}")

# Create a combined group label
adata.obs['deg_group'] = adata.obs['cell_type_level2'].astype(str)
print(f"deg_group values: {sorted(adata.obs['deg_group'].unique())}")

# Run DEG: each macrophage subtype vs Monocyte
# sc.tl.rank_genes_groups operates on the AnnData directly
print("Running rank_genes_groups...")
sc.tl.rank_genes_groups(
    adata, groupby='deg_group', groups=mac_subtypes,
    reference='Monocyte', method='wilcoxon',
    n_genes=adata.n_vars, use_raw=False,
)

# ── Gene symbol mapping (var_names are Ensembl IDs) ──
ensembl_to_symbol = dict(zip(adata.var_names, adata.var['feature_name']))
# Some feature_names are NaN or empty — filter those
ensembl_to_symbol = {k: str(v) for k, v in ensembl_to_symbol.items() if str(v) not in ('nan', '', 'None')}
# Build reverse mapping (symbol → list of Ensembl IDs, for GO set mapping)
symbol_to_ensembl = {}
for eid, sym in ensembl_to_symbol.items():
    sym_upper = sym.upper()
    if sym_upper not in symbol_to_ensembl:
        symbol_to_ensembl[sym_upper] = []
    symbol_to_ensembl[sym_upper].append(eid)

all_symbols = set(ensembl_to_symbol.values())
print(f"  Genes with symbols: {len(ensembl_to_symbol)}/{adata.n_vars}")

# Extract DEG results
deg_results = []
for subtype in mac_subtypes:
    result = sc.get.rank_genes_groups_df(adata, group=subtype, key='rank_genes_groups')
    result['comparison'] = f'{subtype}_vs_Monocyte'
    deg_results.append(result)

deg_df = pd.concat(deg_results, ignore_index=True)
deg_df = deg_df.rename(columns={'names': 'gene', 'logfoldchanges': 'log2FC'})
deg_df['abs_log2FC'] = deg_df['log2FC'].abs()
# Map Ensembl to symbol
deg_df['symbol'] = deg_df['gene'].map(ensembl_to_symbol)
deg_df.to_csv(OUT_DIR / "monocyte_vs_macrophage_deg.csv", index=False)

# Top DEGs (with gene symbols)
print("\nTop up-regulated in Macrophage subtypes vs Monocyte:")
for subtype in mac_subtypes:
    sub_df = deg_df[(deg_df['comparison'] == f'{subtype}_vs_Monocyte') & (deg_df['log2FC'] > 0.5)]
    top5 = sub_df.head(5)
    top5_names = [f"{r['symbol']}" if pd.notna(r['symbol']) else r['gene'] for _, r in top5.iterrows()]
    print(f"  {subtype}: {', '.join(top5_names)}")

# ── GO Enrichment (using gene symbols) ──
print("\n=== GO Enrichment ===")
GO_SETS = {
    'Inflammatory response': ['IL1B', 'TNF', 'IL6', 'CCL2', 'CCL3', 'CCL4', 'CCL5',
                               'CXCL8', 'CXCL2', 'CXCL10', 'CXCL9', 'CXCL11',
                               'IL18', 'IL10', 'IFNG', 'NFKB1', 'RELA', 'TLR2', 'TLR4'],
    'Phagocytosis': ['CD14', 'CD68', 'FCGR1A', 'FCGR2A', 'FCGR3A', 'MRC1', 'CD163',
                     'MARCO', 'SCARB1', 'MSR1', 'CLEC7A', 'ITGAM', 'ITGAX'],
    'Lipid metabolism': ['APOE', 'APOC1', 'APOC2', 'PLIN2', 'LPL', 'LIPA', 'ABCA1',
                         'ABCG1', 'SCARB2', 'FABP4', 'FABP5', 'ACSL1', 'CD36',
                         'OLR1', 'LDLR', 'VLDLR', 'SOAT1', 'SOAT2'],
    'MHC-II antigen presentation': ['HLA-DRA', 'HLA-DRB1', 'HLA-DRB5', 'HLA-DQA1',
                                     'HLA-DQA2', 'HLA-DQB1', 'HLA-DPA1', 'HLA-DPB1',
                                     'HLA-DMA', 'HLA-DMB', 'CD74', 'CIITA'],
    'Trained immunity / IFN': ['STAT1', 'IRF1', 'IRF7', 'IRF8', 'MX1', 'MX2',
                                'ISG15', 'IFIT1', 'IFIT2', 'IFIT3', 'OAS1', 'OAS2',
                                'GBP1', 'GBP2', 'GBP5', 'RSAD2', 'IFITM1', 'IFITM2',
                                'IFITM3', 'IFI44', 'IFI44L'],
    'Foam cell differentiation': ['CD36', 'MSR1', 'OLR1', 'ABCA1', 'ABCG1', 'LPL',
                                   'APOE', 'PLIN2', 'NR1H3', 'PPARG', 'SCARB1',
                                   'FABP4', 'ACSL1', 'SREBF1', 'LDLR'],
    'Chemotaxis': ['CCR1', 'CCR2', 'CCR5', 'CXCR4', 'CX3CR1', 'CSF1R',
                   'CCL2', 'CCL3', 'CCL4', 'CCL5', 'CXCL2', 'CXCL8',
                   'CXCL10', 'CXCL12', 'CX3CL1', 'CSF1'],
    'TREM1 signaling': ['TREM1', 'TYROBP', 'HMGB1', 'SYK', 'CARD9', 'BCL10', 'MALT1',
                        'NLRC4', 'NLRP3', 'IL1B', 'TNF', 'CCL2', 'CXCL8'],
    'ECM remodeling': ['FN1', 'COL1A1', 'COL1A2', 'COL3A1', 'COL4A1', 'COL4A2',
                       'MMP2', 'MMP9', 'MMP14', 'TIMP1', 'TIMP2', 'TIMP3',
                       'SPP1', 'TGFB1', 'CTGF'],
    'FAO / Lipid oxidation': ['CPT1A', 'CPT1B', 'CPT2', 'ACADM', 'ACADL', 'ACADS',
                               'ACAA2', 'HADHA', 'HADHB', 'EHHADH', 'ACOX1', 'ACOX2',
                               'PPARA', 'PPARD', 'PPARGC1A', 'SIRT1', 'EP300'],
    'Glycolysis': ['HK1', 'HK2', 'HK3', 'GPI', 'PFKL', 'PFKM', 'PFKP',
                    'ALDOA', 'ALDOC', 'GAPDH', 'PGK1', 'PGAM1', 'ENO1', 'ENO2',
                    'PKM', 'LDHA', 'LDHB', 'SLC2A1', 'SLC2A3', 'HIF1A'],
    'Epigenetic remodeling': ['HDAC1', 'HDAC2', 'HDAC3', 'HDAC4', 'HDAC5', 'HDAC6',
                               'HDAC7', 'HDAC8', 'HDAC9', 'HDAC10', 'HDAC11',
                               'SIRT1', 'SIRT2', 'SIRT3', 'SIRT4', 'SIRT5',
                               'EP300', 'CREBBP', 'KAT2A', 'KAT2B', 'KDM1A', 'KDM5A'],
}

all_symbols_upper = set(s.upper() for s in all_symbols)
N_total = len(all_symbols_upper)

enrich_results = []
for subtype in mac_subtypes:
    sub_df = deg_df[deg_df['comparison'] == f'{subtype}_vs_Monocyte']
    # Use gene SYMBOLS for enrichment (not Ensembl IDs)
    up_mask = (sub_df['log2FC'] > 0.3) & (sub_df['pvals_adj'] < 0.05) & sub_df['symbol'].notna()
    down_mask = (sub_df['log2FC'] < -0.3) & (sub_df['pvals_adj'] < 0.05) & sub_df['symbol'].notna()
    sub_genes_up = set(sub_df.loc[up_mask, 'symbol'].str.upper())
    sub_genes_down = set(sub_df.loc[down_mask, 'symbol'].str.upper())

    for gs_name, gs_genes in GO_SETS.items():
        gs_upper = set(g.upper() for g in gs_genes)
        gs_present = gs_upper & all_symbols_upper
        if len(gs_present) < 3:
            continue

        for direction, deg_genes in [('up', sub_genes_up), ('down', sub_genes_down)]:
            if len(deg_genes) == 0:
                continue
            overlap = deg_genes & gs_present
            k = len(overlap)
            K = len(gs_present)
            n = len(deg_genes)
            p_val = hypergeom.sf(k - 1, N_total, K, n) if k > 0 else 1.0

            enrich_results.append({
                'subtype': subtype, 'gene_set': gs_name, 'direction': direction,
                'overlap': k, 'genes': ','.join(sorted(overlap)),
                'pval': p_val, 'gs_size': K, 'n_deg': n,
            })

enrich_df = pd.DataFrame(enrich_results)
if len(enrich_df) > 0:
    enrich_df.to_csv(OUT_DIR / "go_enrichment_monocyte_vs_macrophage.csv", index=False)
    print("\nTop enriched gene sets (up in Mac vs Monocyte):")
    for subtype in mac_subtypes:
        sub_enrich = enrich_df[(enrich_df['subtype'] == subtype) & (enrich_df['direction'] == 'up')]
        sub_enrich = sub_enrich.sort_values('pval')
        top3 = sub_enrich.head(3)
        for _, row in top3.iterrows():
            print(f"  {subtype:35s} | {row['gene_set']:30s}: overlap={row['overlap']}/{row['gs_size']}, p={row['pval']:.2e}")
else:
    print("  WARNING: No enrichment results — check symbol matching")

# ── Key gene stats (match by symbol) ──
print("\n=== Key gene log2FC (Mac vs Monocyte) ===")
key_genes = ['TREM1', 'PLIN2', 'APOE', 'SPP1', 'TNF', 'IL1B', 'CCL2',
             'HLA-DRA', 'HLA-DRB1', 'CD68', 'CD14', 'FN1', 'TYROBP',
             'CPT1A', 'SIRT1', 'HMGB1', 'TREM2', 'HK2', 'EP300']
for gene in key_genes:
    for subtype in mac_subtypes:
        sub_df = deg_df[(deg_df['comparison'] == f'{subtype}_vs_Monocyte') & (deg_df['symbol'].str.upper() == gene.upper())]
        if len(sub_df) > 0:
            r = sub_df.iloc[0]
            if abs(r['log2FC']) > 0.2:
                print(f"  {gene:10s} {subtype:35s}: log2FC={r['log2FC']:+.2f}, padj={r['pvals_adj']:.1e}")

print(f"\nResults saved to {OUT_DIR}")