"""GSE155512 (carotid) CellChat L-R analysis for cross-bed validation."""
import scanpy as sc, numpy as np, pandas as pd
from pathlib import Path
import networkx as nx, warnings
warnings.filterwarnings('ignore')

RES_DIR = Path("D:/openclaw_workspace/atherosclerosis_atlas/data/validation/GSE155512")
OUT_DIR = Path("D:/openclaw_workspace/atherosclerosis_atlas/results/validation")
OUT_DIR.mkdir(exist_ok=True, parents=True)

# ── Load data ──
print("Loading GSE155512...")
adata = sc.read_h5ad(RES_DIR / "GSE155512_annotated.h5ad")
print(f"  {adata.n_obs} cells, {adata.n_vars} genes")

# ── Curated immune L-R database ──
LR_PAIRS = [
    ('HMGB1', 'TREM1'), ('HSPA1A', 'TREM1'),
    ('APOE', 'TREM2'), ('APOE', 'LRP1'),
    ('CCL2', 'CCR2'), ('CCL3', 'CCR1'), ('CCL3', 'CCR5'),
    ('CCL4', 'CCR5'), ('CCL5', 'CCR5'), ('CCL5', 'CCR1'),
    ('CXCL2', 'CXCR2'), ('CXCL8', 'CXCR1'), ('CXCL8', 'CXCR2'),
    ('CXCL12', 'CXCR4'), ('CX3CL1', 'CX3CR1'),
    ('IL1B', 'IL1R1'), ('IL1B', 'IL1R2'),
    ('TNF', 'TNFRSF1A'), ('TNF', 'TNFRSF1B'),
    ('IL6', 'IL6R'), ('IL6', 'IL6ST'),
    ('IL10', 'IL10RA'), ('IL10', 'IL10RB'),
    ('IFNG', 'IFNGR1'), ('IFNG', 'IFNGR2'),
    ('IL18', 'IL18R1'),
    ('ICAM1', 'ITGAL'), ('ICAM1', 'ITGAM'), ('ICAM1', 'ITGB2'),
    ('VCAM1', 'ITGA4'), ('VCAM1', 'ITGB1'),
    ('SELL', 'CD34'),
    ('HMGB1', 'TLR2'), ('HMGB1', 'TLR4'), ('HMGB1', 'AGER'),
    ('S100A8', 'TLR4'), ('S100A9', 'TLR4'),
    ('ANXA1', 'FPR1'), ('ANXA1', 'FPR2'),
    ('SPP1', 'ITGAV'), ('SPP1', 'ITGB3'), ('SPP1', 'CD44'),
    ('C1QA', 'C1QBP'), ('C3', 'C3AR1'),
    ('DLL1', 'NOTCH1'), ('DLL4', 'NOTCH1'), ('JAG1', 'NOTCH1'),
    ('CSF1', 'CSF1R'), ('TGFB1', 'TGFBR1'), ('TGFB1', 'TGFBR2'),
    ('VEGFA', 'FLT1'), ('VEGFA', 'KDR'),
    ('HLA-DRA', 'CD4'), ('HLA-DRB1', 'CD4'),
    ('CD80', 'CD28'), ('CD80', 'CTLA4'), ('CD86', 'CD28'),
    ('CD274', 'PDCD1'),
    ('CD40', 'CD40LG'), ('CD70', 'CD27'),
    ('FASLG', 'FAS'), ('TNFSF10', 'TNFRSF10A'),
    ('LIGHT', 'TNFRSF14'),
    ('COL1A1', 'ITGB1'), ('FN1', 'ITGB1'), ('FN1', 'ITGA5'),
    ('LAMA4', 'ITGA6'), ('LAMA4', 'ITGB1'),
]

lr_valid = [(l, r) for l, r in LR_PAIRS if l in adata.var_names and r in adata.var_names]
print(f"Valid L-R pairs: {len(lr_valid)}/{len(LR_PAIRS)}")

all_ligands = sorted(set(l for l, r in lr_valid))
all_receptors = sorted(set(r for l, r in lr_valid))
all_lr_genes = sorted(set(all_ligands + all_receptors))

# ── Cell type groups ──
ct_col = 'cell_type'
cell_types = sorted(adata.obs[ct_col].unique())
print(f"Cell types: {cell_types}")

# ── Fraction expressing ──
print("\nComputing fraction expressing...")
frac_df = pd.DataFrame(index=cell_types, columns=all_lr_genes)
for ct in cell_types:
    mask = adata.obs[ct_col] == ct
    ct_data = adata[mask]
    for gene in all_lr_genes:
        frac_df.loc[ct, gene] = (ct_data[:, gene].X.toarray() > 0).mean()

# ── Communication probability matrix ──
print("Computing communication probability...")
comm_mat = pd.DataFrame(0.0, index=cell_types, columns=cell_types)
comm_details = []

for sender in cell_types:
    for receiver in cell_types:
        probs = []
        for ligand, receptor in lr_valid:
            l_frac_s = frac_df.loc[sender, ligand]
            r_frac_r = frac_df.loc[receiver, receptor]
            if l_frac_s > 0.05 and r_frac_r > 0.05:
                prob = np.sqrt(l_frac_s * r_frac_r)
                probs.append(prob)
                if prob > 0.1:
                    comm_details.append({
                        'sender': sender, 'receiver': receiver,
                        'ligand': ligand, 'receptor': receptor,
                        'prob': prob, 'l_frac': l_frac_s, 'r_frac': r_frac_r
                    })
        comm_mat.loc[sender, receiver] = np.sum(probs)

print(f"\nCommunication matrix ({len(cell_types)}x{len(cell_types)}):")
print(comm_mat.round(2))

# ── Network centrality ──
print("\n=== Network Centrality ===")
G = nx.DiGraph()
for sender in cell_types:
    for receiver in cell_types:
        w = comm_mat.loc[sender, receiver]
        if w > 0.1:
            G.add_edge(sender, receiver, weight=w)

n_nodes = len(G.nodes)
deg_cent = {}
for n in G.nodes:
    in_w = sum(d['weight'] for _, _, d in G.in_edges(n, data=True))
    out_w = sum(d['weight'] for _, _, d in G.out_edges(n, data=True))
    deg_cent[n] = (in_w + out_w) / (n_nodes - 1) if n_nodes > 1 else 0

bet_cent = nx.betweenness_centrality(G, weight='weight')
out_deg = {n: sum(d['weight'] for _, _, d in G.out_edges(n, data=True)) for n in G.nodes}
in_deg = {n: sum(d['weight'] for _, _, d in G.in_edges(n, data=True)) for n in G.nodes}

cent_df = pd.DataFrame({
    'degree_centrality': deg_cent,
    'betweenness': bet_cent,
    'outgoing_strength': out_deg,
    'incoming_strength': in_deg,
}).sort_values('degree_centrality', ascending=False)

print(cent_df.round(4))

# ── Hub cell validation ──
mac_deg = cent_df.loc['Macrophage', 'degree_centrality']
mac_bet = cent_df.loc['Macrophage', 'betweenness']
deg_rank = cent_df['degree_centrality'].rank(ascending=False)
n_types = len(cell_types)

print(f"\n=== Validation: Hub Centrality (Carotid) ===")
print(f"Macrophage degree centrality: {mac_deg:.4f} (rank {int(deg_rank['Macrophage'])}/{n_types})")
print(f"Macrophage betweenness: {mac_bet:.4f}")

# ── Key L-R pair consistency ──
print(f"\n=== L-R Consistency ===")
key_pairs = [('TREM1', 'TYROBP'), ('HMGB1', 'TREM1'), ('APOE', 'TREM2'),
             ('SPP1', 'CD44'), ('ANXA1', 'FPR1'), ('CSF1', 'CSF1R'),
             ('CCL2', 'CCR2'), ('IL1B', 'IL1R1'), ('TNF', 'TNFRSF1A')]

for lig, rec in key_pairs:
    lig_in_mac = frac_df.loc['Macrophage', lig] if lig in frac_df.columns else 0
    rec_in_mac = frac_df.loc['Macrophage', rec] if rec in frac_df.columns else 0
    lig_in_others = frac_df.drop('Macrophage')[lig].mean() if lig in frac_df.columns else 0
    rec_in_others = frac_df.drop('Macrophage')[rec].mean() if rec in frac_df.columns else 0
    print(f"  {lig}-{rec}: Mac={lig_in_mac:.1%}L/{rec_in_mac:.1%}R | Others={lig_in_others:.1%}L/{rec_in_others:.1%}R")

# ── Top L-R pairs ──
print(f"\nTop L-R pairs from Macrophage (outgoing):")
mac_outgoing = [d for d in comm_details if d['sender'] == 'Macrophage']
mac_outgoing.sort(key=lambda x: x['prob'], reverse=True)
for d in mac_outgoing[:10]:
    print(f"  Mac -> {d['receiver']:15s}: {d['ligand']}-{d['receptor']:10s} prob={d['prob']:.3f}")

print(f"\nTop L-R pairs to Macrophage (incoming):")
mac_incoming = [d for d in comm_details if d['receiver'] == 'Macrophage']
mac_incoming.sort(key=lambda x: x['prob'], reverse=True)
for d in mac_incoming[:10]:
    print(f"  {d['sender']:15s} -> Mac: {d['ligand']}-{d['receptor']:10s} prob={d['prob']:.3f}")

# ── Cross-dataset comparison with GSE131778 (coronary) ──
print("\n=== Cross-Dataset Comparison (Coronary vs Carotid) ===")
gse131778_path = OUT_DIR / "gse131778_network_centrality.csv"
if gse131778_path.exists():
    cent_coronary = pd.read_csv(gse131778_path, index_col=0)
    print(f"Loaded GSE131778 coronary centrality ({len(cent_coronary)} cell types)")

    common_types = sorted(set(cent_coronary.index) & set(cent_df.index))
    print(f"Common cell types: {common_types}")

    comp_df = pd.DataFrame({
        'coronary_degree': [cent_coronary.loc[ct, 'degree_centrality'] for ct in common_types],
        'carotid_degree': [cent_df.loc[ct, 'degree_centrality'] for ct in common_types],
        'coronary_betweenness': [cent_coronary.loc[ct, 'betweenness'] for ct in common_types],
        'carotid_betweenness': [cent_df.loc[ct, 'betweenness'] for ct in common_types],
    }, index=common_types)

    # Spearman correlation of degree centrality ranks
    from scipy.stats import spearmanr
    rho, pval = spearmanr(comp_df['coronary_degree'], comp_df['carotid_degree'])
    print(f"\nCross-dataset degree centrality correlation: rho={rho:.4f}, p={pval:.4f}")
    print(comp_df.round(4))
    comp_df.to_csv(OUT_DIR / "gse155512_cross_dataset_centrality.csv")
else:
    print("GSE131778 results not found — skipping cross-dataset comparison")

# ── Save results ──
cent_df.to_csv(OUT_DIR / "gse155512_network_centrality.csv")
pd.DataFrame(comm_details).to_csv(OUT_DIR / "gse155512_lr_details.csv", index=False)
comm_mat.to_csv(OUT_DIR / "gse155512_communication_matrix.csv")
print(f"\nResults saved to {OUT_DIR}")