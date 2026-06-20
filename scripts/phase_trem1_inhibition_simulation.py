"""TREM1 virtual knockout: network perturbation simulation (Fig5)."""
import scanpy as sc, numpy as np, pandas as pd
from pathlib import Path
import networkx as nx, warnings
warnings.filterwarnings('ignore')

RES_DIR = Path("D:/openclaw_workspace/atherosclerosis_atlas/results")
FIG5_DIR = Path("D:/openclaw_workspace/atherosclerosis_atlas/results/fig5_trem1")
FIG5_DIR.mkdir(exist_ok=True, parents=True)

# ── Load myeloid subset ──
print("Loading myeloid data...")
adata = sc.read_h5ad(RES_DIR / "myeloid_mono_mac.h5ad")
print(f"  {adata.n_obs} cells, {adata.n_vars} genes")

# ── Gene symbol ↔ Ensembl ID mapping ──
symbol_to_eid = {}
eid_to_symbol = {}
for eid, sym in zip(adata.var_names, adata.var['feature_name']):
    sym_clean = str(sym).upper()
    if sym_clean not in ('NAN', '', 'NONE'):
        symbol_to_eid[sym_clean] = eid
        eid_to_symbol[eid] = sym_clean
print(f"  Genes with symbols: {len(symbol_to_eid)}")

# ── Curated immune L-R database ──
LR_PAIRS_SYMBOL = [
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
    ('VCAM1', 'ITGA4'), ('VCAM1', 'ITGB1'), ('SELL', 'CD34'),
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
    ('FASLG', 'FAS'), ('TNFSF10', 'TNFRSF10A'), ('LIGHT', 'TNFRSF14'),
    ('COL1A1', 'ITGB1'), ('FN1', 'ITGB1'), ('FN1', 'ITGA5'),
    ('LAMA4', 'ITGA6'), ('LAMA4', 'ITGB1'),
]

# Convert to Ensembl IDs
lr_valid = []
for l_sym, r_sym in LR_PAIRS_SYMBOL:
    l_eid = symbol_to_eid.get(l_sym.upper())
    r_eid = symbol_to_eid.get(r_sym.upper())
    if l_eid and r_eid:
        lr_valid.append((l_eid, r_eid, l_sym, r_sym))

print(f"Valid L-R pairs: {len(lr_valid)}/{len(LR_PAIRS_SYMBOL)}")

all_ligands_eid = sorted(set(l for l, r, ls, rs in lr_valid))
all_receptors_eid = sorted(set(r for l, r, ls, rs in lr_valid))
all_lr_genes_eid = sorted(set(all_ligands_eid + all_receptors_eid))

# EID → symbol lookup
eid_to_lr_symbol = {}
for l_eid, r_eid, l_sym, r_sym in lr_valid:
    eid_to_lr_symbol[l_eid] = l_sym
    eid_to_lr_symbol[r_eid] = r_sym

# ── Cell types ──
cell_types = sorted(adata.obs['cell_type_level2'].unique())
print(f"Cell types: {cell_types}")

# ── Fraction expressing ──
def compute_frac(adata_subset):
    """Compute fraction expressing for L-R genes per cell type using sparse ops."""
    frac = pd.DataFrame(0.0, index=cell_types, columns=all_lr_genes_eid)
    for ct in cell_types:
        mask = adata_subset.obs['cell_type_level2'] == ct
        if mask.sum() == 0:
            continue
        # Sparse: compute fraction directly without densifying
        ct_X = adata_subset[mask, all_lr_genes_eid].X
        # csr matrix: (X > 0).mean(axis=0)
        frac_row = np.array((ct_X > 0).mean(axis=0)).flatten()
        frac.loc[ct] = frac_row
    return frac

# ── Communication matrix ──
def compute_communication(frac_df, min_frac=0.05):
    comm_mat = pd.DataFrame(0.0, index=cell_types, columns=cell_types)
    comm_details = []
    for sender in cell_types:
        for receiver in cell_types:
            probs = []
            for l_eid, r_eid, l_sym, r_sym in lr_valid:
                l_frac = frac_df.loc[sender, l_eid]
                r_frac = frac_df.loc[receiver, r_eid]
                if l_frac > min_frac and r_frac > min_frac:
                    prob = np.sqrt(l_frac * r_frac)
                    probs.append(prob)
                    if prob > 0.1:
                        comm_details.append({
                            'sender': sender, 'receiver': receiver,
                            'ligand_eid': l_eid, 'receptor_eid': r_eid,
                            'ligand': l_sym, 'receptor': r_sym,
                            'prob': prob, 'l_frac': l_frac, 'r_frac': r_frac
                        })
            comm_mat.loc[sender, receiver] = np.sum(probs)
    return comm_mat, comm_details

# ── Network metrics ──
def compute_network(comm_mat, min_weight=0.2):
    G = nx.DiGraph()
    for s in cell_types:
        for r in cell_types:
            w = comm_mat.loc[s, r]
            if w > min_weight:
                G.add_edge(s, r, weight=w)
    n = max(len(G.nodes), 1)
    deg_cent = {}
    for node in cell_types:
        if node in G.nodes:
            in_w = sum(d['weight'] for _, _, d in G.in_edges(node, data=True))
            out_w = sum(d['weight'] for _, _, d in G.out_edges(node, data=True))
            deg_cent[node] = (in_w + out_w) / (n - 1) if n > 1 else 0
        else:
            deg_cent[node] = 0.0
    bet_cent = nx.betweenness_centrality(G, weight='weight')
    cent_df = pd.DataFrame({
        'cell_type': cell_types,
        'degree_centrality': [deg_cent[c] for c in cell_types],
        'betweenness': [bet_cent.get(c, 0) for c in cell_types],
    }).set_index('cell_type')
    return G, cent_df

# ═══════════════════════════════════════════
# BASELINE: Full network
# ═══════════════════════════════════════════
print("\n=== BASELINE: Full Network ===")
frac_base = compute_frac(adata)
comm_base, details_base = compute_communication(frac_base)
G_base, cent_base = compute_network(comm_base)
cent_base = cent_base.sort_values('degree_centrality', ascending=False)
print(cent_base.round(4))

# Identify hub cell type
hub_ct = cent_base.index[0]
print(f"Hub cell type: {hub_ct} (degree={cent_base.iloc[0]['degree_centrality']:.4f})")

# ═══════════════════════════════════════════
# TREM1 INHIBITION: Zero TREM1 in fraction matrix (avoids copying AnnData)
# ═══════════════════════════════════════════
print("\n=== TREM1 INHIBITION ===")
trem1_eid = symbol_to_eid.get('TREM1')
print(f"  TREM1 Ensembl ID: {trem1_eid}")

# Simulate KO by zeroing TREM1 fraction in all cell types
frac_ko = frac_base.copy()
if trem1_eid and trem1_eid in frac_ko.columns:
    frac_ko[trem1_eid] = 0.0
    print(f"  TREM1 fraction zeroed in all {len(cell_types)} cell types")
else:
    print(f"  WARNING: TREM1 not in fraction matrix columns")

comm_ko, details_ko = compute_communication(frac_ko)
G_ko, cent_ko = compute_network(comm_ko)
cent_ko = cent_ko.sort_values('degree_centrality', ascending=False)
print(cent_ko.round(4))

# ═══════════════════════════════════════════
# Network Perturbation Analysis
# ═══════════════════════════════════════════
print("\n=== Network Perturbation (Δ after TREM1 inhibition) ===")
delta = pd.DataFrame({
    'degree_baseline': cent_base['degree_centrality'],
    'degree_ko': cent_ko['degree_centrality'],
    'betweenness_baseline': cent_base['betweenness'],
    'betweenness_ko': cent_ko['betweenness'],
}).fillna(0)
delta['delta_degree'] = delta['degree_ko'] - delta['degree_baseline']
delta['delta_pct'] = (delta['delta_degree'] / delta['degree_baseline'].replace(0, np.nan)) * 100
delta = delta.sort_values('delta_degree')
print(delta.round(4))

# ═══════════════════════════════════════════
# Communication Loss by cell-type pair
# ═══════════════════════════════════════════
print("\n=== Communication Δ (top losses) ===")
comm_delta = comm_ko - comm_base
losses = []
for sender in cell_types:
    for receiver in cell_types:
        d = comm_delta.loc[sender, receiver]
        if d < -0.01:
            losses.append({'sender': sender, 'receiver': receiver,
                          'delta': d, 'baseline': comm_base.loc[sender, receiver],
                          'ko': comm_ko.loc[sender, receiver]})
losses.sort(key=lambda x: x['delta'])
for lp in losses[:15]:
    print(f"  {lp['sender']:35s} -> {lp['receiver']:35s}: Δ={lp['delta']:.3f} (base={lp['baseline']:.3f} → ko={lp['ko']:.3f})")

# ═══════════════════════════════════════════
# Specific TREM1-mediated pairs lost
# ═══════════════════════════════════════════
print("\n=== TREM1-Mediated L-R Pairs (lost in KO) ===")
trem1_base = [d for d in details_base if d['ligand'] == 'TREM1' or d['receptor'] == 'TREM1']
trem1_ko = [d for d in details_ko if d['ligand'] == 'TREM1' or d['receptor'] == 'TREM1']
print(f"Baseline TREM1 pairs: {len(trem1_base)}, KO TREM1 pairs: {len(trem1_ko)}")
for d in sorted(trem1_base, key=lambda x: x['prob'], reverse=True):
    print(f"  {d['sender']:35s} -> {d['receiver']:35s}: {d['ligand']}-{d['receptor']:8s} prob={d['prob']:.3f}")

# ═══════════════════════════════════════════
# Hub Disruption Index
# ═══════════════════════════════════════════
print("\n=== Hub Disruption Index ===")
hdi = delta[['degree_baseline', 'degree_ko', 'delta_degree', 'delta_pct']].copy()
hdi['hub_disruption'] = -hdi['delta_degree']
hdi = hdi.sort_values('hub_disruption', ascending=False)
print(hdi.round(4))

# ═══════════════════════════════════════════
# Pathway-level analysis
# ═══════════════════════════════════════════
print("\n=== Pathway-Level Δ ===")
PATHWAYS = {
    'TREM1 signaling': ['HMGB1-TREM1', 'HSPA1A-TREM1'],
    'TREM2/Lipid': ['APOE-TREM2', 'APOE-LRP1'],
    'Chemokine': ['CCL2-CCR2', 'CCL3-CCR1', 'CCL3-CCR5', 'CCL4-CCR5', 'CCL5-CCR5',
                  'CCL5-CCR1', 'CXCL2-CXCR2', 'CXCL8-CXCR1', 'CXCL8-CXCR2',
                  'CXCL12-CXCR4', 'CX3CL1-CX3CR1'],
    'Inflammatory cytokine': ['IL1B-IL1R1', 'IL1B-IL1R2', 'TNF-TNFRSF1A', 'TNF-TNFRSF1B',
                              'IL6-IL6R', 'IL6-IL6ST', 'IL10-IL10RA', 'IL10-IL10RB',
                              'IFNG-IFNGR1', 'IFNG-IFNGR2', 'IL18-IL18R1'],
    'DAMP/TLR': ['HMGB1-TLR2', 'HMGB1-TLR4', 'HMGB1-AGER', 'S100A8-TLR4', 'S100A9-TLR4',
                 'ANXA1-FPR1', 'ANXA1-FPR2'],
    'Adhesion/Migration': ['ICAM1-ITGAL', 'ICAM1-ITGAM', 'ICAM1-ITGB2',
                           'VCAM1-ITGA4', 'VCAM1-ITGB1', 'SELL-CD34'],
    'Foamy (SPP1)': ['SPP1-ITGAV', 'SPP1-ITGB3', 'SPP1-CD44'],
    'Growth factor': ['CSF1-CSF1R', 'TGFB1-TGFBR1', 'TGFB1-TGFBR2', 'VEGFA-FLT1', 'VEGFA-KDR'],
    'MHC-II': ['HLA-DRA-CD4', 'HLA-DRB1-CD4'],
    'Complement': ['C1QA-C1QBP', 'C3-C3AR1'],
    'Notch': ['DLL1-NOTCH1', 'DLL4-NOTCH1', 'JAG1-NOTCH1'],
    'Checkpoint': ['CD80-CD28', 'CD80-CTLA4', 'CD86-CD28', 'CD274-PDCD1'],
    'ECM': ['COL1A1-ITGB1', 'FN1-ITGB1', 'FN1-ITGA5', 'LAMA4-ITGA6', 'LAMA4-ITGB1'],
}

pathway_delta = []
for pw_name, pw_pairs in PATHWAYS.items():
    baseline_sum = sum(d['prob'] for d in details_base
                       if f"{d['ligand']}-{d['receptor']}" in pw_pairs)
    ko_sum = sum(d['prob'] for d in details_ko
                 if f"{d['ligand']}-{d['receptor']}" in pw_pairs)
    delta_val = ko_sum - baseline_sum
    pct = (delta_val / baseline_sum * 100) if baseline_sum > 0.01 else 0
    pathway_delta.append({
        'pathway': pw_name, 'baseline': baseline_sum, 'ko': ko_sum,
        'delta': delta_val, 'pct_change': pct
    })

pw_df = pd.DataFrame(pathway_delta).sort_values('delta')
print(pw_df.round(3).to_string())

# ═══════════════════════════════════════════
# Save results
# ═══════════════════════════════════════════
delta.to_csv(FIG5_DIR / "trem1_ko_centrality_delta.csv")
pw_df.to_csv(FIG5_DIR / "trem1_ko_pathway_delta.csv", index=False)
comm_delta.to_csv(FIG5_DIR / "trem1_ko_communication_delta.csv")
pd.DataFrame(details_base).to_csv(FIG5_DIR / "trem1_baseline_lr_details.csv", index=False)
pd.DataFrame(details_ko).to_csv(FIG5_DIR / "trem1_ko_lr_details.csv", index=False)
print(f"\nResults saved to {FIG5_DIR}")
