"""Fig5 fixes: bootstrap CI for centrality loss + network topology + pathway compensation."""
import scanpy as sc, numpy as np, pandas as pd
from pathlib import Path
import networkx as nx
import warnings
warnings.filterwarnings('ignore')

RES_DIR = Path("D:/openclaw_workspace/atherosclerosis_atlas/results")
FIG5_DIR = Path("D:/openclaw_workspace/atherosclerosis_atlas/results/fig5_trem1")

# ── Load data ──
print("Loading myeloid data...")
adata = sc.read_h5ad(RES_DIR / "myeloid_mono_mac.h5ad")

# ── Gene symbol mapping ──
symbol_to_eid = {}
for eid, sym in zip(adata.var_names, adata.var['feature_name']):
    sym_clean = str(sym).upper()
    if sym_clean not in ('NAN', '', 'NONE'):
        symbol_to_eid[sym_clean] = eid

# ── L-R database ──
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

lr_valid = []
for l_sym, r_sym in LR_PAIRS_SYMBOL:
    l_eid = symbol_to_eid.get(l_sym.upper())
    r_eid = symbol_to_eid.get(r_sym.upper())
    if l_eid and r_eid:
        lr_valid.append((l_eid, r_eid, l_sym, r_sym))

all_lr_genes_eid = sorted(set([l for l, r, ls, rs in lr_valid] + [r for l, r, ls, rs in lr_valid]))
cell_types = sorted(adata.obs['cell_type_level2'].unique())
trem1_eid = symbol_to_eid.get('TREM1')

# ── Fraction + communication + network functions ──
def compute_frac(adata_subset):
    frac = pd.DataFrame(0.0, index=cell_types, columns=all_lr_genes_eid)
    for ct in cell_types:
        mask = adata_subset.obs['cell_type_level2'] == ct
        if mask.sum() == 0: continue
        ct_X = adata_subset[mask, all_lr_genes_eid].X
        frac.loc[ct] = np.array((ct_X > 0).mean(axis=0)).flatten()
    return frac

def compute_network_metrics(frac_df):
    comm_mat = pd.DataFrame(0.0, index=cell_types, columns=cell_types)
    for sender in cell_types:
        for receiver in cell_types:
            probs = []
            for l_eid, r_eid, l_sym, r_sym in lr_valid:
                l_frac = frac_df.loc[sender, l_eid]
                r_frac = frac_df.loc[receiver, r_eid]
                if l_frac > 0.05 and r_frac > 0.05:
                    probs.append(np.sqrt(l_frac * r_frac))
            comm_mat.loc[sender, receiver] = np.sum(probs)

    # Build graph
    G = nx.DiGraph()
    for s in cell_types:
        for r in cell_types:
            w = comm_mat.loc[s, r]
            if w > 0.2:
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
    return comm_mat, G, deg_cent

# ═══════════════════════════════════════════
# 1. BOOTSTRAP CI for centrality loss (fraction-level, memory-efficient)
# ═══════════════════════════════════════════
print("=== Bootstrap CI for centrality loss ===")

# Pre-compute per-cell-type cell counts and per-gene expression counts
ct_masks = {}
ct_n_cells = {}
ct_gene_counts = {}  # ct -> {gene: count of cells expressing}
for ct in cell_types:
    mask = adata.obs['cell_type_level2'] == ct
    ct_masks[ct] = mask
    n_cells = mask.sum()
    ct_n_cells[ct] = n_cells
    ct_X = adata[mask, all_lr_genes_eid].X
    gene_counts = np.array((ct_X > 0).sum(axis=0)).flatten()
    ct_gene_counts[ct] = dict(zip(all_lr_genes_eid, gene_counts))
    print(f"  {ct}: {n_cells} cells")

np.random.seed(42)
n_boot = 1000
boot_deltas = {ct: [] for ct in cell_types}

for i in range(n_boot):
    if i % 200 == 0: print(f"  Bootstrap {i}/{n_boot}...")
    # Sample fraction matrix by drawing from binomial(n_cells, p_hat)/n_cells for each gene
    frac_b = pd.DataFrame(0.0, index=cell_types, columns=all_lr_genes_eid)
    for ct in cell_types:
        n_cells = ct_n_cells[ct]
        if n_cells == 0: continue
        for eid in all_lr_genes_eid:
            k_obs = ct_gene_counts[ct].get(eid, 0)
            p_hat = k_obs / n_cells
            # Bootstrap: draw from Binomial(n_cells, p_hat)
            k_boot = np.random.binomial(n_cells, p_hat)
            frac_b.loc[ct, eid] = k_boot / n_cells

    # Baseline network
    _, _, deg_b = compute_network_metrics(frac_b)

    # TREM1 KO
    frac_k = frac_b.copy()
    if trem1_eid in frac_k.columns:
        frac_k[trem1_eid] = 0.0
    _, _, deg_k = compute_network_metrics(frac_k)

    for ct in cell_types:
        boot_deltas[ct].append(deg_k[ct] - deg_b[ct])

# Compute bootstrap CIs
boot_ci = {}
for ct in cell_types:
    deltas = np.array(boot_deltas[ct])
    ci_lo = np.percentile(deltas, 2.5)
    ci_hi = np.percentile(deltas, 97.5)
    mean_delta = np.mean(deltas)
    p_emp = (deltas >= 0).mean()
    boot_ci[ct] = {'mean_delta': mean_delta, 'ci_lo': ci_lo, 'ci_hi': ci_hi, 'p_emp': p_emp}
    sig_mark = '**' if p_emp < 0.01 else ('*' if p_emp < 0.05 else 'ns')
    print(f"  {ct:35s}: Δ={mean_delta:.3f}, 95%CI=[{ci_lo:.3f}, {ci_hi:.3f}], p_boot={p_emp:.3f} {sig_mark}")

boot_df = pd.DataFrame(boot_ci).T
boot_df.to_csv(FIG5_DIR / "trem1_ko_bootstrap_ci.csv")

# ═══════════════════════════════════════════
# 2. NETWORK TOPOLOGY (baseline vs KO)
# ═══════════════════════════════════════════
print("\n=== Network Topology ===")
frac_base = compute_frac(adata)
comm_base, G_base, deg_base = compute_network_metrics(frac_base)
frac_ko = frac_base.copy()
if trem1_eid in frac_ko.columns:
    frac_ko[trem1_eid] = 0.0
comm_ko, G_ko, deg_ko = compute_network_metrics(frac_ko)

# Save edge lists for visualization
def save_edges(G, path):
    edges = []
    for u, v, d in G.edges(data=True):
        edges.append({'source': u, 'target': v, 'weight': d['weight']})
    pd.DataFrame(edges).to_csv(path, index=False)

save_edges(G_base, FIG5_DIR / "network_edges_baseline.csv")
save_edges(G_ko, FIG5_DIR / "network_edges_ko.csv")

# Node metrics
def save_nodes(deg_cent, comm_mat, path):
    nodes = []
    for ct in cell_types:
        nodes.append({
            'cell_type': ct,
            'degree_centrality': deg_cent.get(ct, 0),
            'total_outgoing': comm_mat.loc[ct].sum(),
            'total_incoming': comm_mat[ct].sum(),
        })
    pd.DataFrame(nodes).to_csv(path, index=False)

save_nodes(deg_base, comm_base, FIG5_DIR / "network_nodes_baseline.csv")
save_nodes(deg_ko, comm_ko, FIG5_DIR / "network_nodes_ko.csv")

print(f"Baseline: {G_base.number_of_nodes()} nodes, {G_base.number_of_edges()} edges")
print(f"KO:       {G_ko.number_of_nodes()} nodes, {G_ko.number_of_edges()} edges")
print(f"Edge loss: {G_base.number_of_edges() - G_ko.number_of_edges()}")

# ═══════════════════════════════════════════
# 3. PATHWAY COMPENSATION QUANTIFICATION
# ═══════════════════════════════════════════
print("\n=== Pathway Compensation ===")
# Group L-R pairs by pathway
PATHWAYS = {
    'TREM1': ['HMGB1-TREM1', 'HSPA1A-TREM1'],
    'TREM2/Lipid': ['APOE-TREM2', 'APOE-LRP1'],
    'Chemokine': ['CCL2-CCR2', 'CCL3-CCR1', 'CCL3-CCR5', 'CCL4-CCR5', 'CCL5-CCR5',
                  'CCL5-CCR1', 'CXCL2-CXCR2', 'CXCL8-CXCR1', 'CXCL8-CXCR2',
                  'CXCL12-CXCR4', 'CX3CL1-CX3CR1'],
    'Inflammatory cytokine': ['IL1B-IL1R1', 'IL1B-IL1R2', 'TNF-TNFRSF1A', 'TNF-TNFRSF1B',
                              'IL6-IL6R', 'IL6-IL6ST', 'IL10-IL10RA', 'IL10-IL10RB',
                              'IFNG-IFNGR1', 'IFNG-IFNGR2', 'IL18-IL18R1'],
    'DAMP/TLR': ['HMGB1-TLR2', 'HMGB1-TLR4', 'HMGB1-AGER', 'S100A8-TLR4', 'S100A9-TLR4',
                 'ANXA1-FPR1', 'ANXA1-FPR2'],
    'Adhesion': ['ICAM1-ITGAL', 'ICAM1-ITGAM', 'ICAM1-ITGB2', 'VCAM1-ITGA4', 'VCAM1-ITGB1', 'SELL-CD34'],
    'Foamy (SPP1)': ['SPP1-ITGAV', 'SPP1-ITGB3', 'SPP1-CD44'],
    'Growth factor': ['CSF1-CSF1R', 'TGFB1-TGFBR1', 'TGFB1-TGFBR2', 'VEGFA-FLT1', 'VEGFA-KDR'],
    'MHC-II': ['HLA-DRA-CD4', 'HLA-DRB1-CD4'],
    'Complement': ['C1QA-C1QBP', 'C3-C3AR1'],
    'Notch': ['DLL1-NOTCH1', 'DLL4-NOTCH1', 'JAG1-NOTCH1'],
    'Checkpoint': ['CD80-CD28', 'CD80-CTLA4', 'CD86-CD28', 'CD274-PDCD1'],
    'ECM': ['COL1A1-ITGB1', 'FN1-ITGB1', 'FN1-ITGA5', 'LAMA4-ITGA6', 'LAMA4-ITGB1'],
}

# Compute per-pathway total communication (sum over all sender→receiver)
def pathway_totals(frac_df):
    pw_totals = {}
    for pw_name, pw_pairs in PATHWAYS.items():
        total = 0.0
        for sender in cell_types:
            for receiver in cell_types:
                for l_eid, r_eid, l_sym, r_sym in lr_valid:
                    pair_str = f"{l_sym}-{r_sym}"
                    if pair_str in pw_pairs:
                        l_frac = frac_df.loc[sender, l_eid]
                        r_frac = frac_df.loc[receiver, r_eid]
                        if l_frac > 0.05 and r_frac > 0.05:
                            total += np.sqrt(l_frac * r_frac)
        pw_totals[pw_name] = total
    return pw_totals

pw_base = pathway_totals(frac_base)
pw_ko = pathway_totals(frac_ko)

# Compensation analysis
print(f"{'Pathway':25s} {'Baseline':>8s} {'KO':>8s} {'Δ':>8s} {'%Δ':>8s} {'Compensation':>14s}")
comp_data = []
for pw_name in pw_base:
    b = pw_base[pw_name]
    k = pw_ko.get(pw_name, 0)
    d = k - b
    pct = (d / b * 100) if b > 0.01 else 0

    # Compensation: if other pathways increase when TREM1 is removed
    if pw_name != 'TREM1' and d > 0:
        comp_note = f"+{d:.3f} (compensatory)"
    elif pw_name == 'TREM1':
        comp_note = "TARGET ABLATED"
    elif d < 0:
        comp_note = f"{d:.3f} (minor loss)"
    else:
        comp_note = "unchanged"

    comp_data.append({'pathway': pw_name, 'baseline': b, 'ko': k,
                      'delta': d, 'pct_change': pct, 'compensation': comp_note})
    print(f"  {pw_name:25s} {b:8.3f} {k:8.3f} {d:+8.3f} {pct:+7.1f}% {comp_note:>14s}")

comp_df = pd.DataFrame(comp_data)
comp_df.to_csv(FIG5_DIR / "pathway_compensation.csv", index=False)

# Loss absorbed by other pathways?
trem1_loss = abs(pw_base['TREM1'] - pw_ko.get('TREM1', 0))
non_trem1_delta = sum(max(0, pw_ko.get(pw, 0) - pw_base[pw]) for pw in pw_base if pw != 'TREM1')
print(f"\nTREM1 pathway loss: {trem1_loss:.3f}")
print(f"Compensatory gain in other pathways: {non_trem1_delta:.3f}")
print(f"Compensation ratio: {non_trem1_delta/trem1_loss*100:.1f}% of TREM1 loss absorbed")

print(f"\nResults saved to {FIG5_DIR}")