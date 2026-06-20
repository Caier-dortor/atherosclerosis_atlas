"""Apply all supplementary figure review fixes. Run once to update all scripts."""
import re

# ============================================================
# FigS4 fixes (Issues S3, S7, S14, S15)
# ============================================================
path = r'D:\openclaw_workspace\atherosclerosis_atlas\scripts\_gen_figS4_merged.py'
with open(path, 'r', encoding='utf-8') as f:
    content = f.read()

# S3: Fix N.D. legend text
old_nd = "N.D. = Not Detected (gene absent in dataset)"
new_nd = "N.D. = Not Detected (gene absent in respective dataset; CCL2/CSF1/SPP1/HMGB1 absent in Coronary)"
content = content.replace(old_nd, new_nd)

# S7: Add "ns" to non-significant p-value in Panel E
old_e_title = "ax_e.set_title(f'(e) Cross-dataset centrality: Spearman rho = {rho:.2f} (p = {pval:.2f})',\n               fontsize=8, fontweight='bold')"
new_e_title = "ns_flag = ', ns' if pval >= 0.05 else ''\nax_e.set_title(f'(e) Cross-dataset centrality: Spearman rho = {rho:.2f} (p = {pval:.2f}{ns_flag})',\n               fontsize=8, fontweight='bold')"
content = content.replace(old_e_title, new_e_title)

# S14: Fix conclusion text — remove femoral claim from FigS4
old_conclusion = '"Conclusion: PLIN2+/TREM1+ hub cell mechanism\\n"\n    "independently validated in 2 external cohorts\\n"\n    "across 2 vascular beds (coronary + carotid).\\n"\n    "FAO-SIRT1-MHC-II axis is vascular-bed-specific\\n"\n    "(femoral), not universal."'
new_conclusion = '"Conclusion: Macrophage hub architecture and TREM1/\\n"\n    "HMGB1 signaling are conserved across coronary and\\n"\n    "carotid arteries, independently validating the\\n"\n    "primary atlas findings in 2 external cohorts."'
content = content.replace(old_conclusion, new_conclusion)

# S15: Ensure y-axis label fits
old_y_c = "ax_c.set_ylabel('Weighted degree centrality', fontsize=7)"
new_y_c = "ax_c.set_ylabel('Weighted degree\\ncentrality', fontsize=7)"
content = content.replace(old_y_c, new_y_c)

with open(path, 'w', encoding='utf-8') as f:
    f.write(content)
print("FigS4 fixes applied: S3 (N.D. legend), S7 (ns flag), S14 (conclusion), S15 (label)")

# ============================================================
# FigS3A fixes (Issues S1, S2, S12)
# ============================================================
path2 = r'D:\openclaw_workspace\atherosclerosis_atlas\scripts\phase7_pseudotime_trajectory.py'
with open(path2, 'r', encoding='utf-8') as f:
    content2 = f.read()

# S1: Add fallback genes for Panel C to ensure it's never empty
old_genes = """genes_to_plot = [g for g in ['CD14', 'FCGR3A', 'IL1B', 'TNF', 'TREM2', 'APOE', 'HMOX1', 'SPP1'] if g in gene_trends]
n_genes = len(genes_to_plot)"""

new_genes = """genes_to_plot = [g for g in ['TREM1', 'PLIN2', 'CD14', 'CD68', 'HLA-DRA', 'PPARG',
                   'IL1B', 'TNF', 'TREM2', 'APOE', 'HMOX1', 'SPP1', 'FCGR3A'] if g in gene_trends]
if len(genes_to_plot) == 0:
    print('WARNING: No requested genes in gene_trends, using all available')
    genes_to_plot = list(gene_trends.keys())[:8]
n_genes = len(genes_to_plot)
print(f'Panel C: plotting {n_genes} genes along pseudotime')"""

content2 = content2.replace(old_genes, new_genes)

# S12: Fix PAGA labels — larger, better placement
old_paga = "sc.pl.paga(mye, ax=ax2, show=False, node_size_scale=1.5, edge_width_scale=0.5,\n           title='PAGA: Myeloid Lineage Topology')\nax2.set_title('PAGA: Myeloid Lineage\\nTopology', fontsize=10, fontweight='bold')"
new_paga = "sc.pl.paga(mye, ax=ax2, show=False, node_size_scale=1.8, edge_width_scale=0.8,\n           title='PAGA: Myeloid Lineage Topology', fontsize=8,\n           node_label_size=6)\nax2.set_title('PAGA: Myeloid Lineage\\nTopology', fontsize=10, fontweight='bold')"
content2 = content2.replace(old_paga, new_paga)

# S2: Fix mathtext to avoid LaTeX garbage in PDF — use dejavusans fallback
old_mpl = """mpl.rcParams.update({
    'font.family': 'serif',
    'font.serif': ['Times New Roman', 'Times', 'DejaVu Serif', 'serif'],
    'mathtext.fontset': 'custom',
    'mathtext.rm': 'Times New Roman',
    'mathtext.it': 'Times New Roman:italic',
    'mathtext.bf': 'Times New Roman:bold',
    'font.size': 8,"""

new_mpl = """mpl.rcParams.update({
    'font.family': 'serif',
    'font.serif': ['Times New Roman', 'Times', 'DejaVu Serif', 'serif'],
    'mathtext.fontset': 'dejavuserif',
    'pdf.fonttype': 42,
    'font.size': 8,"""

# Only replace the specific instance in phase7_pseudotime_trajectory.py
if old_mpl in content2:
    content2 = content2.replace(old_mpl, new_mpl)
    print("FigS3A: Fixed mathtext settings for clean PDF output")

with open(path2, 'w', encoding='utf-8') as f:
    f.write(content2)
print("FigS3A fixes applied: S1 (fallback genes), S12 (PAGA labels), S2 (mathtext)")

print("\n=== All supplementary fixes applied ===")