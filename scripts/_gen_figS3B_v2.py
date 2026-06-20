"""Fig S3B v2: DEG + functional gene set enrichment with p-value annotations."""
import matplotlib as mpl
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import numpy as np, pandas as pd
from pathlib import Path
from matplotlib.patches import Patch
from matplotlib.lines import Line2D

mpl.rcParams.update({
    "font.family": "serif",
    "font.serif": ["Times New Roman", "Times", "DejaVu Serif", "serif"],
    "mathtext.fontset": "custom",
    "mathtext.rm": "Times New Roman", "mathtext.it": "Times New Roman:italic",
    "mathtext.bf": "Times New Roman:bold",
    "svg.fonttype": "none", "pdf.fonttype": 42, "font.size": 7,
    "axes.spines.right": True, "axes.spines.top": True,
    "axes.linewidth": 0.8, "axes.titleweight": "bold", "axes.titlesize": 8,
    "legend.frameon": False,
})

DATA_DIR = Path("D:/openclaw_workspace/atherosclerosis_atlas/results/figS3")
OUT_DIR = Path("D:/openclaw_workspace/atherosclerosis_atlas/results/figures")
OUT_DIR.mkdir(exist_ok=True, parents=True)

deg = pd.read_csv(DATA_DIR / "monocyte_vs_macrophage_deg.csv")
enrich = pd.read_csv(DATA_DIR / "go_enrichment_monocyte_vs_macrophage.csv")

# ── Key genes for dotplot ──
KEY_GENES = ['TREM1', 'PLIN2', 'SPP1', 'APOE', 'TREM2', 'CD68', 'CD14',
             'HLA-DRA', 'HLA-DRB1', 'TNF', 'IL1B', 'CCL2',
             'FN1', 'TYROBP', 'CPT1A', 'SIRT1', 'HMGB1', 'HK2']
CT_SHORT = ['Hub', 'Foamy', 'Inflam', 'HMOX1+', 'Other']
CT_FULL = ['PLIN2+/TREM1+ Macrophage', 'TREM2+/Foamy Macrophage',
           'Inflammatory Macrophage', 'HMOX1+ Macrophage', 'Other Macrophage']
CT_COLORS = {
    'PLIN2+/TREM1+ Macrophage': '#D55E00', 'TREM2+/Foamy Macrophage': '#0072B2',
    'Inflammatory Macrophage': '#E69F00', 'HMOX1+ Macrophage': '#009E73',
    'Other Macrophage': '#999999',
}

# ── Build dotplot matrix ──
dot_rows = []
for ct in CT_FULL:
    comp = f'{ct}_vs_Monocyte'
    sub = deg[deg['comparison'] == comp]
    for gene in KEY_GENES:
        row = sub[sub['symbol'].str.upper() == gene.upper()]
        if len(row) > 0:
            r = row.iloc[0]
            dot_rows.append({
                'gene': gene, 'subtype': ct,
                'log2FC': r['log2FC'],
                'neg_log10_padj': -np.log10(max(r['pvals_adj'], 1e-300)),
                'padj': r['pvals_adj']
            })
dot_df = pd.DataFrame(dot_rows)

# ── Figure: 2 panels ──
fig = plt.figure(figsize=(18, 9))
gs = gridspec.GridSpec(1, 2, figure=fig, wspace=0.32, left=0.07, right=0.97, top=0.92, bottom=0.10)

# ═════ Panel A: DEG dotplot with p-value annotation ═════
ax_a = fig.add_subplot(gs[0, 0])
fc_mat = pd.DataFrame(0.0, index=KEY_GENES, columns=CT_FULL)
sig_mat = pd.DataFrame(0.0, index=KEY_GENES, columns=CT_FULL)
padj_mat = pd.DataFrame(1.0, index=KEY_GENES, columns=CT_FULL)
for _, r in dot_df.iterrows():
    fc_mat.loc[r['gene'], r['subtype']] = r['log2FC']
    sig_mat.loc[r['gene'], r['subtype']] = r['neg_log10_padj']
    padj_mat.loc[r['gene'], r['subtype']] = r['padj']

y_genes = list(reversed(KEY_GENES))
x_ct = CT_FULL

for j, ct in enumerate(x_ct):
    for i, gene in enumerate(y_genes):
        fc = fc_mat.loc[gene, ct]
        sig = sig_mat.loc[gene, ct]
        padj = padj_mat.loc[gene, ct]
        if abs(fc) > 0.1:
            color = '#CC0000' if fc > 0 else '#2166AC'
            size = min(sig * 8, 120)
            alpha = min(0.3 + sig / 50, 1.0)
            # Border thickness for significance
            lw = 0.8 if padj < 0.01 else (0.4 if padj < 0.05 else 0.15)
            ax_a.scatter(j, i, s=size, c=color, edgecolors='black', linewidth=lw, alpha=alpha, zorder=3)
        else:
            ax_a.scatter(j, i, s=4, c='#DDDDDD', zorder=2)

ax_a.set_xticks(range(len(x_ct)))
ax_a.set_xticklabels(CT_SHORT, rotation=30, ha='right', fontsize=7.5)
ax_a.set_yticks(range(len(y_genes)))
ax_a.set_yticklabels(y_genes, fontsize=7.5, fontstyle='italic')
ax_a.set_title('(a) DEG: Macrophage subtypes vs Monocyte', fontsize=9, fontweight='bold')

# Legend
legend_elements = [
    Patch(color='#CC0000', label='Up in Mac vs Monocyte'),
    Patch(color='#2166AC', label='Down in Mac vs Monocyte'),
    Line2D([0], [0], marker='o', color='w', markerfacecolor='gray', markersize=6, label='padj < 0.05', markeredgecolor='black', markeredgewidth=0.4),
    Line2D([0], [0], marker='o', color='w', markerfacecolor='gray', markersize=10, label='padj < 0.01', markeredgecolor='black', markeredgewidth=0.8),
]
ax_a.legend(handles=legend_elements, fontsize=6, loc='lower left', ncol=2)
for spine in ax_a.spines.values():
    spine.set_visible(True); spine.set_linewidth(0.8)
ax_a.tick_params(direction='out', top=False, right=False, length=3, width=0.7)

# ═════ Panel B: Functional gene set enrichment ═════
ax_b = fig.add_subplot(gs[0, 1])
up_enrich = enrich[enrich['direction'] == 'up'].copy()
up_enrich['neg_log10_p'] = -np.log10(up_enrich['pval'].clip(lower=1e-50))

# Rename "gene_set" in display
up_enrich['display_set'] = up_enrich['gene_set']

gsets = sorted(up_enrich['display_set'].unique())
ct_order = CT_FULL
gs_y_pos = {gs: i for i, gs in enumerate(gsets)}
ct_x_pos = {ct: i for i, ct in enumerate(ct_order)}

for _, row in up_enrich.iterrows():
    x = ct_x_pos[row['subtype']]
    y = gs_y_pos[row['display_set']]
    # Size = overlap ratio
    sz = max(20, row['overlap'] / row['gs_size'] * 180)
    # Color = significance
    neg_log_p = row['neg_log10_p']
    if neg_log_p > 5:
        clr = '#CC0000'
    elif neg_log_p > 2:
        clr = '#E69F00'
    else:
        clr = '#4472C4'
    lw = 0.8 if neg_log_p > 5 else (0.4 if neg_log_p > 2 else 0.15)
    ax_b.scatter(x, y, s=sz, c=clr, edgecolors='black', linewidth=lw, alpha=0.85, zorder=3)
    # Annotate p-value
    p_val = row["pval"]; p_str = f"{p_val:.1e}" if p_val < 0.01 else f"{p_val:.3f}"
    ax_b.annotate(p_str, xy=(x+0.25, y), fontsize=4.5, va='center', color='#333333')

ax_b.set_xticks(range(len(ct_order)))
ax_b.set_xticklabels(CT_SHORT, rotation=30, ha='right', fontsize=7.5)
ax_b.set_yticks(range(len(gsets)))
ax_b.set_yticklabels(gsets, fontsize=7)
ax_b.set_title('(b) Functional gene set enrichment (up in Mac vs Monocyte)', fontsize=9, fontweight='bold')
ax_b.legend(handles=[
    Patch(color='#CC0000', label='p < 1e-5'),
    Patch(color='#E69F00', label='p < 0.01'),
    Patch(color='#4472C4', label='p < 0.05'),
], fontsize=6, loc='lower right')
for spine in ax_b.spines.values():
    spine.set_visible(True); spine.set_linewidth(0.8)
ax_b.tick_params(direction='out', top=False, right=False, length=3, width=0.7)

# Global title
fig.suptitle('Supplementary Fig. S3B  Monocyte-to-Macrophage Conversion: Functional Evidence',
             fontsize=10, fontweight='bold', y=0.98)

# Export
for fmt in ['png', 'svg', 'pdf']:
    path = OUT_DIR / f"FigS3B_degenrich_mono_vs_mac.{fmt}"
    fig.savefig(path, dpi=300 if fmt == 'png' else None, bbox_inches='tight', facecolor='white')
    print(f"Saved: {path}")
plt.close(fig)
print("Done.")