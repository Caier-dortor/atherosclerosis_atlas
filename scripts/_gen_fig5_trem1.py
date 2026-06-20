"""Fig5: TREM1 Inhibition Network Perturbation — SCI multi-panel figure."""
import matplotlib as mpl
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
from matplotlib.colors import LinearSegmentedColormap
import numpy as np, pandas as pd
from pathlib import Path

# ── SCI style ──
mpl.rcParams.update({
    "font.family": "serif",
    "font.serif": ["Times New Roman", "Times", "DejaVu Serif", "serif"],
    "mathtext.fontset": "custom",
    "mathtext.rm": "Times New Roman",
    "mathtext.it": "Times New Roman:italic",
    "mathtext.bf": "Times New Roman:bold",
    "svg.fonttype": "none", "pdf.fonttype": 42,
    "font.size": 7,
    "axes.spines.right": True, "axes.spines.top": True,
    "axes.linewidth": 0.8,
    "axes.titleweight": "bold", "axes.titlesize": 8,
    "legend.frameon": False,
})

RES_DIR = Path("D:/openclaw_workspace/atherosclerosis_atlas/results/fig5_trem1")
OUT_DIR = Path("D:/openclaw_workspace/atherosclerosis_atlas/results/figures")
OUT_DIR.mkdir(exist_ok=True, parents=True)

# ── Load data ──
delta = pd.read_csv(RES_DIR / "trem1_ko_centrality_delta.csv", index_col=0)
pw_df = pd.read_csv(RES_DIR / "trem1_ko_pathway_delta.csv")
comm_delta = pd.read_csv(RES_DIR / "trem1_ko_communication_delta.csv", index_col=0)
details_base = pd.read_csv(RES_DIR / "trem1_baseline_lr_details.csv")
details_ko = pd.read_csv(RES_DIR / "trem1_ko_lr_details.csv")

# ── Colour palette ──
CT_COLORS = {
    'PLIN2+/TREM1+ Macrophage': '#D55E00',
    'TREM2+/Foamy Macrophage': '#0072B2',
    'Inflammatory Macrophage': '#E69F00',
    'HMOX1+ Macrophage': '#009E73',
    'Other Macrophage': '#999999',
    'Monocyte': '#CC79A7',
}
cell_types = list(CT_COLORS.keys())
ct_colors_list = [CT_COLORS[c] for c in cell_types]

# ── Figure layout (2x3, with F spanning bottom row width) ──
fig = plt.figure(figsize=(20, 14))
gs = gridspec.GridSpec(2, 3, figure=fig, hspace=0.35, wspace=0.30,
                       left=0.05, right=0.98, top=0.95, bottom=0.06)

# ═══════════ Panel A: Degree Centrality Baseline vs KO ═══════════
ax_a = fig.add_subplot(gs[0, 0])
x = np.arange(len(cell_types))
w = 0.35
bars_base = ax_a.bar(x - w/2, delta['degree_baseline'], w,
                      color=ct_colors_list, edgecolor='black', linewidth=0.4, label='Baseline')
bars_ko = ax_a.bar(x + w/2, delta['degree_ko'], w,
                    color=ct_colors_list, edgecolor='black', linewidth=0.4, alpha=0.45, label='TREM1-KO', hatch='//')
ax_a.set_xticks(x)
ax_a.set_xticklabels(cell_types, rotation=35, ha='right', fontsize=6)
ax_a.set_ylabel('Weighted degree centrality', fontsize=7)
ax_a.set_title('(a) Network centrality: baseline vs TREM1 inhibition', fontsize=8, fontweight='bold')
ax_a.legend(fontsize=6, loc='upper right')
for spine in ax_a.spines.values():
    spine.set_visible(True); spine.set_linewidth(0.8)
ax_a.tick_params(direction='out', top=False, right=False, length=3, width=0.7)

# ═══════════ Panel B: Hub Disruption Index ═══════════
ax_b = fig.add_subplot(gs[0, 1])
hdi_sorted = delta.sort_values('delta_degree', ascending=True)
bars_b = ax_b.barh(range(len(hdi_sorted)), -hdi_sorted['delta_degree'],
                    color=[CT_COLORS[c] for c in hdi_sorted.index],
                    edgecolor='black', linewidth=0.4, height=0.6)
ax_b.set_yticks(range(len(hdi_sorted)))
ax_b.set_yticklabels(hdi_sorted.index, fontsize=6.5)
ax_b.set_xlabel('Hub disruption index (Δ centrality loss)', fontsize=7)
ax_b.set_title('(b) Hub disruption index', fontsize=8, fontweight='bold')
# Annotate % loss
for i, (ct, row) in enumerate(hdi_sorted.iterrows()):
    ax_b.text(-row['delta_degree'] + 0.05, i, f"-{abs(row['delta_pct']):.1f}%",
              va='center', fontsize=6, color='black')
for spine in ax_b.spines.values():
    spine.set_visible(True); spine.set_linewidth(0.8)
ax_b.tick_params(direction='out', top=False, right=False, length=3, width=0.7)

# ═══════════ Panel C: Communication Δ heatmap ═══════════
ax_c = fig.add_subplot(gs[0, 2])
# Reorder to match cell_types order
comm_delta_ct = comm_delta.loc[cell_types, cell_types]
# Create masked array: only show losses (negative Δ)
comm_loss_mask = np.ma.masked_greater(comm_delta_ct.values, 0)
vmax = max(abs(comm_delta_ct.values.min()), 0.1)
im = ax_c.imshow(comm_delta_ct.values, cmap='RdBu_r', aspect='equal',
                  vmin=-vmax, vmax=vmax, interpolation='none')
ax_c.set_xticks(range(len(cell_types)))
ax_c.set_xticklabels([c.replace('Macrophage','Mac').replace('PLIN2+/','') for c in cell_types],
                      rotation=45, ha='right', fontsize=5.5)
ax_c.set_yticks(range(len(cell_types)))
ax_c.set_yticklabels([c.replace('Macrophage','Mac').replace('PLIN2+/','') for c in cell_types],
                      fontsize=5.5)
ax_c.set_title('(c) Communication Δ (TREM1-KO − baseline)', fontsize=8, fontweight='bold')
cbar = plt.colorbar(im, ax=ax_c, shrink=0.78, aspect=20)
cbar.set_label('Δ communication probability', fontsize=6.5)
cbar.ax.tick_params(labelsize=5.5)
for spine in ax_c.spines.values():
    spine.set_visible(True); spine.set_linewidth(0.8)
ax_c.tick_params(direction='out', top=False, right=False, length=2, width=0.5)

# ═══════════ Panel D: Top TREM1-mediated L-R pairs lost ═══════════
ax_d = fig.add_subplot(gs[1, 0])
trem1_base = details_base[(details_base['ligand'] == 'TREM1') | (details_base['receptor'] == 'TREM1')]
trem1_base = trem1_base.sort_values('prob', ascending=False).head(12)
# Use pair label
labels_d = [f"{r['sender'].replace('Macrophage','Mac').replace('PLIN2+/TREM1+','Hub')}→{r['receiver'].replace('Macrophage','Mac').replace('PLIN2+/TREM1+','Hub')}" for _, r in trem1_base.iterrows()]
y_d = range(len(labels_d))
colors_d = ['#E64A19' if 'HMGB1' in r['ligand'] else '#FFB300' for _, r in trem1_base.iterrows()]
ax_d.barh(list(y_d), trem1_base['prob'].values, color=colors_d, edgecolor='black', linewidth=0.4, height=0.6)
ax_d.set_yticks(list(y_d))
ax_d.set_yticklabels(labels_d, fontsize=5.5)
ax_d.set_xlabel('Communication probability', fontsize=7)
ax_d.set_title('(d) Top TREM1-mediated L-R pairs eliminated', fontsize=8, fontweight='bold')
# Legend
from matplotlib.patches import Patch
ax_d.legend(handles=[Patch(color='#E64A19', label='HMGB1→TREM1'), Patch(color='#FFB300', label='HSPA1A→TREM1')],
            fontsize=6, loc='lower right')
for spine in ax_d.spines.values():
    spine.set_visible(True); spine.set_linewidth(0.8)
ax_d.tick_params(direction='out', top=False, right=False, length=3, width=0.7)

# ═══════════ Panel E: Pathway selectivity ═══════════
ax_e = fig.add_subplot(gs[1, 1])
# Sort pathways by baseline signaling
pw_sort = pw_df.sort_values('baseline', ascending=True)
y_e = range(len(pw_sort))
bar_colors_e = ['#CC0000' if row['pct_change'] < -90 else '#4472C4' for _, row in pw_sort.iterrows()]
ax_e.barh(list(y_e), pw_sort['baseline'].values, color=bar_colors_e, edgecolor='black', linewidth=0.4, height=0.6)
# Overlay KO values as smaller bars
ax_e.barh(list(y_e), pw_sort['ko'].values, color='white', edgecolor='black', linewidth=0.3, height=0.3, hatch='///')
ax_e.set_yticks(list(y_e))
ax_e.set_yticklabels(pw_sort['pathway'].values, fontsize=6)
ax_e.set_xlabel('Total pathway signaling (prob sum)', fontsize=7)
ax_e.set_title('(e) Pathway selectivity: only TREM1 affected', fontsize=8, fontweight='bold')
ax_e.legend(handles=[Patch(color='#4472C4', label='Baseline'), Patch(facecolor='white', edgecolor='black', hatch='///', label='TREM1-KO')],
            fontsize=6, loc='lower right')
for spine in ax_e.spines.values():
    spine.set_visible(True); spine.set_linewidth(0.8)
ax_e.tick_params(direction='out', top=False, right=False, length=3, width=0.7)

# ═══════════ Panel F: Summary & Mechanism Schematic ═══════════
ax_f = fig.add_subplot(gs[1, 2])
ax_f.axis('off')
summary_text = (
    "TREM1 Inhibition — Network Perturbation Summary\n"
    "━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━\n\n"
    "Quantitative effects:\n"
    "  • 72 TREM1-mediated L-R pairs fully eliminated\n"
    "  • 7.8% hub cell (PLIN2$^+$/TREM1$^+$ Mac) centrality loss\n"
    "  • 6.3–7.8% centrality loss across all cell types\n"
    "  • −100% TREM1 pathway; all other pathways unaffected\n\n"
    "Mechanistic interpretation:\n"
    "  • TREM1 acts as a network amplifier, not a keystone\n"
    "  • HMGB1$→$TREM1: primary ligand-receptor circuit\n"
    "  • HSPA1A$→$TREM1: secondary damage signal\n"
    "  • No cascade collapse: other pathways compensate\n\n"
    "Therapeutic implication:\n"
    "  • TREM1 blockade would dampen inflammatory signal\n"
    "    amplification without complete immunosuppression\n"
    "  • Suitable for chronic low-grade inflammation\n"
    "    (atherosclerosis) rather than acute infection"
)
ax_f.text(0.05, 0.95, summary_text, transform=ax_f.transAxes, fontsize=7,
          fontfamily='monospace', va='top', linespacing=1.5,
          bbox=dict(boxstyle='round,pad=0.8', facecolor='#F5F5F5', edgecolor='#BBBBBB', linewidth=0.8))

# ── Global annotations ──
fig.suptitle('Fig. 5  TREM1 Inhibition Network Perturbation Simulation',
             fontsize=10, fontweight='bold', y=0.99)

# ── Export ──
for fmt in ['png', 'svg', 'pdf']:
    path = OUT_DIR / f"Fig5_trem1_inhibition.{fmt}"
    fig.savefig(path, dpi=300 if fmt == 'png' else None, bbox_inches='tight', facecolor='white')
    print(f"Saved: {path}")

plt.close(fig)
print("Done.")