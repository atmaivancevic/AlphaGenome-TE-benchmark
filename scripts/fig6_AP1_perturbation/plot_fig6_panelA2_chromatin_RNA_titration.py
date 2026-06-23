"""
Fig 6 Panel A2 — AP1 motif titration: AG-predicted downstream chromatin +
target gene (ATG12) RNA. Companion to plot_fig6_panelA1_TF_titration.py
(which handles the upstream TF binding panel).

Three stacked subplots sharing one x-axis (intact AP1 motifs 0-200), each
with its own raw y-axis (no normalisation, no twin axes):
  1. Top:    H3K27ac (green) + H3K4me1 (purple) sharing one y-axis
             (~0–3500; H3K4me1 occupies the lower portion since its max
             is ~28% of H3K27ac's)
  2. Middle: ATAC (light blue) own y-axis (~0–2.5)
  3. Bottom: ATG12 RNA-seq (dark blue) own y-axis (~0.158–0.170)

The "all chromatin marks saturate at ~50 motifs" pattern is visible
across all three subplots independently — strongest possible visual
support for the training-distribution-ceiling claim.

H3K4me1 and ATAC actually drop BELOW their WT (20 motif) baselines past
~100 motifs — AG predicts overshooting AP1 density is detrimental to those
marks, an even stronger version of the training-distribution ceiling than
H3K27ac (which plateaus high rather than dropping).

Usage:
    python scripts/fig6_AP1_perturbation/plot_fig6_panelA2_chromatin_RNA_titration.py \\
        --chromatin-csv \\
            results/AG_perturbation_LTR10_ATG12_AP1/predict_chromatin_at_element.csv \\
            results/AG_perturbation_LTR10_ATG12_AP1/predict_chromatin_at_element_INS.csv \\
        --rna-csv \\
            results/AG_perturbation_LTR10_ATG12_AP1/predict_RNA_at_ATG12_gene.csv \\
            results/AG_perturbation_LTR10_ATG12_AP1/predict_RNA_at_ATG12_gene_INS.csv \\
        --output figures/FIG6_FINAL/panelA2_chromatin_RNA_titration.pdf
"""
import argparse
from pathlib import Path
import matplotlib
matplotlib.use('Agg')
matplotlib.rcParams['font.family'] = 'DejaVu Sans'
matplotlib.rcParams['pdf.fonttype'] = 42
import matplotlib.pyplot as plt
import pandas as pd

ap = argparse.ArgumentParser()
ap.add_argument('--chromatin-csv', nargs='+', required=True)
ap.add_argument('--rna-csv', nargs='+', required=True)
ap.add_argument('--output', required=True)
ap.add_argument('--width',  type=float, default=5.5, help='Figure width inches')
ap.add_argument('--height', type=float, default=7.0, help='Figure height inches')
args = ap.parse_args()

PALETTE = {
    'H3K27ac':   '#5BA75A',  # green
    'H3K4me1':   '#7B4FA1',  # purple
    'ATAC':      '#56B4E9',  # light blue
    'ATG12 RNA': '#1F4E79',  # dark blue
}
FS_AXIS, FS_TICK, FS_LEGEND, FS_WT = 10, 9, 8, 8
LW, MS = 2.0, 6

TITRATION = {
    'WT':                    20,
    'LTR10.ATG12_scram05':   15,
    'LTR10.ATG12_scram10':   10,
    'LTR10.ATG12_scram15':    5,
    'LTR10.ATG12_scram20':    0,
    'LTR10.ATG12_add05':     25,
    'LTR10.ATG12_add10':     30,
    'LTR10.ATG12_add20':     40,
    'LTR10.ATG12_add30':     50,
    'LTR10.ATG12_add40':     60,
    'LTR10.ATG12_add50':     70,
    'LTR10.ATG12_add60':     80,
    'LTR10.ATG12_add70':     90,
    'LTR10.ATG12_add80':    100,
    'LTR10.ATG12_add90':    110,
    'LTR10.ATG12_add100':   120,
    'LTR10.ATG12_add110':   130,
    'LTR10.ATG12_add120':   140,
    'LTR10.ATG12_add130':   150,
    'LTR10.ATG12_add140':   160,
    'LTR10.ATG12_add150':   170,
    'LTR10.ATG12_add160':   180,
    'LTR10.ATG12_add170':   190,
    'LTR10.ATG12_add180':   200,
    # Note: add190/200/210/220 alleles exist in the variant tab (210-240
    # motifs) but are excluded — at the user's call we cap x-axis at 200
    # motifs since the plateau is fully established by 100.
}

def raw_curve(csv_paths, track_substring, keep_stranded_only=False):
    dfs = [pd.read_csv(p) for p in csv_paths]
    df = pd.concat(dfs, ignore_index=True)
    df = df[df['track_name'].str.contains(track_substring, case=False, na=False)]
    if keep_stranded_only:
        df = df[~df['track_name'].str.contains('polyA', case=False, na=False)]
    per_var = df.groupby('variant_id').agg(
        sig_REF=('signal_mean_REF', 'mean'),
        sig_ALT=('signal_mean_ALT', 'mean'),
    )
    sig_REF = per_var['sig_REF'].mean()
    out = {20: sig_REF}
    for allele, motifs in TITRATION.items():
        if allele == 'WT': continue
        if allele in per_var.index:
            out[motifs] = float(per_var.loc[allele, 'sig_ALT'])
    return out, sig_REF

h3k_curve,    h3k_ref    = raw_curve(args.chromatin_csv, 'H3K27ac')
h3k4me1_curve, h3k4me1_ref = raw_curve(args.chromatin_csv, 'H3K4me1')
atac_curve,   atac_ref   = raw_curve(args.chromatin_csv, 'ATAC')
rna_curve,    rna_ref    = raw_curve(args.rna_csv, 'total RNA-seq', keep_stranded_only=True)

print(f"H3K27ac    REF(WT)={h3k_ref:.2f}    {dict(sorted(h3k_curve.items()))}")
print(f"H3K4me1    REF(WT)={h3k4me1_ref:.2f}    {dict(sorted(h3k4me1_curve.items()))}")
print(f"ATAC       REF(WT)={atac_ref:.4f}    {dict(sorted(atac_curve.items()))}")
print(f"ATG12 RNA  REF(WT)={rna_ref:.4f}    {dict(sorted(rna_curve.items()))}")

# ── Three stacked subplots, no twin axes, no broken axis ───────────────────
fig, axes = plt.subplots(
    3, 1, figsize=(args.width, args.height), sharex=True,
    gridspec_kw={'height_ratios': [1.5, 1.0, 1.0], 'hspace': 0.15},
)
ax_hist, ax_atac, ax_rna = axes

MARKERS = {'H3K27ac': 'o', 'H3K4me1': 'X', 'ATAC': 'v', 'ATG12 RNA': '^'}

def plot_curve(ax, label, curve):
    xs = sorted(curve.keys()); ys = [curve[x] for x in xs]
    ax.plot(xs, ys, marker=MARKERS[label], color=PALETTE[label],
            linewidth=LW, markersize=MS,
            markeredgecolor='white', markeredgewidth=0.5,
            label=f'Predicted {label}')

# Subplot 1: H3K27ac + H3K4me1 (shared y-axis)
plot_curve(ax_hist, 'H3K27ac', h3k_curve)
plot_curve(ax_hist, 'H3K4me1', h3k4me1_curve)
# Subplot 2: ATAC (own y-axis)
plot_curve(ax_atac, 'ATAC', atac_curve)
# Subplot 3: ATG12 RNA (own y-axis, tiny range)
plot_curve(ax_rna, 'ATG12 RNA', rna_curve)

# Tighten ylims per subplot
hist_max = max(max(h3k_curve.values()), max(h3k4me1_curve.values()))
ax_hist.set_ylim(0, hist_max * 1.10)
atac_min, atac_max = min(atac_curve.values()), max(atac_curve.values())
ax_atac.set_ylim(atac_min - (atac_max - atac_min) * 0.10, atac_max * 1.10)
rna_min, rna_max = min(rna_curve.values()), max(rna_curve.values())
pad = (rna_max - rna_min) * 0.15
ax_rna.set_ylim(rna_min - pad, rna_max + pad)

# y-axis labels
ax_hist.set_ylabel('Predicted histone ChIP\npeak signal', fontsize=FS_AXIS)
ax_atac.set_ylabel('Predicted ATAC\npeak signal', fontsize=FS_AXIS)
ax_rna.set_ylabel('Predicted ATG12\nRNA-seq signal', fontsize=FS_AXIS)

# x-axis on the bottom only
ax_rna.set_xlabel('Intact AP1 motifs in LTR10.ATG12', fontsize=FS_AXIS)
ax_rna.set_xticks([0, 20, 40, 60, 80, 100, 120, 140, 160, 180, 200])
ax_rna.set_xlim(-5, 210)

# Clean spines + ticks
for ax in axes:
    ax.spines['top'].set_visible(False); ax.spines['right'].set_visible(False)
    ax.tick_params(axis='both', labelsize=FS_TICK)

# Per-subplot legends (each subplot has its own series, no shared legend needed)
ax_hist.legend(fontsize=FS_LEGEND, loc='upper right', frameon=False)
ax_atac.legend(fontsize=FS_LEGEND, loc='upper right', frameon=False)
ax_rna.legend(fontsize=FS_LEGEND, loc='upper right', frameon=False)

# WT marker on all 3 subplots
for ax in axes:
    ax.axvline(20, color='grey', linestyle=':', linewidth=0.8, alpha=0.7)
ax_hist.text(19, ax_hist.get_ylim()[1] * 0.97, 'wild-type\n(20 motifs)',
             fontsize=FS_WT, color='grey', ha='right', va='top')

plt.tight_layout()
out = Path(args.output); out.parent.mkdir(parents=True, exist_ok=True)
plt.savefig(out, dpi=300, bbox_inches='tight')
print(f"Wrote {out}")
