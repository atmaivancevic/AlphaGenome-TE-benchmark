"""
Fig 6 Panel A1 — AP1 motif titration: AG-predicted TF binding (FOSL1 + JUND) on a
shared y-axis. Both climb near-linearly with motif count through 200 motifs (AG
predicts TF binding as motif-counting), unlike the saturating chromatin in Panel A2.
X = intact AP1 motifs (0-200, scramble + addition alleles); Y = raw AG ChIP signal.

Example usage:
python scripts/fig6_AP1_perturbation/plot_fig6_panelA1_TF_titration.py \
    --chromatin-csv \
        results/AG_perturbation_LTR10_ATG12_AP1/predict_chromatin_at_element.csv \
        results/AG_perturbation_LTR10_ATG12_AP1/predict_chromatin_at_element_INS.csv \
    --output figures/FIG6_FINAL/panelA1_TF_titration.pdf
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
ap.add_argument('--chromatin-csv', nargs='+', required=True,
                help='One or more chromatin prediction CSVs from predict_variant_tracks.py')
ap.add_argument('--output', required=True)
ap.add_argument('--width',  type=float, default=5.5, help='Figure width inches')
ap.add_argument('--height', type=float, default=3.8, help='Figure height inches')
args = ap.parse_args()

PALETTE = {'FOSL1': '#D7301F', 'JUND': '#F0A030'}
FS_AXIS, FS_TICK, FS_LEGEND, FS_WT = 10, 9, 9, 8
LW, MS = 2.0, 7

# Allele → intact AP1 motif count (scramble + addition titration)
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
    # add190-220 alleles (210-240 motifs) exist in the tab but are excluded;
    # x-axis capped at 200 (plateau already established).
}

def raw_curve(csv_paths, track_substring):
    """Return ({motif_count: alt_signal_mean}, ref_signal_mean) for the
    given track. WT (motif=20) is plotted at the REF prediction."""
    dfs = [pd.read_csv(p) for p in csv_paths]
    df = pd.concat(dfs, ignore_index=True)
    df = df[df['track_name'].str.contains(track_substring, case=False, na=False)]
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

fosl1_curve, fosl1_ref = raw_curve(args.chromatin_csv, 'FOSL1')
jund_curve,  jund_ref  = raw_curve(args.chromatin_csv, 'JUND')

print(f"FOSL1  REF(WT)={fosl1_ref:.2f}    {dict(sorted(fosl1_curve.items()))}")
print(f"JUND   REF(WT)={jund_ref:.2f}    {dict(sorted(jund_curve.items()))}")

# ── Plot — single subplot, shared y-axis ────────────────────────────────────
fig, ax = plt.subplots(figsize=(args.width, args.height))

MARKERS = {'FOSL1': 's', 'JUND': 'D'}

def plot_curve(ax, label, curve):
    xs = sorted(curve.keys())
    ys = [curve[x] for x in xs]
    ax.plot(xs, ys, marker=MARKERS[label], color=PALETTE[label],
            linewidth=LW, markersize=MS,
            markeredgecolor='white', markeredgewidth=0.5,
            label=f'Predicted {label}')

plot_curve(ax, 'FOSL1', fosl1_curve)
plot_curve(ax, 'JUND',  jund_curve)

# WT marker
ax.axvline(20, color='grey', linestyle=':', linewidth=0.8, alpha=0.7)
ax.text(19, ax.get_ylim()[1] * 0.95, 'wild-type\n(20 motifs)',
        fontsize=FS_WT, color='grey', ha='right', va='top')

ax.set_xlabel('Intact AP1 motifs in LTR10.ATG12', fontsize=FS_AXIS)
ax.set_ylabel('Predicted TF ChIP peak signal', fontsize=FS_AXIS)
ax.set_xticks([0, 20, 40, 60, 80, 100, 120, 140, 160, 180, 200])
ax.set_xlim(-5, 210)
ax.tick_params(axis='both', labelsize=FS_TICK)
ax.spines['top'].set_visible(False); ax.spines['right'].set_visible(False)
ax.legend(fontsize=FS_LEGEND, loc='upper left', frameon=False)

plt.tight_layout()
out = Path(args.output); out.parent.mkdir(parents=True, exist_ok=True)
plt.savefig(out, dpi=300, bbox_inches='tight')
print(f"Wrote {out}")
