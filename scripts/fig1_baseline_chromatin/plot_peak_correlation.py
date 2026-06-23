"""
Compare experimental vs predicted ChIP-seq peak signal.

Takes two narrowPeak files (experimental and predicted) with identical
peak coordinates, and produces a scatter plot of signalValue (column 7)
with Pearson and Spearman correlations.

Two render modes:
  default      — standalone 6x6 in panel with axis titles + boxed stats.
  --figure-mode — small square panel sized for a multi-panel figure grid:
                  no axis titles (add shared labels in the layout tool),
                  small fonts, and fixed internal margins so every panel
                  tiles and aligns. Panel side set by --panel-pt.

Usage:
    python scripts/fig1_baseline_chromatin/plot_peak_correlation.py \
        --experimental data/encode_h3k27ac/peaks/merged/EFO_0002824_HCT116.narrowPeak.gz \
        --predicted results/AG_predicted_h3k27ac_batched/EFO_0002824_HCT116/EFO_0002824_HCT116.narrowPeak.gz \
        --label HCT116 \
        --output figures/HCT116_exp_vs_pred.pdf

    # Fig 1A grid panel (111 pt square, no axis titles):
    python scripts/fig1_baseline_chromatin/plot_peak_correlation.py ... --figure-mode --panel-pt 111
"""

import os, argparse, gzip
import numpy as np
import matplotlib
import matplotlib.pyplot as plt
from matplotlib.ticker import MaxNLocator
from scipy import stats

matplotlib.rcParams['font.family'] = 'Helvetica'
matplotlib.rcParams['pdf.fonttype'] = 42
matplotlib.rcParams['ps.fonttype'] = 42

SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
# Script lives in scripts/fig1_baseline_chromatin/, so repo root is two levels up.
PROJECT_ROOT = os.path.dirname(os.path.dirname(SCRIPT_DIR))
FIGURES_DIR = os.path.join(PROJECT_ROOT, 'figures')
os.makedirs(FIGURES_DIR, exist_ok=True)

parser = argparse.ArgumentParser()
parser.add_argument('--experimental', required=True,
                    help='Experimental narrowPeak file (gzipped ok)')
parser.add_argument('--predicted', required=True,
                    help='Predicted narrowPeak file (gzipped ok)')
parser.add_argument('--label', required=True,
                    help='Cell type label for title (e.g. HCT116)')
parser.add_argument('--output', default=None,
                    help='Output PDF path')
parser.add_argument('--width', type=float, default=6.0, help='Figure width in inches (default mode)')
parser.add_argument('--height', type=float, default=6.0, help='Figure height in inches (default mode)')
parser.add_argument('--figure-mode', action='store_true',
                    help='Render a small square grid panel: no axis titles, '
                         'small fonts, fixed margins so panels tile/align.')
parser.add_argument('--panel-pt', type=float, default=125.0,
                    help='Panel side length in points (figure-mode only)')
parser.add_argument('--highlight-tsv', default=None,
                    help='TSV of outlier candidates from find_panel_outliers.py '
                         'OR the annotated supp_table_fig1_BC_outliers.tsv. Rows '
                         'matching --label by `biosample` are drawn as larger red '
                         '(over-pred) / blue (under-pred) dots over the base '
                         'scatter. The `direction` column may be the raw '
                         'B_over/C_under keys or the human-readable '
                         '"AG over-/under-predicts ..." strings.')
parser.add_argument('--label-genes', nargs='*', default=None,
                    help='Gene name(s) to label on the plot (e.g. NFKBIA PCDH7). '
                         'Matched against the highlight TSV `nearest_gene` column '
                         '(requires the annotated supp table). Each matching dot '
                         'gets an italic gene-name label with a thin leader line.')
args = parser.parse_args()

# ── Load peaks ──────────────────────────────────────────────────────────────

def load_narrowpeak(path):
    opener = gzip.open if path.endswith('.gz') else open
    coords, signals = [], []
    with opener(path, 'rt') as f:
        for line in f:
            if line.startswith('#') or line.strip() == '':
                continue
            cols = line.strip().split('\t')
            coords.append((cols[0], int(cols[1]), int(cols[2])))
            signals.append(float(cols[6]))
    return coords, np.array(signals)

exp_coords, exp_sig = load_narrowpeak(args.experimental)
pred_coords, pred_sig = load_narrowpeak(args.predicted)

# ── Match peaks by coordinate ──────────────────────────────────────────────

exp_dict = {coord: sig for coord, sig in zip(exp_coords, exp_sig)}
pred_dict = {coord: sig for coord, sig in zip(pred_coords, pred_sig)}

common = sorted(set(exp_dict) & set(pred_dict))
n_exp_only = len(exp_dict) - len(common)
n_pred_only = len(pred_dict) - len(common)

exp_sig = np.array([exp_dict[c] for c in common])
pred_sig = np.array([pred_dict[c] for c in common])

print(f"Matched {len(common)} peaks for {args.label}")
if n_exp_only > 0:
    print(f"  {n_exp_only} peaks in experimental only (skipped)")
if n_pred_only > 0:
    print(f"  {n_pred_only} peaks in predicted only (skipped)")
print(f"  Experimental signal: min={exp_sig.min():.1f}, max={exp_sig.max():.1f}, median={np.median(exp_sig):.1f}")
print(f"  Predicted signal:    min={pred_sig.min():.1f}, max={pred_sig.max():.1f}, median={np.median(pred_sig):.1f}")

# ── Correlations ────────────────────────────────────────────────────────────

r_pearson, p_pearson = stats.pearsonr(exp_sig, pred_sig)
r_spearman, p_spearman = stats.spearmanr(exp_sig, pred_sig)

print(f"\n  Pearson r  = {r_pearson:.3f}  (p = {p_pearson:.2e})")
print(f"  Spearman ρ = {r_spearman:.3f}  (p = {p_spearman:.2e})")

DOT_COLOR = '#E67E22'   # orange — genome-wide baseline benchmark
HIGHLIGHT_COLORS = {'B_over': '#D62728', 'C_under': '#1F4E79'}  # red / dark blue

def _norm_direction(val):
    """Accept either the raw 'B_over'/'C_under' keys (find_panel_outliers
    output) or the human-readable 'AG over-/under-predicts ...' strings
    (annotated supp table) so either file works as --highlight-tsv."""
    if val in HIGHLIGHT_COLORS:
        return val
    v = val.lower()
    if 'over' in v:
        return 'B_over'
    if 'under' in v:
        return 'C_under'
    return None

# ── Load highlight set ──────────────────────────────────────────────────────

label_genes = set(args.label_genes or [])
highlight = {'B_over': [], 'C_under': []}
gene_labels = []   # (exp, pred, gene_name, direction_key) for --label-genes
if args.highlight_tsv:
    import csv
    with open(args.highlight_tsv) as f:
        reader = csv.DictReader(f, delimiter='\t')
        for row in reader:
            if row['biosample'] != args.label:
                continue
            d = _norm_direction(row['direction'])
            if d is None:
                continue
            coord = (row['chrom'], int(row['start']), int(row['end']))
            e = exp_dict.get(coord)
            p = pred_dict.get(coord)
            if e is None or p is None:
                continue
            highlight[d].append((e, p))
            gene = (row.get('nearest_gene') or '').strip()
            if gene and gene in label_genes:
                gene_labels.append((e, p, gene, d))
    n_b, n_c = len(highlight['B_over']), len(highlight['C_under'])
    print(f"  Highlighting {n_b} B (over-pred) + {n_c} C (under-pred) outliers")
    if gene_labels:
        print(f"  Labelling genes: {', '.join(sorted({g for _, _, g, _ in gene_labels}))}")

def annotate_genes(ax, fontsize):
    """Label each --label-genes dot unambiguously: ring the specific dot in
    black (so it stands out from its same-colored highlight neighbours), then
    point an arrow at it from a gene-name label offset into nearby clear space
    (over-pred dots → up-left; under-pred dots → up). White bbox keeps the
    label legible over the point cloud."""
    if not gene_labels:
        return
    # Per-direction label placement (kept inside the axes):
    #   over-pred (red, upper-left) → up-and-right into the plot interior;
    #   under-pred (blue, lower-right) → directly below the dot, short arrow.
    PLACE = {
        'B_over':  dict(dx=-10, dy=28,  ha='right',  va='bottom'),
        'C_under': dict(dx=0,   dy=-22, ha='center', va='top'),
    }
    for e, p, gene, d in gene_labels:
        ax.scatter([e], [p], s=170, facecolors='none', edgecolors='black',
                   linewidths=1.3, zorder=7)
        pl = PLACE[d]
        ax.annotate(
            gene, xy=(e, p), xytext=(pl['dx'], pl['dy']), textcoords='offset points',
            ha=pl['ha'], va=pl['va'], fontsize=fontsize, fontstyle='italic',
            color='#111111', zorder=8,
            bbox=dict(boxstyle='round,pad=0.2', fc='white', ec='none', alpha=0.85),
            arrowprops=dict(arrowstyle='->', color='black', lw=0.9,
                            shrinkA=2, shrinkB=4))

# ── Plot ────────────────────────────────────────────────────────────────────

if args.figure_mode:
    # Small square panel for a multi-panel grid. Fixed add_axes rectangle +
    # no tight bbox => every panel PDF is exactly panel_pt square with the
    # data box at identical coordinates, so panels tile and align.
    side_in = args.panel_pt / 72.0
    fig = plt.figure(figsize=(side_in, side_in))
    ax = fig.add_axes([0.21, 0.16, 0.75, 0.68])   # [left, bottom, w, h]

    FS_TITLE, FS_TICK, FS_STATS, FS_N = 10.0, 8.0, 8.0, 8.0
    MARKER_SIZE, LINE_W = 2.0, 1.0

    ax.scatter(exp_sig, pred_sig, s=MARKER_SIZE, alpha=0.35, color=DOT_COLOR,
               edgecolors='none', zorder=3, rasterized=True)

    slope, intercept = np.polyfit(exp_sig, pred_sig, 1)
    x_line = np.linspace(exp_sig.min(), exp_sig.max(), 100)
    ax.plot(x_line, slope * x_line + intercept, '--', color='#333333',
            linewidth=LINE_W, zorder=2)

    for direction, pts in highlight.items():
        if not pts:
            continue
        xs, ys = zip(*pts)
        ax.scatter(xs, ys, s=MARKER_SIZE * 6, color=HIGHLIGHT_COLORS[direction],
                   edgecolors='white', linewidths=0.4, zorder=5)

    annotate_genes(ax, FS_STATS)

    ax.text(0.05, 0.96,
            f'r = {r_pearson:.2f}\nρ = {r_spearman:.2f}',
            transform=ax.transAxes, fontsize=FS_STATS, va='top', linespacing=1.3)
    ax.text(0.95, 0.06, f'n = {len(exp_sig):,}',
            transform=ax.transAxes, fontsize=FS_N, ha='right', va='bottom')

    ax.set_title(args.label, fontsize=FS_TITLE, fontweight='bold', pad=3)
    ax.tick_params(axis='both', labelsize=FS_TICK, length=2, pad=1.5, width=0.6)
    ax.xaxis.set_major_locator(MaxNLocator(4))
    ax.yaxis.set_major_locator(MaxNLocator(4))
    for s in ('top', 'right'):
        ax.spines[s].set_visible(False)
    for s in ('left', 'bottom'):
        ax.spines[s].set_linewidth(0.6)

    output_path = args.output or os.path.join(FIGURES_DIR, f'{args.label}_exp_vs_pred.pdf')
    plt.savefig(output_path, dpi=600)   # no bbox_inches='tight' — keep fixed size

else:
    FS_TITLE, FS_AXIS, FS_TICK, FS_STATS, FS_N = 13, 11, 10, 10, 10
    MARKER_SIZE = 15

    fig, ax = plt.subplots(figsize=(args.width, args.height))

    ax.scatter(exp_sig, pred_sig, s=MARKER_SIZE, alpha=0.4, color=DOT_COLOR,
               edgecolors='none', zorder=3, rasterized=True)

    slope, intercept = np.polyfit(exp_sig, pred_sig, 1)
    x_line = np.linspace(exp_sig.min(), exp_sig.max(), 100)
    ax.plot(x_line, slope * x_line + intercept, '--', color='#333333',
            linewidth=1.5, zorder=2)

    for direction, pts in highlight.items():
        if not pts:
            continue
        xs, ys = zip(*pts)
        ax.scatter(xs, ys, s=MARKER_SIZE * 5, color=HIGHLIGHT_COLORS[direction],
                   edgecolors='white', linewidths=0.6, zorder=5)

    annotate_genes(ax, FS_STATS)

    ax.text(0.05, 0.95,
            f'Pearson r = {r_pearson:.2f}\n'
            f'Spearman ρ = {r_spearman:.2f}',
            transform=ax.transAxes, fontsize=FS_STATS, va='top',
            bbox=dict(boxstyle='round,pad=0.3', fc='white', ec='#CCCCCC', alpha=0.9))

    ax.text(0.95, 0.05, f'n = {len(exp_sig)}',
            transform=ax.transAxes, fontsize=FS_N, ha='right', va='bottom')

    ax.set_xlabel('Experimental H3K27ac signal (ENCODE)', fontsize=FS_AXIS)
    ax.set_ylabel('Predicted H3K27ac signal (AlphaGenome)', fontsize=FS_AXIS)
    ax.set_title(f'{args.label}', fontsize=FS_TITLE, fontweight='bold')
    ax.tick_params(axis='both', labelsize=FS_TICK)
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)

    plt.tight_layout()
    output_path = args.output or os.path.join(FIGURES_DIR, f'{args.label}_exp_vs_pred.pdf')
    plt.savefig(output_path, dpi=300, bbox_inches='tight')

print(f"\nSaved: {output_path}")
