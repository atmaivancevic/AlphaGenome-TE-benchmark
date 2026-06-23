"""
Fig 1 Panel B/C browser shot: experimental ENCODE H3K27ac vs AlphaGenome-
predicted H3K27ac at a single biosample/locus, with a GENCODE gene track.

Three stacked tracks sharing the genomic x-axis (top -> bottom):
  1. ENCODE H3K27ac        — fold-change bigWig signal (dark red #8B1A1F)
  2. AlphaGenome H3K27ac   — AG-predicted signal, same biosample (light red)
  3. Gene track            — GENCODE v46 protein-coding, longest transcript/gene

The outlier peak span is highlighted with an orange band across all tracks
(Fig 4/5 convention). The two signal tracks keep INDEPENDENT y-axes — the
ENCODE fold-change scale and the AG output scale are not comparable, so the
story is the *shape* contrast (does one signal sit where the other is absent),
not absolute height. Each signal axis prints a single top tick = its own max.

This is a standalone matplotlib renderer (no AG-SDK plot_components), so the
locked assay palette is applied directly — no pikepdf recolor pass needed.

Usage:
  python scripts/fig1_baseline_chromatin/plot_panel_BC_browser.py \
      --label NFKBIA --biosample GM12878 --ontology EFO:0002784 \
      --bigwig data/encode_h3k27ac/bigwigs/human/ENCSR000AKC_ENCFF469WVA.bigWig \
      --coords chr14:35395000-35415000 --peak 35404062-35406008 \
      --output figures/FIG1_FINAL/Fig1BC_browser/NFKBIA_GM12878.pdf
"""

import os, sys, argparse, re
import numpy as np
import pandas as pd
import pyBigWig
import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle

mpl.rcParams['font.family']     = 'sans-serif'
mpl.rcParams['font.sans-serif'] = ['Helvetica', 'Arial', 'DejaVu Sans']
mpl.rcParams['pdf.fonttype']    = 42
BOLD_FONT = 'Arial'   # Mac Helvetica.ttc exposes no bold to matplotlib; Arial Bold substitutes

# Locked assay palette (reference_assay_color_palette). Experimental H3K27ac
# is the canonical dark red; the AG prediction gets a lighter tint of the same
# hue so the pair reads as "same mark, measured vs predicted".
EXP_COLOR = '#8B1A1F'
AG_COLOR  = '#C77478'
EXON_COLOR = '#000000'

p = argparse.ArgumentParser(description=__doc__,
                            formatter_class=argparse.RawDescriptionHelpFormatter)
p.add_argument('--label', required=True, help='Locus label for title + filename (e.g. NFKBIA)')
p.add_argument('--biosample', required=True, help='Biosample display name (e.g. GM12878)')
p.add_argument('--ontology', required=True, help='AG ontology curie (e.g. EFO:0002784)')
p.add_argument('--bigwig', required=True, help='ENCODE H3K27ac fold-change bigWig path')
p.add_argument('--coords', required=True, help='Visible window as chr:start-end')
p.add_argument('--peak', required=True, help='Outlier peak span to highlight as start-end')
p.add_argument('--output', required=True)
p.add_argument('--gene-feather', default='data/gencode.v46.annotation.feather')
p.add_argument('--exp-ymax', type=float, default=None, help='Override ENCODE track y-max')
p.add_argument('--ag-ymax',  type=float, default=None, help='Override AG track y-max')
p.add_argument('--scale-to-peaks', nargs=2, default=None,
               metavar=('EXP_NARROWPEAK', 'PRED_NARROWPEAK'),
               help='Fix the ENCODE and AG y-axis maxima to the genome-wide max '
                    'signalValue (col 7) of these matched experimental + predicted '
                    'narrowPeak files for this biosample — i.e. the SAME axis range '
                    'as the Panel A scatter. Track height then reads as "how big is '
                    'this peak relative to the whole cell line", which is what makes '
                    'over- vs under-prediction legible. Explicit --exp-ymax / '
                    '--ag-ymax take precedence over this.')
p.add_argument('--fig-width', type=float, default=7.0)
p.add_argument('--no-title', action='store_true')
p.add_argument('--title', default=None,
               help='Override the panel title text (default "<label> (<biosample>)").')
p.add_argument('--width-pt', type=float, default=None,
               help='Exact output width in points. With --height-pt, forces the '
                    'PDF to that exact size (fixed margins, no tight-bbox crop) so '
                    'it drops into a figure layout slot. Use with --no-title.')
p.add_argument('--height-pt', type=float, default=None,
               help='Exact output height in points (see --width-pt).')
args = p.parse_args()
EXACT = args.width_pt is not None and args.height_pt is not None

m = re.match(r'^(chr\S+):(\d+)-(\d+)$', args.coords)
if not m:
    sys.exit(f'ERROR: --coords must be chr:start-end (got {args.coords!r})')
chrom, win_start, win_end = m.group(1), int(m.group(2)), int(m.group(3))
m2 = re.match(r'^(\d+)-(\d+)$', args.peak)
if not m2:
    sys.exit(f'ERROR: --peak must be start-end (got {args.peak!r})')
peak_start, peak_end = int(m2.group(1)), int(m2.group(2))

os.makedirs(os.path.dirname(args.output), exist_ok=True)
SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
REPO = os.path.dirname(os.path.dirname(SCRIPT_DIR))

# Optional: fix both y-axis caps to the per-biosample scatter range (max col-7
# signalValue across all matched peaks). Explicit --exp-ymax/--ag-ymax win.
if args.scale_to_peaks:
    import gzip
    def _max_signalvalue(fp):
        op = gzip.open if fp.endswith('.gz') else open
        m = 0.0
        with op(fp, 'rt') as f:
            for line in f:
                c = line.rstrip('\n').split('\t')
                if len(c) >= 7:
                    m = max(m, float(c[6]))
        return m
    exp_cap = _max_signalvalue(args.scale_to_peaks[0])
    ag_cap  = _max_signalvalue(args.scale_to_peaks[1])
    if args.exp_ymax is None: args.exp_ymax = exp_cap
    if args.ag_ymax  is None: args.ag_ymax  = ag_cap
    print(f'Cell-line scatter range: ENCODE ymax={exp_cap:.1f}, AG ymax={ag_cap:.1f}')

# ── AlphaGenome prediction (1 Mb context, single biosample H3K27ac) ─────────
from alphagenome.data import genome
from alphagenome.models import dna_client
key_path = os.path.join(REPO, 'scripts', 'my_api_key.txt')
with open(key_path) as f:
    model = dna_client.create(f.read().strip())

interval = genome.Interval(chromosome=chrom, start=win_start, end=win_end, strand='+')
print(f'Predicting AG H3K27ac at {args.coords} ({args.ontology}) ...')
pred = model.predict_interval(
    interval=interval.resize(dna_client.SEQUENCE_LENGTH_1MB),
    requested_outputs={dna_client.OutputType.CHIP_HISTONE},
    ontology_terms=[args.ontology],
)
hist = pred.chip_histone
meta = hist.metadata
keep = [i for i in meta.index
        if meta.at[i, 'ontology_curie'] == args.ontology
        and meta.at[i, 'histone_mark'] == 'H3K27ac']
if len(keep) != 1:
    sys.exit(f'ERROR: expected 1 H3K27ac track for {args.ontology}, got {len(keep)}. '
             f'Marks available: {sorted(set(meta["histone_mark"]))}')
ag_td = hist.select_tracks_by_index(keep)
ag_vals = np.asarray(ag_td.values)[:, 0]
ag_iv = ag_td.interval
ag_res = int(round((ag_iv.end - ag_iv.start) / len(ag_vals)))
bin_starts = ag_iv.start + np.arange(len(ag_vals)) * ag_res
sel = np.where((bin_starts + ag_res > win_start) & (bin_starts < win_end))[0]
xs = bin_starts[sel].astype(int)        # left edge of each visible AG bin
ag_y = ag_vals[sel]
xc = xs + ag_res / 2.0                   # bin centers (genomic coords)

# ── ENCODE bigWig at the IDENTICAL bins (per-bp mean per AG bin) ─────────────
bw = pyBigWig.open(args.bigwig)
read0, read1 = int(xs[0]), int(xs[-1] + ag_res)
per_bp = np.array(bw.values(chrom, read0, read1), dtype=float)
bw.close()
per_bp = per_bp.reshape(len(sel), ag_res)
enc_y = np.nan_to_num(np.nanmean(per_bp, axis=1))
print(f'  {len(sel)} bins @ {ag_res} bp; ENCODE max={enc_y.max():.1f}  AG max={ag_y.max():.1f}')

# ── Gene models: longest protein-coding transcript per gene in the window ───
gdf = pd.read_feather(args.gene_feather)
ex = gdf[(gdf['Chromosome'] == chrom) & (gdf['gene_type'] == 'protein_coding')
         & (gdf['Feature'] == 'exon')
         & (gdf['Start'] < win_end) & (gdf['End'] > win_start)]
genes = []  # (name, strand, tx_start, tx_end, [(exon_start, exon_end), ...])
for gname, gg in ex.groupby('gene_name'):
    # pick the transcript with the widest genomic span (canonical-ish)
    tx_span = gg.groupby('transcript_id').apply(
        lambda d: d['End'].max() - d['Start'].min(), include_groups=False)
    best_tx = tx_span.idxmax()
    t = gg[gg['transcript_id'] == best_tx]
    exons = sorted(zip(t['Start'].tolist(), t['End'].tolist()))
    genes.append((gname, t['Strand'].iloc[0], int(t['Start'].min()), int(t['End'].max()), exons))

# Greedy row-packing so overlapping genes don't collide on one line.
genes.sort(key=lambda g: g[2])
rows = []  # each row = rightmost end so far
gene_row = {}
for g in genes:
    placed = False
    for ri, rend in enumerate(rows):
        if g[2] > rend + (win_end - win_start) * 0.02:   # 2% gap for the label
            rows[ri] = g[3]; gene_row[g[0]] = ri; placed = True; break
    if not placed:
        rows.append(g[3]); gene_row[g[0]] = len(rows) - 1
n_gene_rows = max(1, len(rows))

# ── Figure: 3 stacked panels, heights weighted (signal tracks taller) ───────
gene_ratio = 0.16 + 0.24 * n_gene_rows
hspace = 0.18
if EXACT:
    figsize = (args.width_pt / 72.0, args.height_pt / 72.0)
    # Figure-panel proportions: taller gene annotation (also opens up the
    # gene-model -> chromosome-axis gap, which lives inside the gene panel)
    # with only modest inter-panel spacing so the gene model sits close under
    # the signal tracks.
    gene_ratio = 0.45 + 0.32 * n_gene_rows
    hspace = 0.22
else:
    figsize = (args.fig_width, 1.3 + 1.0 + gene_ratio)
fig, (ax_e, ax_a, ax_g) = plt.subplots(
    3, 1, sharex=True, figsize=figsize,
    gridspec_kw={'height_ratios': [1.0, 1.0, gene_ratio], 'hspace': hspace})

def _signal(ax, y, color, ymax, label):
    ax.fill_between(xc, 0, y, color=color, linewidth=0, zorder=3)
    hi = ymax if ymax is not None else (float(np.max(y)) * 1.05 if np.max(y) > 0 else 1.0)
    ax.set_ylim(0, hi)
    ax.set_yticks([hi]); ax.set_yticklabels([f'{hi:.0f}' if hi >= 10 else f'{hi:.1f}'])
    ax.tick_params(axis='y', length=0, labelsize=7, pad=2)
    for s in ('top', 'right', 'bottom'):
        ax.spines[s].set_visible(False)
    ax.set_ylabel(label, rotation=0, ha='right', va='center', fontsize=9)
    ax.yaxis.set_label_coords(-0.015, 0.5)

_signal(ax_e, enc_y, EXP_COLOR, args.exp_ymax, 'ENCODE\nH3K27ac')
_signal(ax_a, ag_y,  AG_COLOR,  args.ag_ymax,  'AG\nH3K27ac')

# Gene track. Rows stack top-down (row 0 at top); each gene labelled just
# below its model. ylim inverted (larger value first) so row 0 sits at top.
# Small top padding (-0.3) keeps the gene model close under the signal tracks;
# the larger bottom padding leaves the gene-model -> chromosome-axis gap.
ax_g.set_ylim(n_gene_rows - 0.55, -0.18)
for s in ('top', 'right', 'left'):
    ax_g.spines[s].set_visible(False)
ax_g.set_yticks([])
ARROW_SPACING = (win_end - win_start) / 22.0   # ~22 arrowhead slots across the window
for gname, strand, ts, te, exons in genes:
    row = gene_row[gname]
    gx0, gx1 = max(ts, win_start), min(te, win_end)
    ax_g.plot([gx0, gx1], [row, row], color=EXON_COLOR, linewidth=0.8, zorder=2)
    for es, ee in exons:
        ax_g.add_patch(Rectangle((es, row - 0.10), ee - es, 0.20,
                                 facecolor=EXON_COLOR, edgecolor='none', zorder=3))
    # Repeated arrowheads along the intron line, pointing 5'->3' (UCSC style):
    # '>' on + strand, '<' on - strand. Placed on intron segments only (the
    # gaps between in-window exons) so they don't sit on top of exon blocks.
    amark = '>' if strand == '+' else '<'
    exs = sorted((es, ee) for es, ee in exons)
    introns, prev = [], gx0
    for es, ee in exs:
        if es > prev:
            introns.append((prev, min(es, gx1)))
        prev = max(prev, ee)
        if prev >= gx1:
            break
    if prev < gx1:
        introns.append((prev, gx1))
    for s0, s1 in introns:
        if s1 - s0 < ARROW_SPACING * 0.4:    # skip introns too short to hold one
            continue
        n = max(1, int(round((s1 - s0) / ARROW_SPACING)))
        for k in range(n):
            xa = s0 + (s1 - s0) * (k + 0.5) / n
            ax_g.plot(xa, row, marker=amark, color=EXON_COLOR,
                      markersize=3.2, markeredgewidth=0, zorder=3, clip_on=True)
    mid = (gx0 + gx1) / 2
    ax_g.text(mid, row + 0.22, gname, ha='center', va='top',
              fontsize=7, style='italic', clip_on=True)

# Orange highlight: drawn after layout as ONE figure-level band (see
# draw_highlight) spanning the full vertical extent of all three panels so it
# reads as a single continuous vertical stripe, bridging the inter-panel gaps.
def draw_biosample_header():
    # Bold biosample name above the top track, right-aligned to the same
    # x-anchor as the track labels (matches the Fig 4/5 LTR10 browser shots,
    # where the biosample group name sits as a bold header above its tracks).
    pos = ax_e.get_position()
    x = pos.x0 + (-0.015) * pos.width          # = the track-label (ylabel) anchor
    fig.text(x, pos.y1 + 0.034, args.biosample, ha='right', va='bottom',
             fontweight='bold', family=BOLD_FONT, fontsize=9)

def draw_highlight():
    fig.canvas.draw()   # finalize axes positions + transforms
    inv = fig.transFigure.inverted()
    x0 = inv.transform(ax_e.transData.transform((peak_start, 0)))[0]
    x1 = inv.transform(ax_e.transData.transform((peak_end, 0)))[0]
    y1 = ax_e.get_position().y1
    y0 = ax_g.get_position().y0
    fig.add_artist(Rectangle((x0, y0), x1 - x0, y1 - y0,
                             transform=fig.transFigure, facecolor='orange',
                             alpha=0.18, edgecolor='none', zorder=10, clip_on=False))

# x-axis: only the bottom panel, labelled in Mb. Round-number ticks pruned
# from both ends so (a) labels read as clean Mb values and (b) no tick sits
# on the exact axis edge (which would clip its label against the fixed canvas).
from matplotlib.ticker import MaxNLocator, FuncFormatter
ax_g.set_xlim(win_start, win_end)
ax_g.xaxis.set_major_locator(MaxNLocator(nbins=5, steps=[1, 2, 2.5, 5, 10], prune='both'))
ax_g.xaxis.set_major_formatter(FuncFormatter(lambda x, _: f'{x/1e6:.3f}'))
ax_g.tick_params(axis='x', labelsize=7, length=2, pad=2)
ax_g.set_xlabel(f'{chrom} position (Mb)', fontsize=9, labelpad=4)

if not args.no_title:
    title_text = args.title if args.title else f'{args.label}  ({args.biosample})'
    fig.suptitle(title_text,
                 fontsize=8 if EXACT else 10,
                 fontweight='bold', family=BOLD_FONT,
                 y=0.965 if EXACT else 0.995)

for ax in (ax_e, ax_a):
    for c in ax.collections:
        c.set_rasterized(True)

if EXACT:
    # Fixed margins (fractions of the exact canvas) so content sits correctly:
    # left = track labels, bottom = x-axis label, top = title (reclaimed when
    # --no-title). Save with an explicit fixed Bbox = the full canvas so the
    # PDF page is EXACTLY width-pt x height-pt — without it the mixed
    # vector/raster PDF renderer snaps the page to the dpi pixel grid and the
    # size drifts by ~0.5 pt.
    from matplotlib.transforms import Bbox
    # Reserve top space for both the centred title (very top) and the bold
    # biosample header that sits just above the top track.
    top = 0.965 if args.no_title else 0.80
    fig.subplots_adjust(left=0.205, right=0.985, top=top, bottom=0.16, hspace=hspace)
    draw_biosample_header()
    draw_highlight()
    fig.savefig(args.output, dpi=200,
                bbox_inches=Bbox([[0, 0], [args.width_pt / 72.0, args.height_pt / 72.0]]))
else:
    draw_biosample_header()
    draw_highlight()
    fig.savefig(args.output, dpi=200, bbox_inches='tight')
plt.close()
print(f'Wrote {args.output}')
