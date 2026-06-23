"""
Generate an AlphaGenome SDK genome-browser-style PDF for one LTR10 enhancer.

Stacks tracks in this order (all HCT116 first; "Normal colon (sigmoid)"
group at the bottom — matches Atma's inspo layout):
  HCT116 group (per-track label = assay name only, no biosample prefix):
    JUND
    FOSL1
    H3K27ac
    H3K4me1
    ATAC-seq
    total RNA-seq + strand
    total RNA-seq − strand (axis INVERTED so peaks point DOWN,
                            butted against + strand at the 0-line)
  Normal colon (sigmoid) group:
    H3K27ac
    H3K4me1

The enhancer span (from supp 10 'pos' + 'SVLEN') is highlighted with two
orange vertical bars via VariantAnnotation. The resulting PDF still uses
the default AG SDK track colours; pipe through
scripts/fig4_5_LTR10_CRISPR_comparison/recolor_AG_screenshot_tracks.py to apply the locked assay palette.

Usage:
  python scripts/fig4_5_LTR10_CRISPR_comparison/plot_AG_browser_shot.py \
      --variant-id LTR10.ATG12 \
      --coords chr5:115820000-115940000 \
      --output figures/FIG5_FINAL/browser_shots/raw/LTR10.ATG12_browser.pdf
"""

import os, sys, argparse, re
import numpy as np
import pandas as pd
import matplotlib as mpl
import matplotlib.pyplot as plt

# Helvetica across all text — title, axis labels, ticklabels, the
# variant-id label drawn after plot. NOTE: Mac's bundled Helvetica.ttc
# only exposes the regular weight to matplotlib (no bold variant), so
# bold-requested text falls back to Arial Bold which is a near-identical
# substitute and DOES have a proper bold weight available.
mpl.rcParams['font.family']     = 'sans-serif'
mpl.rcParams['font.sans-serif'] = ['Helvetica', 'Arial', 'DejaVu Sans']
mpl.rcParams['pdf.fonttype']    = 42

BOLD_FONT = 'Arial'   # see note above on Mac Helvetica bold-fallback

p = argparse.ArgumentParser(description=__doc__,
                            formatter_class=argparse.RawDescriptionHelpFormatter)
p.add_argument('--variant-id', required=True,
               help='LTR10.XXX variant ID (used for highlight + title)')
p.add_argument('--coords', required=True,
               help='Visible window as chr:start-end')
p.add_argument('--enhancer-pos', default=None,
               help='Optional enhancer span as start-end. '
                    'If omitted, derived from supp 10 (pos, SVLEN).')
p.add_argument('--supp10',
               default='supptables/supp_table_10_LTR10_CRISPR_AG_predictions.tsv')
p.add_argument('--output', required=True)
p.add_argument('--hct116-ontology', default='EFO:0002824')
p.add_argument('--sigmoid-ontology', default='UBERON:0001159',
               help='Healthy sigmoid colon (default UBERON:0001159)')
p.add_argument('--fig-width', type=float, default=14)
p.add_argument('--track-height', type=float, default=0.5,
               help='Per-track height in inches (default 0.5 — flatter than AG SDK default)')
p.add_argument('--min-rna-neg-ymax', type=float, default=0.5,
               help='Minimum y-axis upper bound for the - strand RNA panel '
                    '(default 0.5; raise per-variant when faint - strand signal '
                    'still dominates the panel visually)')
p.add_argument('--h3k27ac-ymax', type=float, default=None,
               help='Explicit y-max for BOTH HCT116 and sigmoid colon H3K27ac '
                    'panels. Useful when sigmoid auto-max dwarfs the HCT116 '
                    'peak and you want to cap them at a smaller shared scale.')
p.add_argument('--h3k4me1-ymax', type=float, default=None,
               help='Explicit y-max for BOTH HCT116 and sigmoid colon H3K4me1 panels.')
p.add_argument('--tf-ymax', type=float, default=None,
               help='Explicit y-max for BOTH HCT116 TF ChIP-seq panels (JUND + FOSL1).')
p.add_argument('--atac-ymax', type=float, default=None,
               help='Explicit y-max for the HCT116 ATAC-seq panel.')
p.add_argument('--no-title', action='store_true',
               help='Suppress the figure title (the LTR10.XXX bold label '
                    'already identifies the panel; useful when figure sits '
                    'in a multi-panel layout).')
p.add_argument('--no-jund', action='store_true',
               help='Drop the HCT116 JUND track (keep just FOSL1 for AP1).')
p.add_argument('--no-atac', action='store_true',
               help='Drop the HCT116 ATAC-seq track.')
p.add_argument('--no-h3k4me1', action='store_true',
               help='Drop the H3K4me1 pair (both HCT116 and sigmoid colon).')
p.add_argument('--compact', action='store_true',
               help='Compact-mode layout (smaller fonts, bigger left margin, '
                    'wider label-anchor offset). For small multi-panel figures '
                    'like Fig 5 rows; default layout is sized for the Fig 4 '
                    'standalone ATG12 example.')
args = p.parse_args()

m = re.match(r'^(chr\S+):(\d+)-(\d+)$', args.coords)
if not m:
    sys.exit(f'ERROR: --coords must be chr:start-end (got {args.coords!r})')
chrom, view_start, view_end = m.group(1), int(m.group(2)), int(m.group(3))

if args.enhancer_pos:
    m2 = re.match(r'^(\d+)-(\d+)$', args.enhancer_pos)
    if not m2:
        sys.exit(f'ERROR: --enhancer-pos must be start-end (got {args.enhancer_pos!r})')
    enh_start, enh_end = int(m2.group(1)), int(m2.group(2))
else:
    s10 = pd.read_table(args.supp10, dtype=str, keep_default_na=False)
    s10['variant_id_ff'] = s10['variant_id'].replace('', pd.NA).ffill()
    sub = s10[s10['variant_id_ff'] == args.variant_id]
    if sub.empty:
        sys.exit(f'ERROR: {args.variant_id!r} not found in supp 10')
    pos = int(sub.iloc[0]['pos'])
    svlen = abs(int(sub.iloc[0]['SVLEN']))
    enh_start, enh_end = pos, pos + svlen
    print(f'Enhancer span from supp 10: {chrom}:{enh_start}-{enh_end} ({svlen} bp)')

os.makedirs(os.path.dirname(args.output), exist_ok=True)

SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
from alphagenome.data import genome, gene_annotation, track_data
from alphagenome.data import transcript as transcript_utils
from alphagenome.models import dna_client
from alphagenome.visualization import plot_components
from matplotlib.patches import Rectangle

with open(os.path.join(SCRIPT_DIR, 'my_api_key.txt')) as f:
    api_key = f.read().strip()
model = dna_client.create(api_key)
OutputType = dna_client.OutputType

print('Loading GENCODE v46 annotations...')
gtf = pd.read_feather(
    'https://storage.googleapis.com/alphagenome/reference/gencode/'
    'hg38/gencode.v46.annotation.gtf.gz.feather'
)
gtf_t = gene_annotation.filter_transcript_support_level(
    gene_annotation.filter_protein_coding(gtf), ['1']
)
gtf_l = gene_annotation.filter_to_longest_transcript(gtf_t)
tx_extractor = transcript_utils.TranscriptExtractor(gtf_l)

interval = genome.Interval(chromosome=chrom, start=view_start, end=view_end, strand='+')
print(f'Predicting at {chrom}:{view_start}-{view_end} (1 Mb context)...')
baseline = model.predict_interval(
    interval=interval.resize(dna_client.SEQUENCE_LENGTH_1MB),
    requested_outputs={
        OutputType.CHIP_TF, OutputType.CHIP_HISTONE,
        OutputType.ATAC, OutputType.RNA_SEQ,
    },
    ontology_terms=[args.hct116_ontology, args.sigmoid_ontology],
)

# Robust filter helper. Surfaces metadata columns / values when no track
# matches, so a config mismatch is immediately diagnosable.
def pick(td, label, **filters):
    meta = td.metadata
    keep = list(meta.index)
    for col, val in filters.items():
        if col not in meta.columns:
            sys.exit(f'ERROR ({label}): metadata column {col!r} not found. '
                     f'Columns: {list(meta.columns)}')
        if callable(val):
            keep = [i for i in keep if val(meta.at[i, col])]
        elif isinstance(val, (list, tuple, set)):
            keep = [i for i in keep if meta.at[i, col] in val]
        else:
            keep = [i for i in keep if meta.at[i, col] == val]
    if not keep:
        sys.exit(f'ERROR ({label}): no tracks matched {filters}. '
                 f'Available {list(filters)[0]} values: '
                 f'{sorted(set(meta[list(filters)[0]].astype(str)))[:30]}')
    return td.select_tracks_by_index(keep)

# Strip the "EFO:NNNNNN " / "UBERON:NNNNNN " prefix from the track
# 'name' metadata so the y-axis label reads e.g. "TF ChIP-seq JUND"
# instead of "EFO:0002824 TF ChIP-seq JUND". Compact labels let the
# figure adopt a flatter aspect ratio. Ontology IDs used in this
# script (recorded for provenance — these are the AG SDK queries):
#
#   EFO:0002824   HCT116 (colorectal carcinoma cell line)
#   UBERON:0001159 sigmoid colon (healthy normal-tissue reference)
#
def strip_ontology_prefix(td):
    md = td.metadata.copy()
    if 'name' in md.columns:
        md['name'] = md['name'].astype(str).str.replace(
            r'^(EFO|UBERON):\d+\s+', '', regex=True)
    return track_data.TrackData(
        values=td.values, metadata=md,
        resolution=td.resolution, interval=td.interval,
    )

def set_biosample(td, name):
    """Overwrite biosample_name in metadata so the ylabel reads however
    we want (e.g. 'Sigmoid colon' for the bottom-group tracks)."""
    md = td.metadata.copy()
    md['biosample_name'] = name
    return track_data.TrackData(
        values=td.values, metadata=md,
        resolution=td.resolution, interval=td.interval,
    )

def maybe_strip_assay_prefix(td):
    """In compact mode, drop the 'TF ChIP-seq ' / 'Histone ChIP-seq ' /
    'total ' prefixes so labels read as the bare assay/mark
    ('FOSL1', 'H3K27ac', 'RNA-seq'). Skipped in default mode."""
    if not args.compact:
        return td
    md = td.metadata.copy()
    if 'name' in md.columns:
        md['name'] = (md['name'].astype(str)
                      .str.replace(r'^(TF |Histone )?ChIP-seq\s+', '', regex=True)
                      .str.replace(r'^total\s+', '', regex=True))
    return track_data.TrackData(
        values=td.values, metadata=md,
        resolution=td.resolution, interval=td.interval,
    )

has = lambda sub: (lambda v: sub in str(v))
hct, sig = args.hct116_ontology, args.sigmoid_ontology

def _hct(td_picked): return maybe_strip_assay_prefix(strip_ontology_prefix(td_picked))
def _sig(td_picked): return maybe_strip_assay_prefix(set_biosample(strip_ontology_prefix(td_picked), 'Sigmoid colon'))

fosl1     = _hct(pick(baseline.chip_tf,      'HCT116 FOSL1',   ontology_curie=hct, transcription_factor='FOSL1'))
jund      = _hct(pick(baseline.chip_tf,      'HCT116 JUND',    ontology_curie=hct, transcription_factor='JUND'))
hct_27ac  = _hct(pick(baseline.chip_histone, 'HCT116 H3K27ac', ontology_curie=hct, histone_mark='H3K27ac'))
hct_4me1  = _hct(pick(baseline.chip_histone, 'HCT116 H3K4me1', ontology_curie=hct, histone_mark='H3K4me1'))
hct_atac  = _hct(pick(baseline.atac,         'HCT116 ATAC',    ontology_curie=hct))
hct_rna   = _hct(pick(baseline.rna_seq,      'HCT116 total RNA',
                      ontology_curie=hct, name=has('total RNA-seq'), strand=['+', '-']))
sig_27ac  = _sig(pick(baseline.chip_histone, 'sigmoid H3K27ac', ontology_curie=sig, histone_mark='H3K27ac'))
sig_4me1  = _sig(pick(baseline.chip_histone, 'sigmoid H3K4me1', ontology_curie=sig, histone_mark='H3K4me1'))

# Split stranded total-RNA into two single-track TrackDatas. Each gets
# its own panel with its own y-axis range — the small minus-strand
# signal stays visible at its own scale (vs the mirror approach which
# crushed it under the larger plus-strand range). The minus-strand
# panel is then post-rendered with an inverted y-axis so peaks point
# DOWN — standard genome-browser convention for stranded RNA.
strands = hct_rna.metadata['strand'].astype(str).values
plus_idx  = list(np.where(strands == '+')[0])
minus_idx = list(np.where(strands == '-')[0])
if len(plus_idx) != 1 or len(minus_idx) != 1:
    sys.exit(f'ERROR: expected one + and one - RNA track, '
             f'got {len(plus_idx)} / {len(minus_idx)}')
hct_rna_plus  = hct_rna.select_tracks_by_index(plus_idx)
hct_rna_minus = hct_rna.select_tracks_by_index(minus_idx)

transcripts = tx_extractor.extract(interval)
# Highlight the enhancer span with a single labelled orange box.
enh_interval = genome.Interval(chromosome=chrom, start=enh_start, end=enh_end, strand='+')

TH = args.track_height
T_HCT = lambda td: plot_components.Tracks(
    tdata=td, ylabel_template='{name}',
    filled=True, track_height=TH,
)
T_SIG = lambda td: plot_components.Tracks(
    tdata=td, ylabel_template='{name}',  # biosample name carried by the group header
    filled=True, track_height=TH,
)
# Build tracks list + index map dynamically. Indices are needed
# post-render for axis surgery (RNA invert + butt, y-axis matching,
# group headers, recolor). `idx` maps a stable name to the position
# in fig.get_axes() so flag-based track skipping doesn't break refs.
tracks = [plot_components.TranscriptAnnotation(
    transcripts,
    adaptive_fig_height=not args.compact,  # need this False for fig_height to apply
    fig_height=(0.5 if args.compact else 1.0),
)]
idx = {'transcript': 0}
def _add(name, td, factory):
    tracks.append(factory(td))
    idx[name] = len(tracks) - 1

if not args.no_jund:  _add('jund',  jund,  T_HCT)
_add('fosl1', fosl1, T_HCT)
_add('hct_h3k27ac', hct_27ac, T_HCT)
if not args.no_h3k4me1:  _add('hct_h3k4me1', hct_4me1, T_HCT)
if not args.no_atac:     _add('hct_atac',    hct_atac, T_HCT)
_add('rna_plus',  hct_rna_plus,  T_HCT)
_add('rna_minus', hct_rna_minus, T_HCT)
_add('sig_h3k27ac', sig_27ac, T_SIG)
if not args.no_h3k4me1:  _add('sig_h3k4me1', sig_4me1, T_SIG)

plot_components.plot(
    tracks,
    interval=interval,
    fig_width=args.fig_width,
    title=None if args.no_title else f'AlphaGenome predictions: {args.variant_id}',
)

# Layout constants. Compact mode tightens fonts + widens label anchors
# so multi-panel figures (Fig 5 rows) read at narrow fig-width without
# label/tick overlap; default sizing targets the Fig 4 standalone panel.
if args.compact:
    LABEL_X            = -0.06   # compact mode: tight against data while leaving breathing room before y-ticks
    HEADER_OFFSET      = 0.015
    HEADER_FONTSIZE    = 8
    GENE_LABEL_SIZE    = 8
    TRACK_LABEL_SIZE   = 8
    TICK_LABEL_SIZE    = 7
    LTR10_LABEL_SIZE   = 8
    LTR10_LABEL_OFFSET = 0.005
    GROUP_GAP          = 0.03
    FIG_LEFT_MARGIN    = 0.14   # y-ticks now live inside tracks (compact); axis can hug labels
    XLABEL_PAD         = 5
else:
    LABEL_X            = -0.06
    HEADER_OFFSET      = 0.008
    HEADER_FONTSIZE    = 11
    GENE_LABEL_SIZE    = 12
    TRACK_LABEL_SIZE   = 10  # matplotlib default
    TICK_LABEL_SIZE    = 8
    LTR10_LABEL_SIZE   = 12
    LTR10_LABEL_OFFSET = 0.003
    GROUP_GAP          = 0.03
    FIG_LEFT_MARGIN    = None
    XLABEL_PAD         = 7

fig = plt.gcf()

# In compact mode, push the left subplot edge in so labels have room.
if FIG_LEFT_MARGIN is not None:
    fig.subplots_adjust(left=FIG_LEFT_MARGIN)

# Identify RNA+/RNA- axes by creation-order index in fig.axes,
# looked up via the dynamic `idx` map so flag-driven track skipping
# (--no-jund, --no-atac, etc.) doesn't break the references.
all_axes = list(fig.get_axes())
if 'rna_minus' in idx and len(all_axes) > idx['rna_minus']:
    plus_ax  = all_axes[idx['rna_plus']]
    minus_ax = all_axes[idx['rna_minus']]
    # Enforce a minimum y-max on the - strand so very faint minus-strand
    # signal doesn't fill the whole panel and read as strong. Apply
    # BEFORE inverting (set_ylim takes data-axis order).
    cur_max = float(np.max(np.abs(minus_ax.get_ylim())))
    minus_ax.set_ylim(0, max(args.min_rna_neg_ymax, cur_max))
    minus_ax.invert_yaxis()
    plus_pos  = plus_ax.get_position()
    minus_pos = minus_ax.get_position()
    new_y1   = plus_pos.y0
    new_y0   = new_y1 - minus_pos.height
    minus_ax.set_position([minus_pos.x0, new_y0, minus_pos.width, minus_pos.height])
    # Single "total RNA-seq" label sitting at the shared 0-line (bottom
    # of + panel = top of - panel after butting them together). Same x
    # anchor as the other track labels for alignment.
    plus_ax.set_ylabel('')
    minus_ax.set_ylabel('')
    label_y = new_y1   # the 0-line in figure coords
    label_x = plus_pos.x0 + LABEL_X * plus_pos.width
    fig.text(label_x, label_y,
             'RNA-seq' if args.compact else 'Total RNA-seq',
             ha='right', va='center', fontsize=TRACK_LABEL_SIZE,
             transform=fig.transFigure)

# Match the HCT116 / sigmoid-colon y-axis ranges for each histone mark
# so the cancer-vs-normal comparison is on the same scale (paper
# convention). Use the explicit --h3k27ac-ymax / --h3k4me1-ymax if
# given; otherwise default to the larger of the two auto-maxes.
def _share_ylim(ax_a, ax_b, ymax_override=None):
    if ymax_override is not None:
        hi = ymax_override
    else:
        hi = max(ax_a.get_ylim()[1], ax_b.get_ylim()[1])
    lo = min(ax_a.get_ylim()[0], ax_b.get_ylim()[0], 0.0)
    ax_a.set_ylim(lo, hi)
    ax_b.set_ylim(lo, hi)
if 'hct_h3k27ac' in idx and 'sig_h3k27ac' in idx:
    _share_ylim(all_axes[idx['hct_h3k27ac']], all_axes[idx['sig_h3k27ac']], args.h3k27ac_ymax)
if 'hct_h3k4me1' in idx and 'sig_h3k4me1' in idx:
    _share_ylim(all_axes[idx['hct_h3k4me1']], all_axes[idx['sig_h3k4me1']], args.h3k4me1_ymax)
if 'jund' in idx and 'fosl1' in idx:
    _share_ylim(all_axes[idx['jund']], all_axes[idx['fosl1']], args.tf_ymax)
elif 'fosl1' in idx and args.tf_ymax is not None:
    all_axes[idx['fosl1']].set_ylim(0, args.tf_ymax)
if 'hct_atac' in idx and args.atac_ymax is not None:
    all_axes[idx['hct_atac']].set_ylim(0, args.atac_ymax)

# Insert a gap before the "Normal sigmoid colon" group: shift the two
# sigmoid axes down by GROUP_GAP (set above per compact/default mode).
sig_indices = [idx[k] for k in ('sig_h3k27ac', 'sig_h3k4me1') if k in idx]
for ax_idx in sig_indices:
    pos = all_axes[ax_idx].get_position()
    all_axes[ax_idx].set_position([pos.x0, pos.y0 - GROUP_GAP, pos.width, pos.height])

# Draw a single continuous orange highlight rectangle spanning every data
# axis (genome + tracks), in figure coordinates. Done after plot so it
# bridges the inter-axis gaps cleanly — no per-axis chunks, no doubled
# alpha at axis boundaries (the artifacts of IntervalAnnotation's
# per-axis ymax extension trick).
data_axes = [
    ax for ax in fig.get_axes()
    if abs(ax.get_xlim()[1] - ax.get_xlim()[0]) > 1000
    and ax.get_xlim()[0] >= view_start - 1
    and ax.get_xlim()[1] <= view_end + 1
]
if data_axes:
    inv = fig.transFigure.inverted()
    ref_ax = data_axes[0]
    x0 = inv.transform(ref_ax.transData.transform((enh_start, 0)))[0]
    x1 = inv.transform(ref_ax.transData.transform((enh_end, 0)))[0]
    # Skip the topmost axis (gene-transcript track) so the highlight
    # starts at the top of the first DATA track and the variant-id
    # label has room to sit cleanly in the gap above it.
    sorted_axes = sorted(data_axes, key=lambda a: a.get_position().y1, reverse=True)
    track_axes = sorted_axes[1:] if len(sorted_axes) > 1 else sorted_axes

    # Simplify each track's y-axis: keep only the top "nice" tick that
    # matplotlib auto-selected (e.g. 500, 1000, 30 — not the raw data
    # max which would print messy floats like 712 or 3.17188). No tick
    # marks. For the inverted - strand RNA axis the top tick is the
    # deepest peak value and visually sits at the bottom of the panel.
    def _format_tick(v):
        if abs(v) >= 10:
            return f'{int(round(v))}'
        if abs(v) >= 1:
            return f'{v:.1f}'
        return f'{v:.2f}'.rstrip('0').rstrip('.')

    for ax in track_axes:
        # Track-label font size (matches the compact/default mode).
        ax.yaxis.label.set_fontsize(TRACK_LABEL_SIZE)
        ylim = ax.get_ylim()
        is_inverted = ylim[0] > ylim[1]
        if is_inverted:
            # For the inverted - strand RNA panel, use the exact ylim max
            # (set explicitly via --min-rna-neg-ymax) so the displayed
            # tick is the user-controlled bound, e.g. "−3" for ymax=3.
            top_tick = max(ylim)
        else:
            # For normal axes, use matplotlib's auto-picked "nice" ticks.
            auto_ticks = [t for t in ax.get_yticks()
                          if min(ylim) <= t <= max(ylim) and t != 0]
            top_tick = max(auto_ticks, key=abs) if auto_ticks else max(ylim, key=abs)
        tick_text = _format_tick(-abs(top_tick) if is_inverted else top_tick)
        ax.set_yticks([top_tick])
        ax.set_yticklabels([tick_text])
        ax.tick_params(axis='y', length=0, labelsize=TICK_LABEL_SIZE,
                       pad=(1 if args.compact else 3))

    # Drop the numerical x-tick labels (e.g. "1.1582", "1.1594") and
    # the offset notation ("1e8") on every data axis — keep only the
    # genomic "Chromosome position; interval=..." xlabel at the bottom.
    for ax in track_axes + [sorted_axes[0]]:
        ax.set_xticks([])
        ax.xaxis.get_offset_text().set_visible(False)

    # Push the "Chromosome position;..." xlabel a bit further down so
    # it doesn't crowd the bottom track. labelpad is in points; default
    # is ~4. Find the axis that owns the xlabel (typically the
    # bottommost one) and bump its pad.
    for ax in sorted(track_axes, key=lambda a: a.get_position().y0):
        if ax.get_xlabel():
            ax.xaxis.labelpad = XLABEL_PAD
            ax.xaxis.label.set_fontsize(TRACK_LABEL_SIZE)
            break

    # Gene-name labels on the TranscriptAnnotation track.
    gene_ax = sorted_axes[0]
    for txt in gene_ax.texts:
        txt.set_fontsize(GENE_LABEL_SIZE)

    # Right-align all per-track y-axis labels to the common LABEL_X
    # anchor (defined globally above; matches the RNA "total RNA-seq"
    # label and the bold group headers).
    for ax in track_axes:
        ax.yaxis.label.set_horizontalalignment('right')
        ax.yaxis.set_label_coords(LABEL_X, 0.5)

    # Group headers (bold) above each block, right-aligned to the same
    # x-anchor as the track labels. Top HCT track is JUND if present,
    # else FOSL1. Sigmoid header sits above sigmoid H3K27ac.
    HEADER_OFFSET = 0.008
    def _header_x(ax_pos):
        return ax_pos.x0 + LABEL_X * ax_pos.width
    hct_top_key = 'jund' if 'jund' in idx else 'fosl1'
    if hct_top_key in idx:
        hct_pos = all_axes[idx[hct_top_key]].get_position()
        fig.text(_header_x(hct_pos), hct_pos.y1 + HEADER_OFFSET, 'HCT116',
                 ha='right', va='bottom', fontweight='bold', family=BOLD_FONT,
                 fontsize=HEADER_FONTSIZE, transform=fig.transFigure)
    if 'sig_h3k27ac' in idx:
        sig_pos = all_axes[idx['sig_h3k27ac']].get_position()
        fig.text(_header_x(sig_pos), sig_pos.y1 + HEADER_OFFSET + 0.005,
                 'Sigmoid colon',
                 ha='right', va='bottom', fontweight='bold', family=BOLD_FONT,
                 fontsize=HEADER_FONTSIZE, transform=fig.transFigure)
    top    = max(ax.get_position().y1 for ax in track_axes)
    bottom = min(ax.get_position().y0 for ax in track_axes)
    fig.add_artist(Rectangle(
        (x0, bottom), x1 - x0, top - bottom,
        transform=fig.transFigure, facecolor='orange', alpha=0.2,
        edgecolor=None, zorder=10, clip_on=False,
    ))
    # Compact mode: place the LTR10.XXX label ABOVE the gene track so
    # it doesn't collide with the bold HCT116 header at narrow fig widths
    # (especially when the variant is near the left edge of the window).
    if args.compact and len(sorted_axes) > 0:
        label_y_lt = sorted_axes[0].get_position().y1 + 0.003
    else:
        label_y_lt = top + 0.003
    fig.text((x0 + x1) / 2, label_y_lt, args.variant_id,
             ha='center', va='bottom', fontweight='bold', family=BOLD_FONT,
             fontsize=LTR10_LABEL_SIZE, transform=fig.transFigure)

# Rasterize fills for tractable PDF sizes (and Illustrator compatibility).
for ax in fig.get_axes():
    for c in ax.collections:
        c.set_rasterized(True)

fig.savefig(args.output, dpi=200, bbox_inches='tight')
plt.close()

# Inline recolor: each Tracks() produces one image XObject. With the
# transcript track at axis-index 0 (no image), track-index N maps to
# image /I{N}. Build the Im* → RGB map from our dynamic `idx` dict.
NAVY      = (0x16, 0x33, 0x59)  # TF FOSL1
STEELBLUE = (0x5B, 0x7F, 0xB5)  # TF JUND
DARKRED   = (0x8B, 0x1A, 0x1F)  # H3K27ac
ORANGE    = (0xF0, 0x8C, 0x2D)  # H3K4me1
DARKGREEN = (0x1A, 0x5F, 0x2A)  # ATAC
PURPLE    = (0x7A, 0x33, 0x70)  # RNA (both strands)
TRACK_KIND_COLORS = {
    'jund': STEELBLUE, 'fosl1': NAVY,
    'hct_h3k27ac': DARKRED, 'hct_h3k4me1': ORANGE,
    'hct_atac': DARKGREEN,
    'rna_plus': PURPLE, 'rna_minus': PURPLE,
    'sig_h3k27ac': DARKRED, 'sig_h3k4me1': ORANGE,
}
import pikepdf
def _recolour(img_obj, new_rgb):
    cs = img_obj['/ColorSpace']
    palette = cs[3]
    if isinstance(palette, pikepdf.String):
        raw = bytes(palette)
        cs[3] = pikepdf.String(raw[:3] + bytes(new_rgb))
    else:
        raw = bytes(palette.read_bytes())
        fresh = pdf.make_stream(raw[:3] + bytes(new_rgb))
        img_obj['/ColorSpace'] = pikepdf.Array([cs[0], cs[1], cs[2], fresh])

pdf = pikepdf.open(args.output, allow_overwriting_input=True)
page = pdf.pages[0]
xobjs = page.Resources.XObject
for name, i in idx.items():
    if name == 'transcript' or name not in TRACK_KIND_COLORS:
        continue
    key = f'/I{i}'
    if key in xobjs:
        _recolour(xobjs[key], TRACK_KIND_COLORS[name])
# Strip Illustrator's private artwork cache (PieceInfo et al).
for k in ('/PieceInfo', '/LastModified', '/Thumb'):
    if k in page:
        del page[k]
pdf.save(args.output)
print(f'Wrote {args.output}')
