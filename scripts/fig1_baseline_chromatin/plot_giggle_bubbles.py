"""
Bubble plot of GIGGLE TE enrichment in ENCODE H3K27ac peaks.

Uses the 52 TE subfamilies that passed filters (OR>5, score>=100, overlaps>40)
in at least one biosample, but plots their enrichment across ALL biosamples
from the unfiltered GIGGLE results.

Usage:
    python scripts/fig1_baseline_chromatin/plot_giggle_bubbles.py [--top-TEs N] [--top-biosamples N]
"""

import os, argparse
import pandas as pd
import plotly.express as px

SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
PROJECT_ROOT = os.path.dirname(os.path.dirname(SCRIPT_DIR))

parser = argparse.ArgumentParser()
parser.add_argument('--top-TEs', type=int, default=None,
                    help='Show only top N TE subfamilies by max GIGGLE score')
parser.add_argument('--top-biosamples', type=int, default=None,
                    help='Show only top N biosamples by max GIGGLE score')
parser.add_argument('--keep-TEs', default=None,
                    help='Comma-separated TE subfamilies to keep '
                         '(exact repeat names; e.g. LTR10A,LTR10F,LTR2B). '
                         'Plotted in the order given. Overrides --top-TEs.')
parser.add_argument('--keep-biosamples', default=None,
                    help='Comma-separated biosamples to keep — match against the '
                         'human-readable biosample_name (e.g. HCT116,GM12878,K562). '
                         'Plotted in the order given. Overrides --top-biosamples.')
parser.add_argument('--flip', action='store_true',
                    help='Put biosamples on x-axis, TEs on y-axis')
parser.add_argument('--suffix', default=None,
                    help='Override the auto-suffix used in the output filename')
parser.add_argument('--width', type=int, default=None,
                    help='Override figure width in pixels (≈ pt at 72 DPI). '
                         'Default: auto-sized based on number of TEs.')
parser.add_argument('--height', type=int, default=None,
                    help='Override figure height in pixels (≈ pt at 72 DPI). '
                         'Default: auto-sized based on number of biosamples.')
parser.add_argument('--size-max', type=int, default=30,
                    help='Max bubble diameter for the largest GIGGLE score. '
                         'Lower it to reduce overlap of the big bubbles.')
parser.add_argument('--outdir', default='figures',
                    help='Output directory, relative to the project root '
                         '(default: figures/).')
parser.add_argument('--render-scale', type=float, default=1.0,
                    help='Uniformly scale the rendered canvas, fonts, margins '
                         'and bubbles by this factor. --width/--height stay the '
                         'DESIGN size; the PDF is rendered larger and meant to '
                         'be placed at 1/scale in the layout tool. Use 2.0 to '
                         'dodge the kaleido small-figure export bug (Error 525).')
args = parser.parse_args()
SCALE = args.render_scale

# 52 TE subfamilies that passed filters (OR>5, score>=100, overlaps>40)
filtered = pd.read_csv(os.path.join(PROJECT_ROOT,
    'data/encode_h3k27ac/giggle_results/filtered_OR5_score100_overlaps40.tab'), sep='\t')
keep_tes = filtered['repeat'].unique().tolist()

# Load unfiltered results, keep only those 52 subfamilies (unless --keep-TEs
# specifies an explicit list, in which case use that instead — supports
# picking TE families that didn't pass the global OR/score/overlap filter
# but are biologically interesting for a curated panel).
df = pd.read_csv(os.path.join(PROJECT_ROOT,
    'data/encode_h3k27ac/giggle_results/h3k27ac_vs_TEs_rankedByScore.tab'), sep='\t')

# Drop rows with non-positive scores (no enrichment)
df = df[df['gigglescore'] > 0]

# --- TE order ---
if args.keep_TEs:
    keep = [x.strip() for x in args.keep_TEs.split(',') if x.strip()]
    df = df[df['repeat'].isin(keep)]
    # Sort by max GIGGLE score within the curated subset (strongest first),
    # matching the supp-plot convention.
    te_order = df.groupby('repeat')['gigglescore'].max().sort_values(ascending=False).index.tolist()
else:
    df = df[df['repeat'].isin(keep_tes)]
    te_order = filtered.groupby('repeat')['gigglescore'].max().sort_values(ascending=False).index.tolist()
    if args.top_TEs:
        # Drop MER57E3 — constitutively active in nearly all cell types, not cell-type-specific
        te_order = [t for t in te_order if t != 'MER57E3']
        te_order = te_order[:args.top_TEs]
        df = df[df['repeat'].isin(te_order)]

# Clean biosample names (needed before --keep-biosamples filter)
df['biosample_name'] = df['biosample'].str.replace(r'^[A-Z]+_\d+_', '', regex=True).str.replace('_', ' ')
name_map = df.drop_duplicates('biosample').set_index('biosample')['biosample_name'].to_dict()
inverse_map = {v: k for k, v in name_map.items()}

# --- Biosample order ---
if args.keep_biosamples:
    requested = [x.strip() for x in args.keep_biosamples.split(',') if x.strip()]
    missing = [r for r in requested if r not in inverse_map]
    if missing:
        raise ValueError(f"--keep-biosamples names not found in data: {missing}\n"
                         f"Available examples: {list(inverse_map.keys())[:10]}")
    df = df[df['biosample'].isin([inverse_map[r] for r in requested])]
    # Sort by max GIGGLE score within the curated subset (strongest first),
    # matching the supp-plot convention.
    bs_order = df.groupby('biosample')['gigglescore'].max().sort_values(ascending=False).index.tolist()
    bs_name_order = [name_map[b] for b in bs_order]
else:
    bs_order = df.groupby('biosample')['gigglescore'].max().sort_values(ascending=False).index.tolist()
    if args.top_biosamples:
        bs_order = bs_order[:args.top_biosamples]
        df = df[df['biosample'].isin(bs_order)]
    bs_name_order = [name_map[b] for b in bs_order]

n_tes = len(te_order)
n_bs = len(bs_order)
if args.suffix:
    suffix = args.suffix if args.suffix.startswith('_') else '_' + args.suffix
else:
    parts = []
    if args.top_TEs:        parts.append(f'top{args.top_TEs}')
    if args.top_biosamples: parts.append(f'bs{args.top_biosamples}')
    if args.keep_TEs:       parts.append(f'TEs{len(args.keep_TEs.split(","))}')
    if args.keep_biosamples: parts.append(f'bs{len(args.keep_biosamples.split(","))}')
    suffix = '_' + '_'.join(parts) if parts else '_all'

if args.flip:
    # Biosamples on x, TEs on y
    w_design = args.width  if args.width  else max(800, 200 + n_bs * 28)
    h_design = args.height if args.height else max(600, 80 + n_tes * 28)
    x_col, y_col = 'biosample_name', 'repeat'
    x_label = f'<b>Biosample (n = {n_bs})</b>'
    y_label = f'<b>TE subfamily (n = {n_tes})</b>'
    x_range, y_range = [-1, n_bs], [-1, n_tes]
    cat_orders = {'biosample_name': bs_name_order, 'repeat': te_order}
else:
    # TEs on x, biosamples on y (default)
    w_design = args.width  if args.width  else max(800, 200 + n_tes * 28)
    h_design = args.height if args.height else max(600, 80 + n_bs * 28)
    x_col, y_col = 'repeat', 'biosample_name'
    x_label = f'<b>TE subfamily (n = {n_tes})</b>'
    y_label = f'<b>Biosample (n = {n_bs})</b>'
    x_range, y_range = [-1, n_tes], [-1, n_bs]
    cat_orders = {'biosample_name': bs_name_order, 'repeat': te_order}

# Rendered size = design size x render-scale (place at 1/scale in layout tool).
w, h = w_design * SCALE, h_design * SCALE

SIZE_MAX = args.size_max * SCALE   # max bubble diameter; --size-max to tune overlap

fig = px.scatter(df, x=x_col, y=y_col,
                 size='gigglescore', color='oddsratio',
                 size_max=SIZE_MAX, opacity=1,
                 color_continuous_scale=[(0.0, '#FD9367'), (0.33, '#C3305D'),
                                         (0.67, '#782D65'), (1, '#432967')],
                 labels={'biosample_name': x_label if args.flip else y_label,
                         'gigglescore': 'Enrichment score',
                         'oddsratio': '<b>Odds ratio</b>',
                         'repeat': y_label if args.flip else x_label},
                 category_orders=cat_orders)

# Auto-shrink margins/font when the user explicitly downsizes the figure
# (e.g. for a main-figure mini panel at 260 × 200 pt). Keyed on the DESIGN
# size, then scaled by SCALE so the placed-at-1/scale result is unchanged.
_small = (args.width and args.width < 400) or (args.height and args.height < 300)
_top_margin = (50 if _small else 120) * SCALE
_font_size = (7 if _small else 10) * SCALE
fig.update_layout(template='plotly_white', autosize=False, width=w, height=h,
                  font=dict(size=_font_size, family='Helvetica'),
                  margin=dict(l=10 * SCALE, r=10 * SCALE, t=_top_margin, b=10 * SCALE),
                  coloraxis_colorbar=dict(len=0.75, yanchor='top', y=1, x=1.02))

# No in-figure size legend — Plotly Express renders shape circles at a
# slightly different effective scale from scatter markers (different
# rendering paths), so an in-figure reference would mislead. The size
# legend is built manually in Illustrator using existing data bubbles
# as size references (e.g. H1 + MER41B ~50, placenta + LTR5_Hs ~100,
# placenta + LTR10A ~500, HUES6 + LTR5_Hs ~1000).
fig.update_yaxes(tickfont_size=9 * SCALE, ticks='outside', ticklen=4 * SCALE,
                 showline=True, linecolor='black', linewidth=1 * SCALE,
                 mirror=True, range=y_range,
                 title=dict(font=dict(size=16 * SCALE)))
fig.update_xaxes(tickfont_size=10 * SCALE, tickangle=-45, ticks='outside', ticklen=4 * SCALE,
                 showline=True, linecolor='black', linewidth=1 * SCALE,
                 mirror=True, range=x_range,
                 title=dict(font=dict(size=16 * SCALE)),
                 side='top')

outdir = os.path.join(PROJECT_ROOT, args.outdir)
os.makedirs(outdir, exist_ok=True)
out_pdf = os.path.join(outdir, f'fig1_giggle_bubbles{suffix}.pdf')
out_png = os.path.join(outdir, f'fig1_giggle_bubbles{suffix}.png')
fig.write_image(out_pdf)
fig.write_image(out_png, scale=4)
print(f"{n_tes} TEs x {n_bs} biosamples | {len(df)} associations plotted")
print(f"Saved: {out_pdf}")
print(f"Saved: {out_png} (high-res)")
