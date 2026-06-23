"""
Coordinate scatter of AlphaGenome RNA-seq variant effects for a single LTR10
CRISPRi deletion. X = gene TSS, Y = AG raw_score (log2FC), point colour by
|quantile_score| threshold. Reproduces the layout of
plot_variant_effects_coord_scatterplot() from scripts/wip/old/helper_functions.py
with one substantive change: instead of `track_filter='max_effect'` (pick the
single highest-|q| track per gene), per-gene scores here are the **mean of the
strand-specific RNA tracks only** — i.e. unstranded polyA tracks are dropped
before averaging (AG track aggregation convention, 2026-05-11).

For HCT116 (EFO:0002824), AG emits two tracks per gene:
  - `total RNA-seq`  (track_strand = + or -, matched to gene strand)  ← KEPT
  - `polyA plus RNA-seq`  (track_strand = '.', unstranded)            ← DROPPED
So the per-gene score collapses to the single stranded total RNA-seq track for
HCT116. The same rule generalises to GM12878 (which has two stranded tracks,
so the average is over both).

Inputs:
  --variants     data/LTR10_variants.tab
  --rna-csv      results/AG_LFC_LTR10_CRISPR/<variant>_<cell>.csv
  --gencode      data/gencode.v46.annotation.feather (TSS lookup)

Usage (LTR10.ATG12 in HCT116, matching the legacy call):
    python scripts/fig4_5_LTR10_CRISPR_comparison/plot_variant_coord_scatter.py \\
        --variant LTR10.ATG12 \\
        --cell-line HCT116 \\
        --variants data/LTR10_variants.tab \\
        --rna-csv results/AG_LFC_LTR10_CRISPR/LTR10.ATG12_HCT116.csv \\
        --gencode data/gencode.v46.annotation.feather \\
        --gene-types protein_coding \\
        --color-threshold 0.9 \\
        --ylim -0.3 0.3 \\
        --output figures/FIG4_LTR10/LTR10.ATG12_HCT116_coord_scatter.pdf
"""
import argparse
from pathlib import Path
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker
from matplotlib.patches import Rectangle
from adjustText import adjust_text

ap = argparse.ArgumentParser()
ap.add_argument('--variant', required=True, help='Variant ID (e.g. LTR10.ATG12)')
ap.add_argument('--cell-line', default='HCT116')
ap.add_argument('--variants', required=True, help='LTR10_variants.tab (ID/CHROM/POS/REF/ALT/...)')
ap.add_argument('--rna-csv', required=True, help='Per-variant AG RNA CSV from score_variant_lfc.py')
ap.add_argument('--gencode', required=True, help='gencode.v46.annotation.feather')
ap.add_argument('--gene-types', default='protein_coding',
                help="'all', 'protein_coding', or comma-separated list (default: protein_coding)")
ap.add_argument('--color-threshold', type=float, default=0.9,
                help='|quantile_score| above which a gene is coloured (default 0.9, per Avsec 2026)')
ap.add_argument('--window-mb', type=float, default=10.0,
                help='Plot window width centred on the variant midpoint (default 10 Mb, matching the legacy CRISPRi coord-scatter zoom)')
ap.add_argument('--xlim', type=float, nargs=2, default=None, help='Manual xlim (overrides --window-mb)')
ap.add_argument('--ylim', type=float, nargs=2, default=None, help='Manual ylim (default: auto symmetric)')
ap.add_argument('--title', default=None,
                help='Plot title (default: "<variant> Predicted Effects of Deletion")')
ap.add_argument('--xlab', default='Genome Coordinates (bp)')
ap.add_argument('--ylab', default='Raw Variant Score (Effect Magnitude)')
ap.add_argument('--output', required=True, help='Output PDF path')
args = ap.parse_args()

# ---- Variant ----
vdf = pd.read_table(args.variants, sep='\t')
vrow = vdf[vdf['ID'] == args.variant]
if vrow.empty:
    raise SystemExit(f"Variant {args.variant!r} not found in {args.variants}")
vrow = vrow.iloc[0]
chrom = vrow['CHROM']
ref_len = len(vrow['REF'])
var_start = int(vrow['POS'])
var_end   = var_start + ref_len
var_mid   = (var_start + var_end) // 2

# ---- AG RNA scores: keep strand-specific tracks only, average per gene ----
rna = pd.read_csv(args.rna_csv)
n_total = len(rna)
strand_specific = rna[rna['track_strand'].isin(['+', '-'])].copy()
n_kept = len(strand_specific)
n_tracks_per_gene = strand_specific.groupby('gene_id').size()
print(f"Loaded {n_total} (gene, track) rows for {args.variant} in {args.cell_line}")
print(f"  kept {n_kept} rows from strand-specific RNA tracks "
      f"(dropped {n_total - n_kept} unstranded rows)")
print(f"  tracks per gene after filter: min={n_tracks_per_gene.min()}  "
      f"max={n_tracks_per_gene.max()}  median={int(n_tracks_per_gene.median())}")

per_gene = (strand_specific
            .groupby(['gene_id','gene_name','gene_type','gene_strand'], as_index=False)
            .agg(raw_score=('raw_score','mean'),
                 quantile_score=('quantile_score','mean'),
                 n_tracks=('raw_score','size')))

# Gene-type filter
if args.gene_types != 'all':
    types = [t.strip() for t in args.gene_types.split(',')]
    before = len(per_gene)
    per_gene = per_gene[per_gene['gene_type'].isin(types)].copy()
    print(f"  filtered to gene_type ∈ {types}: {len(per_gene)}/{before} genes")

# ---- TSS lookup from GENCODE ----
gtf = pd.read_feather(args.gencode)
gene_rows = gtf[(gtf['Feature'] == 'gene') & (gtf['Chromosome'] == chrom)].copy()
gene_rows['TSS'] = gene_rows.apply(
    lambda r: r['Start'] if r['Strand'] == '+' else r['End'], axis=1)
gene_rows['gene_id_base'] = gene_rows['gene_id'].str.split('.').str[0]

# Build name-keyed + id-keyed lookup
tss_by_name = dict(zip(gene_rows['gene_name'], gene_rows['TSS']))
tss_by_id   = dict(zip(gene_rows['gene_id_base'], gene_rows['TSS']))

def lookup_tss(row):
    if row['gene_name'] in tss_by_name:
        return tss_by_name[row['gene_name']]
    gid = str(row['gene_id']).split('.')[0]
    return tss_by_id.get(gid, pd.NA)
per_gene['TSS'] = per_gene.apply(lookup_tss, axis=1)
n_no_tss = per_gene['TSS'].isna().sum()
if n_no_tss:
    miss = per_gene[per_gene['TSS'].isna()][['gene_name','gene_id']]
    print(f"  WARNING: {n_no_tss} genes have no GENCODE TSS (dropped):")
    print(miss.to_string(index=False))
per_gene = per_gene.dropna(subset=['TSS']).copy()
per_gene['TSS'] = per_gene['TSS'].astype(int)
print(f"  plotting {len(per_gene)} genes")

# ---- Categorise by colour ----
thr = args.color_threshold
def cat(q):
    if abs(q) < thr: return 'minimal'
    return 'decrease' if q < 0 else 'increase'
per_gene['cat'] = per_gene['quantile_score'].apply(cat)
COLORS = {'decrease': '#FF0000', 'minimal': '#808080', 'increase': '#0000FF'}
SIZES  = {'decrease': 40,        'minimal': 20,        'increase': 40}

# ---- Plot ----
fig, ax = plt.subplots(figsize=(5.5, 4))

ax.scatter(per_gene['TSS'], per_gene['raw_score'],
           c=per_gene['cat'].map(COLORS),
           s=per_gene['cat'].map(SIZES),
           alpha=0.75, edgecolors='none', linewidth=0)

# Label coloured genes, with adjustText collision avoidance + leader lines.
# Labels start at the dot and adjust_text iteratively pushes them apart;
# arrowprops draws a thin connector back to the dot (matches the legacy
# 'ARL14EPL ──•' style in Fig 4D).
texts = []
for _, r in per_gene[per_gene['cat'] != 'minimal'].iterrows():
    texts.append(ax.text(r['TSS'], r['raw_score'], r['gene_name'],
                         fontsize=11, alpha=0.9, ha='center'))
if texts:
    # Aggressive collision avoidance: wider bbox padding, stronger inter-label
    # repulsion, push away from scatter points too so labels don't sit on dots.
    # Leader lines are drawn from the moved label back to the dot.
    adjust_text(
        texts, ax=ax,
        arrowprops=dict(arrowstyle='-', color='black', lw=0.5, alpha=0.6),
        expand=(1.6, 2.0),
        force_text=(0.8, 1.2),
        force_static=(0.6, 0.8),
        force_pull=(0.01, 0.01),
        max_move=80,
        only_move={'text': 'xy', 'static': 'xy'},
    )

# Axis limits
if args.xlim is not None:
    ax.set_xlim(args.xlim[0], args.xlim[1])
else:
    half = int(args.window_mb * 1_000_000 / 2)
    ax.set_xlim(var_mid - half, var_mid + half)

if args.ylim is not None:
    ax.set_ylim(args.ylim[0], args.ylim[1])
else:
    y0, y1 = ax.get_ylim()
    m = max(abs(y0), abs(y1)) * 1.05
    ax.set_ylim(-m, m)

# Variant rectangle at y=0
y0, y1 = ax.get_ylim()
rect_h = (y1 - y0) * 0.045
x0, x1 = ax.get_xlim()
min_display = (x1 - x0) * 0.04
rect_w = max(var_end - var_start, min_display)
ax.add_patch(Rectangle((var_start, -rect_h/2), rect_w, rect_h,
                       facecolor='black', edgecolor='black', linewidth=0, zorder=3))
ax.axhline(0, color='black', linewidth=1, alpha=0.5)

# Labels + ticks
ax.set_xlabel(args.xlab, fontsize=10)
ax.set_ylabel(args.ylab, fontsize=10)
title = args.title if args.title else f'{args.variant} Predicted Effects of Deletion'
ax.set_title(title, fontsize=12, fontweight='bold')
ax.xaxis.set_major_formatter(plt.FuncFormatter(lambda x, _: f'{int(x):,}'))
ax.text(-0.13, -0.08, chrom, transform=ax.transAxes,
        fontsize=15, fontweight='bold', va='top', clip_on=False)
for spine in ax.spines.values():
    spine.set_linewidth(1.5)
ax.tick_params(width=1.5, labelsize=8)
ax.yaxis.set_major_locator(ticker.MultipleLocator(0.1))

# Scale bar (bottom-left)
xlims, ylims = ax.get_xlim(), ax.get_ylim()
scale_len = (xlims[1] - xlims[0]) / args.window_mb * 1.5 if args.window_mb else 1_500_000
scale_len = min(scale_len, (xlims[1] - xlims[0]) * 0.4)
bar_x = xlims[0] + (xlims[1] - xlims[0]) * 0.03
bar_y = ylims[0] + (ylims[1] - ylims[0]) * 0.06
ax.plot([bar_x, bar_x + scale_len], [bar_y, bar_y],
        color='black', linewidth=1.5, solid_capstyle='butt')
ax.text(bar_x, bar_y + (ylims[1] - ylims[0]) * 0.06,
        f'{scale_len/1e6:.1f} Mb', fontsize=10, va='top')

plt.tight_layout()
out = Path(args.output); out.parent.mkdir(parents=True, exist_ok=True)
plt.savefig(out, dpi=300, bbox_inches='tight')
print(f"Wrote {out}")
