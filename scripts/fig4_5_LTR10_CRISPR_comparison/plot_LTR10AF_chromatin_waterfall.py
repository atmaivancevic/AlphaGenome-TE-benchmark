"""
Waterfall of the 650 merged LTR10A/F elements ranked by AG-predicted chromatin
peak score (H3K27ac or FOSL1) — Fig 5 Panel A. Bars sorted descending, coloured
grey (no experimental data), orange (113 with experimental peaks), red (the 8
CRISPR-validated fragments, labelled by gene). 2 of the 6 CRISPR enhancers span
two LTR10A/F fragments, so 6 elements appear as 8 bars. Input: Supp Table 11.

Example usage:
python scripts/fig4_5_LTR10_CRISPR_comparison/plot_LTR10AF_chromatin_waterfall.py \
    --input supptables/supp_table_11_LTR10AF_experimental_vs_AG.tsv \
    --mark H3K27ac \
    --output figures/FIG5_FINAL/panelA_LTR10AF_H3K27ac_waterfall.pdf
"""
import argparse
from pathlib import Path
import matplotlib
matplotlib.use('Agg')
matplotlib.rcParams['font.family'] = 'DejaVu Sans'
matplotlib.rcParams['pdf.fonttype'] = 42
import matplotlib.pyplot as plt
from matplotlib.patches import Patch
import pandas as pd

ap = argparse.ArgumentParser()
ap.add_argument('--input',
                default='supptables/supp_table_11_LTR10AF_experimental_vs_AG.tsv',
                help='Supp Table 11 (650 LTR10A/F elements + experimental + AG predictions)')
ap.add_argument('--mark', default='H3K27ac', choices=['H3K27ac', 'FOSL1'],
                help='Which AG-predicted chromatin mark to rank by (default H3K27ac)')
ap.add_argument('--output', required=True)
# Optional overrides — figure dimensions, which CRISPR labels to annotate,
# and the legend text for the orange candidate tier.
ap.add_argument('--width', type=float, default=8.0, help='Figure width in inches')
ap.add_argument('--height', type=float, default=4.2, help='Figure height in inches')
ap.add_argument('--only-label', nargs='*', default=None,
                help='Only annotate these CRISPR labels (e.g. LTR10.ATG12_1). '
                     'Default: annotate all CRISPR-validated bars.')
ap.add_argument('--label-rename', nargs='*', default=None,
                help='Pairs of old=new label renames (e.g. '
                     'LTR10.ATG12_1="LTR10.ATG12 enhancer")')
ap.add_argument('--candidate-name', default='With experimental peak',
                help='Legend text for the orange (candidate-enhancer) tier')
args = ap.parse_args()

label_rename = {}
if args.label_rename:
    for pair in args.label_rename:
        if '=' not in pair:
            raise SystemExit(f'--label-rename expects K=V pairs, got {pair!r}')
        k, v = pair.split('=', 1)
        label_rename[k] = v

# Font sizing — Fig 5 Panel A defaults.
FS_TITLE   = 11
FS_AXIS    = 10
FS_TICK    = 8
FS_LEGEND  = 8
FS_CRISPR  = 6

df = pd.read_table(args.input, sep='\t')
pred_col = f'AG_predicted_{args.mark}_peak_score'
expt_col = f'expt_{args.mark}_peak_score'

# Categorise
df['is_crispr']       = df['CRISPR_validated'].astype(str).str.startswith('yes', na=False)
df['is_experimental'] = df[expt_col].notna()
df['crispr_label']    = df['CRISPR_validated'].astype(str).str.replace(
    'yes:', '', regex=False).where(df['is_crispr'], '')

n_total  = len(df)
n_exp    = int(df['is_experimental'].sum())
n_crispr = int(df['is_crispr'].sum())
n_grey   = n_total - n_exp

# Sort descending by AG prediction
sorted_df = df.sort_values(pred_col, ascending=False).reset_index(drop=True)

fig, ax = plt.subplots(figsize=(args.width, args.height))

# Layer 1 — all 650 elements in grey
ax.bar(range(n_total), sorted_df[pred_col].values,
       width=1.0, color='#CCCCCC', linewidth=0, zorder=2)
# Layer 2 — 113 candidate (experimental peak data) in orange
idx_exp = sorted_df[sorted_df['is_experimental']].index
ax.bar(idx_exp, sorted_df.loc[idx_exp, pred_col].values,
       width=1.0, color='#FF8C00', alpha=1.0, linewidth=0.3, edgecolor='#CC6600', zorder=3)
# Layer 3 — 8 CRISPR-validated fragments in red
idx_crispr = sorted_df[sorted_df['is_crispr']].index
ax.bar(idx_crispr, sorted_df.loc[idx_crispr, pred_col].values,
       width=1.0, color='#E62020', alpha=1.0, linewidth=0.3, edgecolor='darkred', zorder=4)

# Annotate CRISPR-validated bars with their target gene. Default labels all 8;
# --only-label restricts to a subset and draws a leader arrow (for poster scale).
single_label_mode = args.only_label is not None and len(args.only_label) == 1
for idx_val in idx_crispr:
    row = sorted_df.iloc[idx_val]
    crispr_label = row['crispr_label']
    if args.only_label is not None and crispr_label not in args.only_label:
        continue
    display = label_rename.get(crispr_label, crispr_label)
    if single_label_mode:
        # Leader-arrow annotation: text well above the bar, arrow points down.
        ax.annotate(display, xy=(idx_val, row[pred_col]),
                    xytext=(40, 70), textcoords='offset points',
                    fontsize=FS_CRISPR, fontweight='bold', color='#C0180C',
                    ha='left', va='bottom',
                    arrowprops=dict(arrowstyle='->', color='#C0180C',
                                    lw=1.2, shrinkA=0, shrinkB=2),
                    zorder=5)
    else:
        ax.annotate(display, (idx_val, row[pred_col]),
                    xytext=(0, 8), textcoords='offset points',
                    fontsize=FS_CRISPR, fontweight='bold', color='#C0180C',
                    ha='center', va='bottom', rotation=45, zorder=5)

# Legend
handles = [
    Patch(facecolor='#E62020', edgecolor='darkred',
          label=f'CRISPR-validated (n={n_crispr})'),
    Patch(facecolor='#FF8C00', edgecolor='#CC6600',
          label=f'{args.candidate_name} (n={n_exp - n_crispr})'),
    Patch(facecolor='#CCCCCC',
          label=f'No experimental peak (n={n_grey})'),
]
ax.legend(handles=handles, fontsize=FS_LEGEND, loc='upper right')

ax.set_xlabel(f'LTR10A/F elements (n={n_total}, sorted by predicted {args.mark})',
              fontsize=FS_AXIS)
ax.set_ylabel(f'AlphaGenome-predicted {args.mark} peak score', fontsize=FS_AXIS)
ax.tick_params(axis='both', labelsize=FS_TICK)
ax.set_xlim(-1, n_total)
ax.spines['top'].set_visible(False)
ax.spines['right'].set_visible(False)

plt.tight_layout()
out = Path(args.output); out.parent.mkdir(parents=True, exist_ok=True)
plt.savefig(out, dpi=300, bbox_inches='tight')
print(f"Wrote {out}")
print(f"  n_total={n_total}, n_experimental={n_exp}, n_CRISPR={n_crispr}")
