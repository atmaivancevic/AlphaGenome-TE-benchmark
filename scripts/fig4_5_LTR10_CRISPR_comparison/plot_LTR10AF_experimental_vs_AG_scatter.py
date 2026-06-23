"""
Scatter of experimental ChIP-seq vs AG-predicted peak scores across the 113
LTR10A/F elements with experimental data (H3K27ac or FOSL1) — Fig 5 Panel B.
Both axes are non-negative peak intensities (correlation reported signed only).
Orange = the 113 candidates, red = the 8 CRISPR-validated fragments (labelled).
Pearson r + Spearman rho annotated. Input: Supp Table 11.

Example usage:
python scripts/fig4_5_LTR10_CRISPR_comparison/plot_LTR10AF_experimental_vs_AG_scatter.py \
    --input supptables/supp_table_11_LTR10AF_experimental_vs_AG.tsv \
    --mark H3K27ac \
    --output figures/FIG5_FINAL/panelB_LTR10AF_H3K27ac_scatter.pdf
"""
import argparse
from pathlib import Path
import matplotlib
matplotlib.use('Agg')
matplotlib.rcParams['font.family'] = 'DejaVu Sans'
matplotlib.rcParams['pdf.fonttype'] = 42
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy import stats

ap = argparse.ArgumentParser()
ap.add_argument('--input',
                default='supptables/supp_table_11_LTR10AF_experimental_vs_AG.tsv',
                help='Supp Table 11 (650 LTR10A/F elements + experimental + AG predictions)')
ap.add_argument('--mark', default='H3K27ac', choices=['H3K27ac', 'FOSL1'],
                help='Which chromatin mark to plot (default H3K27ac)')
ap.add_argument('--output', required=True)
args = ap.parse_args()

df = pd.read_table(args.input, sep='\t')
expt_col = f'expt_{args.mark}_peak_score'
pred_col = f'AG_predicted_{args.mark}_peak_score'

# Keep only elements with experimental data
matched = df[df[expt_col].notna()].copy()
matched['is_crispr']    = matched['CRISPR_validated'].astype(str).str.startswith('yes', na=False)
matched['crispr_label'] = matched['CRISPR_validated'].astype(str).str.replace(
    'yes:', '', regex=False).where(matched['is_crispr'], '')

non_crispr = matched[~matched['is_crispr']]
crispr     = matched[matched['is_crispr']]
n_total  = len(matched)
n_crispr = int(matched['is_crispr'].sum())

# Stats signed only (x, y are non-negative peak intensities, so unsigned is identical).
x = matched[expt_col].values.astype(float)
y = matched[pred_col].values.astype(float)
r,  p_r   = stats.pearsonr(x, y)
rho, p_rho = stats.spearmanr(x, y)
print(f"{args.mark} (n={n_total}):  Pearson r = {r:+.3f} (p={p_r:.2e})   "
      f"Spearman ρ = {rho:+.3f} (p={p_rho:.2e})")

# Plot
fig, ax = plt.subplots(figsize=(5.2, 4.2))

ax.scatter(non_crispr[expt_col], non_crispr[pred_col],
           s=28, color='#FF8C00', alpha=0.75, edgecolors='#CC6600', linewidths=0.3,
           zorder=3, label=f'With experimental peak (n={len(non_crispr)})')
ax.scatter(crispr[expt_col], crispr[pred_col],
           s=55, color='#E62020', alpha=1.0, edgecolors='darkred', linewidths=0.5,
           zorder=4, label=f'CRISPR-validated (n={n_crispr})')

for _, row in crispr.iterrows():
    ax.annotate(row['crispr_label'], (row[expt_col], row[pred_col]),
                xytext=(5, 5), textcoords='offset points',
                fontsize=7, fontweight='bold', color='#C0180C', zorder=5)

# Best-fit dashed line
slope, intercept = np.polyfit(x, y, 1)
x_line = np.linspace(x.min(), x.max(), 100)
ax.plot(x_line, slope * x_line + intercept,
        color='#333333', linewidth=1.2, linestyle='--', zorder=2)

# Stats top-left (same placement convention as Fig 3 G/H and Fig 5 B-G)
ax.text(0.05, 0.95,
        f'Pearson r = {r:.2f}, p = {p_r:.1e}\n'
        f'Spearman ρ = {rho:.2f}, p = {p_rho:.1e}',
        transform=ax.transAxes, fontsize=9, va='top',
        bbox=dict(boxstyle='round,pad=0.3', fc='white', ec='#CCCCCC', alpha=0.9))

ax.set_xlabel(f'Experimental {args.mark} peak score', fontsize=10)
ax.set_ylabel(f'AlphaGenome-predicted {args.mark} peak score', fontsize=10)
ax.legend(fontsize=8, loc='lower right')
ax.spines['top'].set_visible(False)
ax.spines['right'].set_visible(False)

plt.tight_layout()
out = Path(args.output); out.parent.mkdir(parents=True, exist_ok=True)
plt.savefig(out, dpi=300, bbox_inches='tight')
print(f"Wrote {out}")
