"""
Fig 6 Panel B — TF substitution specificity at LTR10.ATG12. Grouped bar chart of
AG-predicted H3K27ac + FOSL1 across WT, scram20 (all AP1 scrambled), and the four
TF substitutions (TP53/CTCF/GATA1/HNF1A). Tests whether the predicted output is
AP1-specific: if so, all substitutions collapse to scram20 levels.

Example usage:
python scripts/fig6_AP1_perturbation/plot_fig6_panelB_TF_specificity.py \
    --chromatin-csv results/AG_perturbation_LTR10_ATG12_AP1/predict_chromatin_at_element.csv \
    --output figures/FIG6_FINAL/panelB_TF_specificity.pdf
"""
import argparse
from pathlib import Path
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
matplotlib.rcParams['font.family'] = 'DejaVu Sans'
matplotlib.rcParams['pdf.fonttype'] = 42
import matplotlib.pyplot as plt

ap = argparse.ArgumentParser()
ap.add_argument('--chromatin-csv', required=True)
ap.add_argument('--output', required=True)
args = ap.parse_args()

df = pd.read_csv(args.chromatin_csv)

# Pivot: variant_id × (H3K27ac, FOSL1) — ALT signal
def is_track(s, key): return s.str.contains(key, case=False, na=False)
df['mark'] = np.where(is_track(df['track_name'], 'H3K27ac'), 'H3K27ac',
              np.where(is_track(df['track_name'], 'FOSL1'),   'FOSL1', '?'))
ag = df.pivot_table(index='variant_id', columns='mark', values='signal_mean_ALT', aggfunc='mean')

# Include WT as a row (REF prediction, identical across alleles). Dict assignment
# keeps column-name alignment explicit (the pivot orders columns alphabetically).
wt_h3k = df[df['mark']=='H3K27ac']['signal_mean_REF'].mean()
wt_fosl1 = df[df['mark']=='FOSL1']['signal_mean_REF'].mean()
ag.loc['WT', 'H3K27ac'] = wt_h3k
ag.loc['WT', 'FOSL1']   = wt_fosl1

# Order alleles for plot: WT, scram20 (full AP1 disruption), then 4 TF substitutions
ORDER = ['WT',
         'LTR10.ATG12_scram20',
         'LTR10.ATG12_allTF_TP53',
         'LTR10.ATG12_allTF_CTCF',
         'LTR10.ATG12_allTF_GATA1',
         'LTR10.ATG12_allTF_HNF1A']
LABELS = ['WT', 'scrambled',
          'AP1→TP53', 'AP1→CTCF', 'AP1→GATA1', 'AP1→HNF1A']

ag = ag.reindex(ORDER)
print(ag.to_string(float_format='%.2f'))

# Grouped bar chart
fig, ax = plt.subplots(figsize=(6.5, 4.2))
x = np.arange(len(ORDER))
w = 0.38
ax.bar(x - w/2, ag['H3K27ac'].values, width=w,
       color='#5BA75A', edgecolor='white', linewidth=0.5,
       label='Predicted H3K27ac')
ax.bar(x + w/2, ag['FOSL1'].values, width=w,
       color='#D7301F', edgecolor='white', linewidth=0.5,
       label='Predicted FOSL1')

# WT reference horizontal lines
ax.axhline(wt_h3k, color='#5BA75A', linestyle=':', linewidth=0.7, alpha=0.5)
ax.axhline(wt_fosl1, color='#D7301F', linestyle=':', linewidth=0.7, alpha=0.5)

ax.set_xticks(x)
ax.set_xticklabels(LABELS, rotation=30, ha='right')
ax.set_ylabel('Predicted signal (signal_mean over LTR10.ATG12)', fontsize=10)
ax.set_xlabel('LTR10.ATG12 allele (20 AP1 motifs each, replaced with stated TF consensus)', fontsize=10)
ax.spines['top'].set_visible(False)
ax.spines['right'].set_visible(False)
ax.legend(fontsize=8, loc='upper right', frameon=False)

plt.tight_layout()
out = Path(args.output); out.parent.mkdir(parents=True, exist_ok=True)
plt.savefig(out, dpi=300, bbox_inches='tight')
print(f"Wrote {out}")
