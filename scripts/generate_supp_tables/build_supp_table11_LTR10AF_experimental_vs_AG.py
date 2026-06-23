"""
Build Supp Table 11 — experimental vs AlphaGenome-predicted chromatin peak
scores across all 650 merged LTR10A/F elements. One row per element.

This is the variant-level companion to Supp Table 10's gene-level CRISPR/AG
comparison: where Supp 10 asks "does AG predict the gene-expression effect
of CRISPR-perturbing each of 6 LTR10 enhancers?", Supp 11 asks "does AG
correctly rank the H3K27ac/FOSL1 activity of LTR10 enhancers genome-wide?".

Source data (already built by prior work):
  data/650_LTR10AF_with_experimental_and_predicted.tsv
    - 650 merged LTR10A/F elements (chr1-22 + chrX), peak scores from
      Ivancevic 2024 in-house HCT116 H3K27ac + FOSL1 ChIP-seq (113/650
      elements have experimental peak signal; rest are NA)
    - AG predictions generated with scripts/wip/predict_baseline_ChIP.py
    - 8 elements marked CRISPR_validated (2 fragments each for LTR10.ATG12
      and LTR10.XRCC4 — the merged set splits these long elements; the
      remaining 4 CRISPR-tested elements are single fragments)

This script adds two derived columns (length_bp, family) and writes the
result to supptables/. Notes / crispr_target_gene columns are intentionally
omitted — the CRISPR_validated column already encodes both.

TODO (deferred, tracked in todos):
  - Replace in-house HCT116 H3K27ac peak scores with ENCODE HCT116 H3K27ac
    (matches the data AG was trained on; current in-house calls are
    callable but use slightly different peak-calling settings).
  - Score AG ATAC + H3K4me1 across all 650 elements and append columns
    (currently Supp 11 has only FOSL1 + H3K27ac since those are what
    scripts/wip/predict_baseline_ChIP.py was originally run for).

Usage:
    python scripts/generate_supp_tables/build_supp_table11_LTR10AF_experimental_vs_AG.py \\
        --input data/650_LTR10AF_with_experimental_and_predicted.tsv \\
        --tsv   supptables/supp_table_11_LTR10AF_experimental_vs_AG.tsv
"""
import argparse
from pathlib import Path
import pandas as pd

p = argparse.ArgumentParser()
p.add_argument('--input', default='data/650_LTR10AF_with_experimental_and_predicted.tsv',
               help='Source 650-element table (prior work)')
p.add_argument('--tsv',   required=True, help='Output supp table path')
args = p.parse_args()

df = pd.read_table(args.input, sep='\t')

# length_bp = end - start_1based + 1 (1-based inclusive coords)
df['length_bp'] = (df['end'] - df['start_1based'] + 1).astype(int)
# family = LTR10A or LTR10F from the element_id prefix
df['family'] = df['element_id'].str.split('_').str[0]
# candidate = yes for the 113 elements with experimental peak scores
# (these are the highlighted/yellow points on the Fig 5 / Fig 1
# correlation + waterfall plots — the "candidate" enhancer set picked for
# follow-up in our paper). When the experimental dataset gets swapped to
# ENCODE HCT116 (deferred), this column preserves the original published
# candidate set so we don't lose its identity.
df['candidate'] = df['expt_H3K27ac_peak_score'].notna().map(
    {True: 'yes', False: 'no'})

# Final column order: identity + family/length first, then the two
# label columns (candidate + CRISPR_validated), then experimental scores,
# then AG-predicted (mirrors the Supp 10 layout where ground-truth
# columns precede model predictions).
cols = ['element_id', 'chrom', 'start_1based', 'end', 'length_bp', 'family',
        'candidate', 'CRISPR_validated',
        'expt_FOSL1_peak_score', 'expt_H3K27ac_peak_score',
        'AG_predicted_FOSL1_peak_score', 'AG_predicted_H3K27ac_peak_score']
out = df[cols]
Path(args.tsv).parent.mkdir(parents=True, exist_ok=True)
out.to_csv(args.tsv, sep='\t', index=False)
print(f"Wrote {args.tsv}  ({len(out)} rows x {out.shape[1]} cols)")

# Diagnostics
print(f"\nFamily breakdown:")
print(out['family'].value_counts().to_string())
print(f"\nCRISPR_validated rows ({(out['CRISPR_validated'] != 'no').sum()}):")
print(out[out['CRISPR_validated'] != 'no'][
    ['element_id','chrom','length_bp','CRISPR_validated',
     'expt_H3K27ac_peak_score','AG_predicted_H3K27ac_peak_score']
].to_string(index=False))
print(f"\nExperimental data coverage: {out['expt_H3K27ac_peak_score'].notna().sum()}/{len(out)} elements")
