"""
Build Supp Table 11 — experimental vs AlphaGenome-predicted chromatin peak scores
(H3K27ac, FOSL1) across all 650 merged LTR10A/F elements (one row per element).
Variant-level companion to Supp Table 10: "does AG rank LTR10 enhancer activity
correctly genome-wide?". 113/650 elements have experimental peaks; 8 fragments are
CRISPR-validated. Adds derived length_bp + family columns. Input:
data/650_LTR10AF_with_experimental_and_predicted.tsv.
TODO: swap in-house HCT116 H3K27ac for ENCODE HCT116; add AG ATAC + H3K4me1 columns.

Example usage:
python scripts/generate_supp_tables/build_supp_table11_LTR10AF_experimental_vs_AG.py \
    --input data/650_LTR10AF_with_experimental_and_predicted.tsv \
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
