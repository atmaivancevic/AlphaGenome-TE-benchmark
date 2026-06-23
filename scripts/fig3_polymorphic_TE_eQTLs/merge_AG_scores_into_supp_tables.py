"""
Merge AlphaGenome per-(variant, gene) RNA + per-variant chromatin scores into
Supplementary Tables 4 and 5.

Inputs (two AG scoring masters):
- results/AG_LFC_polymorphic_TE/master_ranked_GM12878.csv
    Per-(variant, gene, RNA-seq-track) raw_score + quantile_score from
    scripts/score_variant_lfc.py. RNA scores are aggregated to one number
    per (variant, gene) by taking the mean across all gene-strand-matched
    RNA-seq tracks AG returned for GM12878 — typically 3: stranded total +
    stranded polyA + unstranded polyA. This averaging follows the
    `--track_filter average` convention from our prior baseline-scoring
    pipeline (scripts/predict_baseline_RNA.py).
- results/AG_chromatin_polymorphic_TE/master_chromatin_GM12878.csv
    Per-(variant, mark) raw_score + quantile_score from
    scripts/score_variant_chromatin.py. GM12878 has exactly one track per
    mark for our 3 chosen marks (H3K27ac, H3K4me1, ATAC), so no
    aggregation is needed.

Outputs (8 new AG columns per supp table, in this order):
- AG_RNA_raw_score, AG_RNA_quantile_score
- AG_H3K27ac_raw_score, AG_H3K27ac_quantile_score
- AG_H3K4me1_raw_score, AG_H3K4me1_quantile_score
- AG_ATAC_raw_score, AG_ATAC_quantile_score

Placement: appended after the existing eQTL columns so AG predictions sit
next to the observed β / p / q values for direct comparison.

Naming convention: AG_<assay>_raw_score is the predicted ALT-REF effect
score for that assay (LFC for RNA; log2-sum DIFF for chromatin marks);
AG_<assay>_quantile_score is AG's percentile rank against a precomputed
genome-wide null distribution in [-1, +1]. Both columns emitted directly
by alphagenome.models.variant_scorers.tidy_scores(). See main_text.md
Methods for the full description.

NA semantics:
- RNA cols in Supp 4: NA if (a) variant not in AG RNA output, (b) MAGE-top
  gene's ENSG isn't in AG's output for this variant (GENCODE v38 vs v46
  mismatch), or (c) variant has no MAGE top gene (n_cis_genes_MAGE = 0).
- RNA cols in Supp 5: NA if AG didn't return that (variant, gene) pair.
- Chromatin cols: NA if variant absent from AG chromatin master.

Usage:
    python scripts/fig3_polymorphic_TE_eQTLs/merge_AG_scores_into_supp_tables.py
"""
import argparse
from pathlib import Path

import pandas as pd

p = argparse.ArgumentParser()
p.add_argument('--variant_class', default='INS',
               help='INS (default; merges into Supp 4 + 5) or DEL (merges into Supp 7 + 8). '
                    'Changes default input + output paths to the corresponding arm.')
p.add_argument('--rna-master',       default=None,
               help='AG GeneMaskLFCScorer master CSV; default depends on --variant_class')
p.add_argument('--chromatin-master', default=None,
               help='AG chromatin master CSV; default depends on --variant_class')
p.add_argument('--top-gene-supp',    default=None,
               help='Top-gene-per-variant supp table; default Supp 4 (INS) or Supp 7 (DEL)')
p.add_argument('--all-pairs-supp',   default=None,
               help='All-variant-gene-pairs supp table; default Supp 5 (INS) or Supp 8 (DEL)')
# Back-compat aliases (drop in a later cleanup):
p.add_argument('--supp4', default=None, help='(legacy alias for --top-gene-supp)')
p.add_argument('--supp5', default=None, help='(legacy alias for --all-pairs-supp)')
args = p.parse_args()

args.variant_class = args.variant_class.upper()
assert args.variant_class in ('INS', 'DEL'), '--variant_class must be INS or DEL'

# Resolve defaults based on variant_class. The DEL arm scores go under
# results/AG_*_polymorphic_TE_DEL/ and merge into Supp 7 + 8.
if args.variant_class == 'INS':
    if args.rna_master is None:       args.rna_master       = 'results/AG_LFC_polymorphic_TE/master_ranked_GM12878.csv'
    if args.chromatin_master is None: args.chromatin_master = 'results/AG_chromatin_polymorphic_TE/master_chromatin_GM12878.csv'
    if args.top_gene_supp is None:    args.top_gene_supp    = args.supp4 or 'supptables/supp_table_4_top_gene_per_variant.tsv'
    if args.all_pairs_supp is None:   args.all_pairs_supp   = args.supp5 or 'supptables/supp_table_5_all_variant_gene_pairs.tsv'
else:  # DEL
    if args.rna_master is None:       args.rna_master       = 'results/AG_LFC_polymorphic_TE_DEL/master_ranked_GM12878.csv'
    if args.chromatin_master is None: args.chromatin_master = 'results/AG_chromatin_polymorphic_TE_DEL/master_chromatin_GM12878.csv'
    if args.top_gene_supp is None:    args.top_gene_supp    = args.supp4 or 'supptables/supp_table_7_top_gene_per_variant_DEL.tsv'
    if args.all_pairs_supp is None:   args.all_pairs_supp   = args.supp5 or 'supptables/supp_table_8_all_variant_gene_pairs_DEL.tsv'
print(f'Variant class: {args.variant_class}')
print(f'  RNA master:   {args.rna_master}')
print(f'  Chrom master: {args.chromatin_master}')
print(f'  Top-gene supp -> {args.top_gene_supp}')
print(f'  All-pairs supp -> {args.all_pairs_supp}')

# RNA columns to add (one pair) + chromatin columns to add (3 pairs).
NEW_AG_COLS = [
    'AG_RNA_raw_score',      'AG_RNA_quantile_score',
    'AG_H3K27ac_raw_score',  'AG_H3K27ac_quantile_score',
    'AG_H3K4me1_raw_score',  'AG_H3K4me1_quantile_score',
    'AG_ATAC_raw_score',     'AG_ATAC_quantile_score',
]
# Two summary columns appended after the score columns
NEW_SUMMARY_COLS = ['AG_RNA_concordance_with_MAGE', 'AG_chromatin_direction']
# Legacy column names that may exist from a previous run; drop them before
# the new merge so we don't leave stale duplicates.
LEGACY_AG_COLS = ['AG_raw_score', 'AG_quantile_score',
                  'AG_RNA_concordance_MAGE']  # earlier candidate name

# Confidence threshold: |quantile_score| >= 0.99 marks "high-confidence"
# predictions, following Avsec et al 2026's "approximately 99th percentile
# of common variants" convention (AG SDK quantile_score is scaled to [-1, +1]
# such that saturation toward ±1.0 represents extreme-tail predictions).
QTL_THRESHOLD = 0.99

# === RNA aggregation: mean across STRAND-MATCHED STRANDED tracks ===========
# Per Liu et al 2026 convention, we use the stranded total + stranded polyA
# RNA-seq tracks only (two strand-specific tracks for the gene's strand),
# dropping the unstranded polyA track (track_strand == '.') from both viz and
# scoring. This:
#   - matches the precedent paper's visualization
#   - avoids the protocol-divergence problem where the unstranded track
#     occasionally disagreed with the two stranded counterparts at small inserts
#   - simplifies the aggregation to a clean "mean of two stranded tracks"
rna = pd.read_csv(args.rna_master)
rna['gene_id_uv'] = rna['gene_id'].str.replace(r'\.\d+$', '', regex=True)
n_total = len(rna)
rna = rna[rna['track_strand'].isin(['+', '-'])].copy()
print(f'RNA master: {n_total:,} rows -> {len(rna):,} after dropping unstranded polyA '
      f'({n_total - len(rna):,} unstranded rows removed)')
rna_agg = (rna.groupby(['variant_id', 'gene_id_uv'], as_index=False)
              .agg(AG_RNA_raw_score=('raw_score', 'mean'),
                   AG_RNA_quantile_score=('quantile_score', 'mean'),
                   _n_tracks=('raw_score', 'size')))
print(f'Aggregated: {len(rna_agg):,} (variant, gene_uv) pairs '
      f'(mean tracks={rna_agg["_n_tracks"].mean():.2f}, expected 2)')

# === Chromatin pivot: one column pair per mark, one row per variant ========
# The master has 3 rows per variant (ATAC + H3K27ac + H3K4me1).
chrom = pd.read_csv(args.chromatin_master)

def mark_label(row):
    """Map AG output_type/histone_mark to our column suffix."""
    if row['output_type'] == 'ATAC':         return 'ATAC'
    if row['output_type'] == 'CHIP_HISTONE': return row['histone_mark']
    return None

chrom['mark'] = chrom.apply(mark_label, axis=1)
chrom_raw = chrom.pivot_table(index='variant_id', columns='mark',
                              values='raw_score', aggfunc='first')
chrom_qtl = chrom.pivot_table(index='variant_id', columns='mark',
                              values='quantile_score', aggfunc='first')
chrom_raw.columns = [f'AG_{m}_raw_score'      for m in chrom_raw.columns]
chrom_qtl.columns = [f'AG_{m}_quantile_score' for m in chrom_qtl.columns]
chrom_wide = chrom_raw.join(chrom_qtl).reset_index()
print(f'Chromatin master: {len(chrom):,} rows -> {len(chrom_wide):,} variants × '
      f'{len(chrom_raw.columns)} marks (= {len(chrom_raw.columns) + len(chrom_qtl.columns)} new cols)')

def classify_rna_concordance(row, beta_col='beta_MAGE'):
    """
    AG_RNA prediction vs MAGE β direction. The `_high_confidence` suffix is
    appended when |AG_RNA_quantile_score| >= 0.99 (Avsec 2026 ~99th-percentile
    convention) — the default label without suffix is still a real match /
    mismatch, just not at the strongest-possible-confidence tier.

    Returns one of:
      match_high_confidence    : signs match, |AG_RNA_quantile_score| >= 0.99
      match                    : signs match, |AG_RNA_quantile_score| <  0.99
      mismatch_high_confidence : signs differ, |AG_RNA_quantile_score| >= 0.99
      mismatch                 : signs differ, |AG_RNA_quantile_score| <  0.99
      AG_RNA_not_available     : AG_RNA_raw_score is NA
      MAGE_beta_not_available  : MAGE β is NA
    """
    mage = row.get(beta_col); ag_raw = row.get('AG_RNA_raw_score')
    ag_qtl = row.get('AG_RNA_quantile_score')
    if pd.isna(ag_raw):  return 'AG_RNA_not_available'
    if pd.isna(mage):    return 'MAGE_beta_not_available'
    same_sign = (mage > 0) == (ag_raw > 0)
    high_conf = abs(ag_qtl) >= QTL_THRESHOLD if pd.notna(ag_qtl) else False
    if same_sign and high_conf:     return 'match_high_confidence'
    if same_sign and not high_conf: return 'match'
    if not same_sign and high_conf: return 'mismatch_high_confidence'
    return 'mismatch'


def classify_chromatin_direction(row):
    """
    Cross-mark summary of the three chromatin marks (H3K27ac + H3K4me1 + ATAC).
    A consistent gain/loss call requires all three marks to agree on sign; the
    `_high_confidence` suffix requires all three |quantile_score| >= 0.99.

    Returns one of:
      consistent_gain_high_confidence
      consistent_gain
      consistent_loss_high_confidence
      consistent_loss
      mixed_direction
    """
    raws = [row['AG_H3K27ac_raw_score'], row['AG_H3K4me1_raw_score'],
            row['AG_ATAC_raw_score']]
    qtls = [row['AG_H3K27ac_quantile_score'], row['AG_H3K4me1_quantile_score'],
            row['AG_ATAC_quantile_score']]
    if any(pd.isna(r) for r in raws): return 'chromatin_not_available'
    all_positive = all(r > 0 for r in raws)
    all_negative = all(r < 0 for r in raws)
    if not (all_positive or all_negative): return 'mixed_direction'
    high_conf = all(pd.notna(q) and abs(q) >= QTL_THRESHOLD for q in qtls)
    if all_positive and high_conf:     return 'consistent_gain_high_confidence'
    if all_positive and not high_conf: return 'consistent_gain'
    if all_negative and high_conf:     return 'consistent_loss_high_confidence'
    return 'consistent_loss'


def merge_supp(path, gene_col, label):
    df = pd.read_csv(path, sep='\t')
    n_before = len(df)
    df = df.drop(columns=[c for c in LEGACY_AG_COLS + NEW_AG_COLS + NEW_SUMMARY_COLS
                          if c in df.columns])

    # RNA join key: (variant_id, gene_id_uv)
    df['_gene_uv'] = df[gene_col].fillna('').str.replace(r'\.\d+$', '', regex=True)
    df = df.merge(
        rna_agg.rename(columns={'gene_id_uv': '_gene_uv'})[
            ['variant_id', '_gene_uv', 'AG_RNA_raw_score', 'AG_RNA_quantile_score']],
        on=['variant_id', '_gene_uv'], how='left'
    ).drop(columns=['_gene_uv'])

    # Chromatin join key: variant_id (one row of chromatin per variant)
    df = df.merge(chrom_wide, on='variant_id', how='left')

    # Summary calls
    df['AG_RNA_concordance_with_MAGE'] = df.apply(classify_rna_concordance, axis=1)
    df['AG_chromatin_direction']       = df.apply(classify_chromatin_direction, axis=1)

    assert len(df) == n_before, f'{label}: row count changed {n_before} -> {len(df)}'

    rna_with   = df['AG_RNA_raw_score'].notna().sum()
    chrom_with = df['AG_H3K27ac_raw_score'].notna().sum()
    print(f'{label}: {n_before:,} rows')
    print(f'  AG_RNA filled:   {rna_with:>5,} ({100*rna_with/n_before:.1f}%)')
    print(f'  AG_chromatin filled: {chrom_with:>5,} ({100*chrom_with/n_before:.1f}%)')
    print(f'  AG_RNA_concordance_with_MAGE breakdown:')
    for v, n in df['AG_RNA_concordance_with_MAGE'].value_counts().items():
        print(f'    {v:<35s} {n:>6,}')
    print(f'  AG_chromatin_direction breakdown:')
    for v, n in df['AG_chromatin_direction'].value_counts().items():
        print(f'    {v:<35s} {n:>6,}')

    # Reorder: existing columns first, then the 8 score columns, then summaries
    score_present = [c for c in NEW_AG_COLS      if c in df.columns]
    sum_present   = [c for c in NEW_SUMMARY_COLS if c in df.columns]
    other = [c for c in df.columns if c not in score_present + sum_present]
    df = df[other + score_present + sum_present]

    df.to_csv(path, sep='\t', index=False)
    print(f'  wrote {path}')

top_gene_label  = 'Supp 4' if args.variant_class == 'INS' else 'Supp 7'
all_pairs_label = 'Supp 5' if args.variant_class == 'INS' else 'Supp 8'
merge_supp(args.top_gene_supp,  gene_col='top_gene_ensg_MAGE', label=top_gene_label)
merge_supp(args.all_pairs_supp, gene_col='gene_ensg_MAGE',     label=all_pairs_label)

print('\nDone.')
