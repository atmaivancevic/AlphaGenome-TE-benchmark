"""
Score TE variants (INS, DEL, or same-length substitution variants with
--allow-sub) with AlphaGenome's chromatin scorer (CenterMaskScorer +
DIFF_LOG2_SUM). Sister to score_variant_lfc.py (the RNA arm). Per (variant, mark)
it reports raw_score (log2 ALT/REF over the center mask) and quantile_score (AG
percentile vs a null, in [-1, +1]) for the active-enhancer marks H3K27ac,
H3K4me1, and ATAC. The 1 Mb window is built as
variant.reference_interval.resize(1Mb), centered on the variant (mask at the
deleted-span midpoint for DELs). Mask widths: 501 bp ATAC, 2001 bp ChIP.
Shared script: also used by Fig 4/5 (LTR10 CRISPR) and Fig 6 (--include-tf).

Example usage (polymorphic TE variants):
python scripts/fig3_polymorphic_TE_eQTLs/score_variant_chromatin.py \
    --variants data/fig3_all_testable_variants.tab \
    --output_dir results/AG_chromatin_polymorphic_TE

Example usage (LTR10 AP1 perturbation substitution variants + TF tracks, Fig 6):
python scripts/fig3_polymorphic_TE_eQTLs/score_variant_chromatin.py \
    --variants data/LTR10_ATG12_AP1_perturbations.tab \
    --allow-sub --include-tf FOSL1 JUND \
    --output_dir results/AG_chromatin_LTR10_ATG12_AP1
"""
import argparse
import sys
from pathlib import Path

import pandas as pd
import alphagenome
from alphagenome.data import genome
from alphagenome.models import dna_client, variant_scorers
from alphagenome.models.dna_output import OutputType

# Marks of interest — filter applied to AG output after scoring
TARGET_HISTONE_MARKS = ['H3K27ac', 'H3K4me1']
# (ATAC has a single track so no mark filter needed.)

# Acceptable bases for a literal sequence ALT/REF. 'N' is tolerated.
_ACGT = set('ACGTN')

p = argparse.ArgumentParser()
p.add_argument('--variants',   required=True, help='tab file: ID CHROM POS REF ALT')
p.add_argument('--cell_line',  default='GM12878', help='used only in output filenames')
p.add_argument('--ontology',   default='EFO:0002784',
               help='Metadata-only — printed for record-keeping. AG score_variant() does NOT '
                    'accept ontology_terms (unlike predict_variant); biosample filtering happens '
                    'post-hoc by --cell_line.')
p.add_argument('--output_dir', default='results/AG_chromatin_polymorphic_TE')
p.add_argument('--api_key',    default='scripts/my_api_key.txt')
p.add_argument('--test',       action='store_true', help='Score only the first variant')
p.add_argument('--allow-sub',  action='store_true',
               help='Allow same-length REF/ALT (SNP / substitution variant). Default is to reject — useful '
                    'when running the Fig 6 AP1-perturbation alleles where REF and ALT are '
                    'both 2,358 bp.')
p.add_argument('--include-tf', nargs='*', default=None,
               help='Optional list of TF names to also score via the CHIP_TF CenterMaskScorer '
                    '(default: histone + ATAC only). Example: --include-tf FOSL1 FOS JUN. '
                    'Filters AG-returned CHIP_TF tracks by track_name substring.')
args = p.parse_args()

with open(args.api_key) as f:
    api_key = f.read().strip()
out_dir = Path(args.output_dir); out_dir.mkdir(parents=True, exist_ok=True)
print(f"AlphaGenome version: {alphagenome.__version__}")
print(f"Cell line: {args.cell_line}  |  Ontology: {args.ontology}")
print(f"Histone marks of interest: {TARGET_HISTONE_MARKS}  +  ATAC")
print("Initializing model...")
dna_model = dna_client.create(api_key)

# Pick AG's recommended CenterMaskScorers for CHIP_HISTONE + ATAC, and CHIP_TF
# when --include-tf is set (Fig 6 FOSL1/JUND). Default widths 2001/501 bp,
# DIFF_LOG2_SUM aggregation = log2 fold-change at the variant locus.
TARGET_OUTPUTS = {OutputType.CHIP_HISTONE, OutputType.ATAC}
if args.include_tf:
    TARGET_OUTPUTS.add(OutputType.CHIP_TF)
selected = [
    s for s in variant_scorers.RECOMMENDED_VARIANT_SCORERS.values()
    if isinstance(s, variant_scorers.CenterMaskScorer)
    and s.requested_output in TARGET_OUTPUTS
    and s.aggregation_type == variant_scorers.AggregationType.DIFF_LOG2_SUM
]
expected = 3 if args.include_tf else 2
if len(selected) != expected:
    sys.exit(f"ERROR: expected {expected} scorers, got {len(selected)}")
print(f"Using {len(selected)} scorers:")
for s in selected:
    print(f"  {type(s).__name__}  output={s.requested_output.name}  "
          f"width={s.width}  agg={s.aggregation_type.name}")

variant_df = pd.read_table(args.variants, sep='\t')
if args.test: variant_df = variant_df.head(1)
print(f"Loaded {len(variant_df)} variants from {args.variants}")

master_rows = []
for i, row in variant_df.iterrows():
    vid = row['ID']
    chrom = row['CHROM']; pos = int(row['POS'])
    var_ref = row['REF']; var_alt = row['ALT']
    csv_path = out_dir / f'{vid}_{args.cell_line}.csv'

    # Checkpoint: skip if already scored
    if csv_path.exists():
        print(f"[{i+1}/{len(variant_df)}] SKIP {vid} (already scored)")
        prev = pd.read_csv(csv_path)
        master_rows.append(prev); continue

    # REF/ALT validation. Two accepted forms: VCF-biallelic (both ACGT(N), ALT not
    # symbolic) and CRISPR synthetic-DEL (ALT == '' or '.', normalized to '').
    if not isinstance(var_alt, str) or not isinstance(var_ref, str):
        print(f"[{i+1}/{len(variant_df)}] SKIP {vid}: ALT/REF is not a string"); continue
    if var_alt == '*' or var_alt.startswith('<') or ',' in var_alt:
        print(f"[{i+1}/{len(variant_df)}] SKIP {vid}: ALT symbolic / multi-allelic ({var_alt!r})"); continue
    if var_ref.startswith('<'):
        print(f"[{i+1}/{len(variant_df)}] SKIP {vid}: REF symbolic ({var_ref!r})"); continue
    if var_alt == '.':
        var_alt = ''
    if var_alt and not set(var_alt.upper()) <= _ACGT:
        print(f"[{i+1}/{len(variant_df)}] SKIP {vid}: non-ACGT bases in ALT"); continue
    if not set(var_ref.upper()) <= _ACGT:
        print(f"[{i+1}/{len(variant_df)}] SKIP {vid}: non-ACGT bases in REF"); continue
    var_ref = var_ref.upper(); var_alt = var_alt.upper()

    # Net length change: + for INS, - for DEL. SNPs / equal-length substitution
    # variants rejected by default (this script targets polymorphic-TE INS/DEL);
    # --allow-sub to score them.
    delta = len(var_alt) - len(var_ref)
    if delta == 0 and not args.allow_sub:
        print(f"[{i+1}/{len(variant_df)}] SKIP {vid}: REF and ALT same length (SNP/SUB); pass --allow-sub to override"); continue

    vtype = "SUB" if delta == 0 else ("INS" if delta > 0 else "DEL")
    print(f"[{i+1}/{len(variant_df)}] {vid}  {chrom}:{pos} {vtype} {abs(delta)} bp")
    variant = genome.Variant(
        chromosome=chrom, position=pos,
        reference_bases=var_ref, alternate_bases=var_alt,
        name=vid,
    )

    # score_variant does NOT take ontology_terms (unlike predict_variant); it scores
    # every biosample the scorer covers, so we filter by biosample_name afterwards.
    try:
        scores = dna_model.score_variant(
            interval=variant.reference_interval.resize(dna_client.SEQUENCE_LENGTH_1MB),
            variant=variant,
            variant_scorers=selected,
            organism=dna_client.Organism.HOMO_SAPIENS,
        )
    except (RuntimeError, ValueError, ConnectionError, TimeoutError) as e:
        # Narrowed from `except Exception`: keep the batch resumable through network /
        # transient AG errors, but let real bugs propagate.
        print(f"    ! AG score_variant failed: {e}; skipping"); continue

    effects = variant_scorers.tidy_scores(scores, match_gene_strand=True)
    if effects is None or len(effects) == 0:
        print(f"    ! AG returned no scoreable rows; skipping"); continue

    # Filter to target biosample + marks: CHIP_HISTONE by requested histone marks,
    # ATAC passes through, CHIP_TF (with --include-tf) by TF name.
    df = effects[effects['biosample_name'] == args.cell_line].copy()
    if len(df) and 'output_type' in df.columns:
        keep_mask = pd.Series(False, index=df.index)
        # ATAC: keep all
        keep_mask |= (df['output_type'] == 'ATAC')
        # CHIP_HISTONE: filter by histone_mark
        if 'histone_mark' in df.columns:
            keep_mask |= (
                (df['output_type'] == 'CHIP_HISTONE') &
                df['histone_mark'].isin(TARGET_HISTONE_MARKS)
            )
        # CHIP_TF: filter by TF name (substring match on track_name)
        if args.include_tf and 'track_name' in df.columns:
            tf_pattern = '|'.join(args.include_tf)
            keep_mask |= (
                (df['output_type'] == 'CHIP_TF') &
                df['track_name'].astype(str).str.contains(tf_pattern, case=False, na=False)
            )
        df = df[keep_mask].copy()

    if len(df) == 0:
        print(f"    ! No rows after biosample+mark filter; skipping"); continue

    df['variant_id'] = vid
    df = df.drop(columns=[c for c in ['scored_interval'] if c in df.columns])
    df = df.sort_values('raw_score').reset_index(drop=True)
    df.to_csv(csv_path, index=False)

    # Report headline rows (defensive: histone_mark / output_type can be NaN
    # for ATAC rows — coerce to strings before formatting)
    print(f"    -> {csv_path}  ({len(df)} rows)")
    for _, r in df.iterrows():
        ot = str(r.get('output_type', '?'))
        mk = '' if pd.isna(r.get('histone_mark', None)) else str(r.get('histone_mark', ''))
        print(f"       {ot:15s} {mk:10s} raw={r['raw_score']:+.4f}  qtl={r['quantile_score']:+.4f}")

    master_rows.append(df)

# Master CSV — concatenated across variants
if master_rows:
    master = pd.concat(master_rows, ignore_index=True).sort_values(['variant_id', 'raw_score'])
    master_path = out_dir / f'master_chromatin_{args.cell_line}.csv'
    master.to_csv(master_path, index=False)
    print(f"\nMaster: {master_path}  ({len(master)} rows from {len(master_rows)} variants)")
