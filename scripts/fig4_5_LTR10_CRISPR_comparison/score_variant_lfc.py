"""
Score TE variants (INS, DEL, CRISPR synthetic-DEL, or same-length substitution
variants with --allow-sub) with AlphaGenome's GeneMaskLFCScorer (RNA log2FC,
Liu 2026 Fig 2A). Input encodings: VCF-biallelic (INS: REF=anchor, ALT=anchor+
insert; DEL: the reverse) and CRISPR synthetic-DEL (REF=element, ALT='.'/'' ->
normalized to ''). The 1 Mb window = variant.reference_interval.resize(1Mb),
centered on the variant (deleted-span midpoint for DELs). Outputs one CSV per
variant + a master ranked CSV; rows are (variant, gene, RNA track) with raw_score
(LFC) and quantile_score (AG percentile vs null, [-1,+1]).
Shared script: used by Fig 3 (eQTL), Fig 4/5 (LTR10 CRISPR), and Fig 6.

Example usage (polymorphic TE variants):
python scripts/fig4_5_LTR10_CRISPR_comparison/score_variant_lfc.py \
    --variants data/fig3_all_testable_variants.tab \
    --output_dir results/AG_LFC_polymorphic_TE

Example usage (CRISPR-validated LTR10 enhancers, HCT116):
python scripts/fig4_5_LTR10_CRISPR_comparison/score_variant_lfc.py \
    --variants data/LTR10_variants.tab \
    --cell_line HCT116 --ontology EFO:0002824 \
    --output_dir results/AG_LFC_LTR10_CRISPR
"""
import argparse, os, sys
from pathlib import Path

import pandas as pd
import alphagenome
from alphagenome.data import genome
from alphagenome.models import dna_client, variant_scorers
from alphagenome.models.dna_output import OutputType

p = argparse.ArgumentParser()
p.add_argument('--variants',   required=True, help='tab file: ID CHROM POS REF ALT')
p.add_argument('--cell_line',  default='GM12878', help='used only in output filenames')
p.add_argument('--ontology',   default='EFO:0002784',
               help='Metadata-only — printed for record-keeping. AG score_variant() does NOT '
                    'accept ontology_terms (unlike predict_variant); biosample filtering happens '
                    'post-hoc by --cell_line.')
p.add_argument('--output_dir', default='results/AG_LFC_polymorphic_TE')
p.add_argument('--api_key',    default='scripts/my_api_key.txt')
p.add_argument('--test',       action='store_true', help='Score only the first variant')
p.add_argument('--allow-sub',  action='store_true',
               help='Allow same-length REF/ALT (SNP / substitution variant). Default is to reject — useful '
                    'when running the Fig 6 AP1-perturbation alleles where REF and ALT are '
                    'both 2,358 bp.')
args = p.parse_args()

with open(args.api_key) as f:
    api_key = f.read().strip()
out_dir = Path(args.output_dir); out_dir.mkdir(parents=True, exist_ok=True)
print(f"AlphaGenome version: {alphagenome.__version__}")
print(f"Cell line: {args.cell_line}  |  Ontology: {args.ontology}")
print("Initializing model...")
dna_model = dna_client.create(api_key)

# Pick the recommended GeneMaskLFCScorer for RNA_SEQ (Liu Fig 2A). Direct
# class+enum comparison is stable across AG SDK versions.
selected = [s for s in variant_scorers.RECOMMENDED_VARIANT_SCORERS.values()
            if isinstance(s, variant_scorers.GeneMaskLFCScorer)
            and s.requested_output == OutputType.RNA_SEQ]
if not selected:
    sys.exit("ERROR: no RNA_SEQ GeneMaskLFCScorer in RECOMMENDED_VARIANT_SCORERS")
print(f"Using {len(selected)} RNA_SEQ scorer(s)")

variant_df = pd.read_table(args.variants, sep='\t')
if args.test: variant_df = variant_df.head(1)
print(f"Loaded {len(variant_df)} variants from {args.variants}")

# Acceptable bases for a literal ALT/REF ('N' tolerated). Symbolic / multi-allelic /
# non-string alleles are rejected before calling AG.
_ACGT = set('ACGTN')

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

    # REF/ALT validation: VCF-biallelic (both ACGT(N), ALT not symbolic) or CRISPR
    # synthetic-DEL (ALT '' or '.'). Rejects NaN, '*', <INS>/<DEL>, multi-allelic, non-ACGT.
    if not isinstance(var_alt, str) or not isinstance(var_ref, str):
        print(f"[{i+1}/{len(variant_df)}] SKIP {vid}: ALT/REF is not a string"); continue
    if var_alt == '*' or var_alt.startswith('<') or ',' in var_alt:
        print(f"[{i+1}/{len(variant_df)}] SKIP {vid}: ALT symbolic / multi-allelic ({var_alt!r})"); continue
    if var_ref.startswith('<'):
        print(f"[{i+1}/{len(variant_df)}] SKIP {vid}: REF symbolic ({var_ref!r})"); continue
    # Normalize CRISPR synthetic-DEL: '.' → '' (alternate_bases='' is AG's
    # canonical full-deletion form, used by batch_score_elements.py).
    if var_alt == '.':
        var_alt = ''
    if var_alt and not set(var_alt.upper()) <= _ACGT:
        print(f"[{i+1}/{len(variant_df)}] SKIP {vid}: non-ACGT bases in ALT"); continue
    if not set(var_ref.upper()) <= _ACGT:
        print(f"[{i+1}/{len(variant_df)}] SKIP {vid}: non-ACGT bases in REF"); continue
    var_ref = var_ref.upper(); var_alt = var_alt.upper()

    # Net length change: + for INS, - for DEL. SNPs / equal-length substitution
    # variants rejected by default; --allow-sub to score them (Fig 6 AP1 alleles).
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
    # 1 Mb window centered on the variant. score_variant does NOT take ontology_terms;
    # it scores every biosample the scorer covers, so we filter by biosample_name after.
    try:
        scores = dna_model.score_variant(
            interval=variant.reference_interval.resize(dna_client.SEQUENCE_LENGTH_1MB),
            variant=variant,
            variant_scorers=selected,
            organism=dna_client.Organism.HOMO_SAPIENS,
        )
    except (RuntimeError, ValueError, ConnectionError, TimeoutError) as e:
        # Narrowed from `except Exception`: keep batch resumable through
        # network glitches / AG transient errors, but let unexpected bugs
        # (AttributeError, ImportError, KeyError on a stable schema, ...)
        # propagate so we notice them instead of silently skipping variants.
        print(f"    ! AG score_variant failed: {e}; skipping"); continue

    # tidy_scores can return None when AG had no scoreable genes for this variant
    # (e.g., gene-deserts, masked-window edges). Skip and move on rather than crash.
    effects = variant_scorers.tidy_scores(scores, match_gene_strand=True)
    if effects is None or len(effects) == 0:
        print(f"    ! AG returned no scoreable genes; skipping"); continue
    df = effects[effects['biosample_name'] == args.cell_line].copy()
    if len(df) == 0:
        print(f"    ! AG returned no rows for biosample={args.cell_line}; skipping"); continue
    df['variant_id'] = vid
    df = df.drop(columns=[c for c in ['scored_interval'] if c in df.columns])
    df = df.sort_values('raw_score').reset_index(drop=True)
    df.to_csv(csv_path, index=False)
    print(f"    -> {csv_path}  ({len(df)} (variant, gene) rows)")
    master_rows.append(df)

# Master ranked CSV (negative scores = downregulation, à la Liu et al)
if master_rows:
    master = pd.concat(master_rows, ignore_index=True).sort_values('raw_score')
    master_path = out_dir / f'master_ranked_{args.cell_line}.csv'
    master.to_csv(master_path, index=False)
    print(f"\nMaster ranked: {master_path}  ({len(master)} rows from {len(master_rows)} variants)")
