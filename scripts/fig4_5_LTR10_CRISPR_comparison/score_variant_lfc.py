"""
Score TE variants (INSERTIONS or DELETIONS) with AlphaGenome's GeneMaskLFCScorer.
Accepts two encodings of the (CHROM, POS, REF, ALT) input row, both valid:

  (1) VCF-biallelic form — used for Schloissnig polymorphic TEs.
      INS: REF = anchor base, ALT = anchor + inserted sequence.
      DEL: REF = anchor + deleted sequence, ALT = anchor only.

  (2) CRISPR synthetic-deletion form — used for CRISPRi-validated CREs
      (e.g. data/LTR10_variants.tab).
      REF = full element sequence, ALT = '.' or '' (no anchor).
      Normalized internally to alternate_bases=''.

**1 Mb scoring window:** built via
    interval = variant.reference_interval.resize(SEQUENCE_LENGTH_1MB)
matching AG's variant_scoring_ui tutorial exactly
(https://www.alphagenomedocs.com/colabs/variant_scoring_ui.html). The
tutorial describes this as "the input interval is derived from the variant
(centered on it)". For an SNV/INS (len(REF)=1) this centers the window on
POS; for a DEL (len(REF)>1) the window centers on the midpoint of the
deleted span, i.e. centered on the affected region rather than the anchor.
This is AG's documented and intended behavior — we follow it across both
encodings (VCF-biallelic and CRISPR synthetic) rather than fork centering.

Output: one CSV per variant in OUTPUT_DIR/, plus a master ranked CSV across
all variants. Each CSV has one row per (variant, gene, RNA-seq track) tuple
with two effect columns:
  - raw_score      = AG GeneMaskLFCScorer effect magnitude (= LFC, Liu Fig 2A)
  - quantile_score = AG-emitted quantile of raw_score against a precomputed
                     genome-wide null distribution (interpretable as a
                     confidence/extremeness score in [-1, +1]; |q|→1 means
                     the predicted effect ranks at the extreme of the null).
Both columns are auto-emitted by AG's variant_scorers.tidy_scores() — the
quantile_score comes from the AnnData layers['quantiles'] tensor.

Following Liu et al 2026 (RHD bioRxiv) Fig 2A for raw_score:
    LFC = log(mean(ALT_RNA_seq over exon mask) + 1e-3)
        - log(mean(REF_RNA_seq over exon mask) + 1e-3)
GeneMaskLFCScorer implements this exactly; the small pseudocount is
AG-internal and does not need to be set by the caller.

Usage:
    python scripts/fig4_5_LTR10_CRISPR_comparison/score_variant_lfc.py \\
        --variants data/fig3_tier1_variants.tab \\
        --cell_line GM12878 \\
        --output_dir results/AG_LFC_polymorphic_TE
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
p.add_argument('--allow-mnv',  action='store_true',
               help='Allow same-length REF/ALT (SNP / MNV). Default is to reject — useful '
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

# Pick the recommended GeneMaskLFCScorer for RNA_SEQ (matches Liu Fig 2A formula).
# Use direct class+enum comparison rather than parsing __name__/__str__ — the
# string forms can drift between AG SDK versions; isinstance / Enum equality
# is stable.
selected = [s for s in variant_scorers.RECOMMENDED_VARIANT_SCORERS.values()
            if isinstance(s, variant_scorers.GeneMaskLFCScorer)
            and s.requested_output == OutputType.RNA_SEQ]
if not selected:
    sys.exit("ERROR: no RNA_SEQ GeneMaskLFCScorer in RECOMMENDED_VARIANT_SCORERS")
print(f"Using {len(selected)} RNA_SEQ scorer(s)")

variant_df = pd.read_table(args.variants, sep='\t')
if args.test: variant_df = variant_df.head(1)
print(f"Loaded {len(variant_df)} variants from {args.variants}")

# Acceptable bases for a literal sequence ALT/REF. 'N' is tolerated (long-read
# assemblies sometimes carry unknown bases). Anything else (symbolic alleles,
# multi-allelic comma lists, non-string NaN) is rejected before AG is called.
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

    # REF/ALT validation. Accept two forms:
    #   - VCF-biallelic: REF and ALT both ACGT(N) sequences, ALT not symbolic.
    #   - CRISPR synthetic-DEL: ALT == '' or '.', meaning "delete REF entirely
    #     with no anchor base". Normalized below to alternate_bases=''.
    # Reject in either form: NaN, '*', symbolic <INS>/<DEL>, multi-allelic ALT,
    # non-ACGT bases in a non-empty allele.
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

    # Net length change: positive for INS, negative for DEL. Reject SNPs /
    # equal-length MNVs by default (AG can score them, but this script is
    # built around polymorphic-TE INS/DEL — fail loudly rather than silently
    # rescore). Pass --allow-mnv to score same-length substitutions like the
    # Fig 6 AP1-perturbation alleles.
    delta = len(var_alt) - len(var_ref)
    if delta == 0 and not args.allow_mnv:
        print(f"[{i+1}/{len(variant_df)}] SKIP {vid}: REF and ALT same length (SNP/MNV); pass --allow-mnv to override"); continue

    vtype = "MNV" if delta == 0 else ("INS" if delta > 0 else "DEL")
    print(f"[{i+1}/{len(variant_df)}] {vid}  {chrom}:{pos} {vtype} {abs(delta)} bp")
    variant = genome.Variant(
        chromosome=chrom, position=pos,
        reference_bases=var_ref, alternate_bases=var_alt,
        name=vid,
    )
    # AG receptive field 1 Mb, centered on the variant. NOTE: score_variant
    # does NOT accept an `ontology_terms` argument (unlike predict_variant);
    # it returns scores across every biosample the scorer covers and we filter
    # by biosample_name post-hoc. This is the canonical AG-SDK convention for
    # variant scoring, also used by scripts/wip/batch_score_elements.py.
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
