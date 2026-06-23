"""
Build Supp Table 3 (polymorphic TE insertion catalog) by joining the filtered
solo TE INS VCF (output of filter_solo_TE_INS_vcf.py, hg38-anchored, SVAN-
annotated, no GTs) with the SVIM-asm hg38 genotyped BCF (per-sample GTs, no
SVAN).

Both files are matched on variant ID (`SvimAsm########`). The genotyped BCF
already provides AC, AN, AF, MAF, NS, F_MISSING, AC_Hom, AC_Het, AC_Hemi, HWE,
ExcHet in INFO — no need to recompute. Per-sample iteration is used only to
build the carrier sample-ID lists (carriers_hom, carriers_het) and the flagship
sample's GT (NA12878 by default).

Multi-allelic sites (multiple comma-separated ALTs at the same position) are
skipped: the genotyped BCF still contains some after the `bcftools merge` step
in the Schloissnig pipeline, and packing N alleles into one Supp Table row makes
per-variant stats ambiguous. Affects ~378 of 29,839 solo MEI records (~1.3%).
A future version will split these via `bcftools norm -m -` so each allele
becomes its own row.

Output is written as a new tab in the SuppTables.xlsx workbook (existing tabs
preserved). Optionally also dumps a TSV alongside.

Usage:
    python scripts/generate_supp_tables/build_supp_table3_polymorphic_TE_insertions.py \\
        --filtered-vcf data/schloissnig_2025/vcf/solo_TE_INS.svim_asm.hg38.SVAN.vcf.gz \\
        --gt-bcf       data/schloissnig_2025/vcf/svim.asm.hg38.bcf \\
        --xlsx         supptables/SuppTables.xlsx \\
        --tsv          supptables/supp_table_3_polymorphic_TE_insertions.tsv
"""
import argparse, gzip, re, subprocess
import pandas as pd

p = argparse.ArgumentParser()
p.add_argument('--filtered-vcf', required=True, help='solo TE INS VCF (output of filter_solo_TE_INS_vcf.py)')
p.add_argument('--gt-bcf', required=True, help='SVIM-asm hg38 genotyped BCF (908 samples)')
p.add_argument('--tsv', required=True, help='Output TSV path (authoritative; xlsx is composed manually at submission)')
p.add_argument('--xlsx', default=None, help='Optional: inject into SuppTables.xlsx (skipped by default; Atma maintains the workbook manually)')
p.add_argument('--sheet-name', default='Suppl Table 3 Polymorphic TE Insertions')
p.add_argument('--na-sample', default='NA12878', help='Flagship sample whose GT is reported as a column')
p.add_argument('--eqtl-testable-tsv', default='supptables/supp_table_4_top_gene_per_variant.tsv',
               help='eQTL results table (Supp Table 4) whose variant_id set defines the '
                    'eqtl_testable_MAGE260 flag. Pass an empty string to leave the flag blank '
                    '(e.g. when building the catalog before the eQTL pipeline has run).')
args = p.parse_args()

# 1. SVAN side: per-variant annotation
svan = {}
with gzip.open(args.filtered_vcf, 'rt') as f:
    for line in f:
        if line.startswith('#'): continue
        cols = line.rstrip('\n').split('\t')
        vid, info = cols[2], cols[7]
        pct = re.search(r'PERC_RESOLVED=([\d.]+)', info)
        fam = re.search(r'FAM_N=([^;]+)', info)
        # SVAN's strand call for the insertion: '+' (sense) / '-' (antisense) / '.'
        # (undetermined). For solo LTR5_Hs insertions strand is often '.' because
        # the LTR's terminal repeats make orientation ambiguous from sequence
        # alone; for short SVA fragments orientation may not be resolvable. We
        # report verbatim from SVAN — documented in Methods + supp legend.
        strand = re.search(r'(?:^|;)STRAND=([^;]+)', info)
        # The length column is `SVLEN`, sourced verbatim from the genotyped BCF's
        # INFO field (see step 2 below). Positive integer for INS records
        # (e.g. 973 for a 973 bp insertion). Column name kept as `SVLEN` rather
        # than `ins_len` so it's directly traceable back to the source VCF/BCF
        # without any renaming step. For all 29,839 solo INS records, SVLEN
        # equals len(ALT) - len(REF) exactly (every record has len(REF)=1, so
        # there are no complex-form records in this set). SVAN's INS_LEN field
        # is anchor-INCLUDED (1 bp larger than SVLEN); we deliberately do not
        # use it.
        svan[vid] = dict(
            variant_id=vid, family=fam.group(1) if fam else None,
            chrom=cols[0], pos=int(cols[1]),
            # SVLEN filled below from genotyped BCF
            strand=strand.group(1) if strand else '.',
            ref=cols[3],
            itype_n='solo',
            filter_svan=cols[6],
            perc_resolved=float(pct.group(1)) if pct else None,
            vcf_source=args.filtered_vcf,
        )
print(f'SVAN solo TE INS: {len(svan):,}')

# 2. Genotyped BCF side: pre-baked stats from INFO + per-sample GTs for carriers
NUM_FIELDS = [
    ('SVLEN', int),
    ('AC', int), ('AN', int), ('AF', float), ('MAF', float),
    ('NS', int), ('F_MISSING', float),
    ('AC_Hom', int), ('AC_Het', int), ('AC_Hemi', int),
    ('HWE', float), ('ExcHet', float),
]

samples = None
matched = 0
multi_allelic_ids = set()
proc = subprocess.Popen(['bcftools', 'view', args.gt_bcf], stdout=subprocess.PIPE, text=True)
for line in proc.stdout:
    if line.startswith('##'): continue
    if line.startswith('#CHROM'):
        samples = line.rstrip('\n').split('\t')[9:]
        if args.na_sample not in samples:
            raise SystemExit(f'{args.na_sample} not in genotyped BCF')
        na_idx = samples.index(args.na_sample)
        continue
    cols = line.rstrip('\n').split('\t')
    rec = svan.get(cols[2])
    if rec is None: continue
    if ',' in cols[4]:  # multi-allelic — skip until we add bcftools norm -m - splitting
        multi_allelic_ids.add(cols[2]); continue
    info = cols[7]
    for k, kind in NUM_FIELDS:
        m = re.search(rf'(?:^|;){k}=([^;]+)', info)
        rec[k] = (kind(m.group(1)) if m else None)
    gt_na = cols[9 + na_idx].split(':')[0]
    rec[f'{args.na_sample}_GT'] = gt_na
    rec[f'{args.na_sample}_alt_dose'] = gt_na.count('1')
    hom, het = [], []
    for s, g in zip(samples, cols[9:]):
        n = g.split(':')[0].count('1')
        if n == 2: hom.append(s)
        elif n == 1: het.append(s)
    rec['n_carriers'] = len(hom) + len(het)
    rec['carriers_hom'] = ','.join(hom)
    rec['carriers_het'] = ','.join(het)
    matched += 1
proc.wait()
if proc.returncode != 0:
    raise SystemExit(f'bcftools view failed (rc={proc.returncode})')

# Drop multi-allelic IDs from the SVAN dict so they don't appear in the output
# with empty GT stats.
for vid in multi_allelic_ids: svan.pop(vid, None)

print(f'Variants joined with genotyped BCF (biallelic): {matched:,}')
print(f'Multi-allelic skipped: {len(multi_allelic_ids):,}')

# Fill missing for variants absent from genotyped BCF (should be zero given
# upstream inspection — kept defensive in case future updates change the set)
for vid, rec in svan.items():
    for k, _ in NUM_FIELDS: rec.setdefault(k, None)
    rec.setdefault(f'{args.na_sample}_GT', '')
    rec.setdefault(f'{args.na_sample}_alt_dose', None)
    rec.setdefault('n_carriers', None)
    rec.setdefault('carriers_hom', '')
    rec.setdefault('carriers_het', '')

# Flag columns linking the full catalog to the Fig 3 eQTL analysis:
#   na12878_00            — NA12878 is hom-ref (matches hg38; the eQTL candidate pool that
#                           mitigates AlphaGenome's GM12878 reference-allele bias).
#   eqtl_testable_MAGE260 — variant actually entered the MAGE-260 cis-eQTL test, i.e. it
#                           appears in Supp Table 4 (NA12878=0/0 + MAF>0.01 + >=5 hom-alt +
#                           >=5 hom-ref carriers in MAGE-260). These are the ~1.3k variants
#                           behind Figure 3; the rest of the catalog is context/resource.
testable_ids = set()
if args.eqtl_testable_tsv:
    try:
        testable_ids = set(pd.read_csv(args.eqtl_testable_tsv, sep='\t',
                                       usecols=['variant_id'])['variant_id'])
        print(f'eQTL-testable variants (from {args.eqtl_testable_tsv}): {len(testable_ids):,}')
    except FileNotFoundError:
        print(f'WARNING: {args.eqtl_testable_tsv} not found; eqtl_testable_MAGE260 left blank')
for vid, rec in svan.items():
    rec['na12878_00'] = (rec.get(f'{args.na_sample}_GT') == '0/0')
    rec['eqtl_testable_MAGE260'] = (vid in testable_ids) if testable_ids else None

cols = ['variant_id','family','chrom','pos','SVLEN','strand','ref',
        'itype_n','filter_svan','perc_resolved',
        'AC','AN','AF','MAF','NS','F_MISSING',
        'AC_Hom','AC_Het','AC_Hemi','HWE','ExcHet',
        f'{args.na_sample}_GT', f'{args.na_sample}_alt_dose',
        'na12878_00','eqtl_testable_MAGE260',
        'n_carriers','carriers_hom','carriers_het','vcf_source']
# Sort: TE family group order (LTR5_Hs > SVA > L1 > Alu > misc), then chrom + pos within family.
family_rank = {'LTR5_Hs': 1, 'SVA': 2, 'L1': 3, 'Alu': 4}
df = pd.DataFrame(list(svan.values()))[cols]
# Build a chrom sort key that orders chr1..chr22 numerically, then chrX, chrY, chrM
def chrom_key(c):
    s = str(c).replace('chr', '')
    return (0, int(s)) if s.isdigit() else (1, {'X':0,'Y':1,'M':2}.get(s, 3))
df['_fam_rank'] = df['family'].map(lambda x: family_rank.get(x, 5))
df['_chr_key']  = df['chrom'].map(chrom_key)
df = df.sort_values(['_fam_rank', '_chr_key', 'pos']).drop(columns=['_fam_rank','_chr_key']).reset_index(drop=True)
print(f'Output rows x cols: {len(df):,} x {len(df.columns)}')

df.to_csv(args.tsv, sep='\t', index=False)
print(f'Wrote TSV: {args.tsv}')

if args.xlsx:
    existing = pd.ExcelFile(args.xlsx)
    sheets = {s: pd.read_excel(args.xlsx, sheet_name=s) for s in existing.sheet_names}
    sheets[args.sheet_name] = df
    with pd.ExcelWriter(args.xlsx, engine='openpyxl', mode='w') as w:
        for s, sdf in sheets.items():
            sdf.to_excel(w, sheet_name=s, index=False)
    print(f'Injected into {args.xlsx} as sheet "{args.sheet_name}"')
