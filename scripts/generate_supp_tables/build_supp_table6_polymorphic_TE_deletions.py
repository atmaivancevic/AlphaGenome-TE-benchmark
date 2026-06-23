"""
Build Supp Table 6 (polymorphic TE deletion catalog): DEL-arm parallel to Supp
Table 3. Join the filtered solo TE DEL VCF with the genotyped BCF; per-variant
stats from BCF INFO, carrier lists + NA12878 GT + flags from per-sample GTs.
SVLEN is negative for DELs (anchor-excluded, VCF convention). NA12878=0/0 means
NA12878 retains the TE (matches hg38); carriers are those who lost it. Multi-allelic
sites skipped. Writes a TSV (xlsx optional via --xlsx).

Example usage:
python scripts/generate_supp_tables/build_supp_table6_polymorphic_TE_deletions.py \
    --filtered-vcf data/schloissnig_2025/vcf/solo_TE_DEL.svim_asm.hg38.SVAN.vcf.gz \
    --gt-bcf       data/schloissnig_2025/vcf/svim.asm.hg38.bcf \
    --tsv          supptables/supp_table_6_polymorphic_TE_deletions.tsv
"""
import argparse, gzip, re, subprocess
import pandas as pd

p = argparse.ArgumentParser()
p.add_argument('--filtered-vcf', required=True, help='solo TE DEL VCF (SVAN-annotated)')
p.add_argument('--gt-bcf', required=True, help='SVIM-asm hg38 genotyped BCF (908 samples)')
p.add_argument('--tsv', required=True, help='Output TSV path (authoritative; xlsx is composed manually at submission)')
p.add_argument('--xlsx', default=None, help='Optional: inject into SuppTables.xlsx (skipped by default; Atma maintains the workbook manually)')
p.add_argument('--sheet-name', default='Suppl Table 6 Polymorphic TE Deletions')
p.add_argument('--na-sample', default='NA12878', help='Flagship sample whose GT is reported as a column')
p.add_argument('--eqtl-testable-tsv', default='supptables/supp_table_7_top_gene_per_variant_DEL.tsv',
               help='DEL eQTL results table (Supp Table 7) whose variant_id set defines the '
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
        # SVAN's strand call for the deletion: '+' / '-' / '.' (undetermined).
        # Reported verbatim from SVAN; documented in Methods + supp legend.
        strand = re.search(r'(?:^|;)STRAND=([^;]+)', info)
        # The length column is `SVLEN`, sourced verbatim from the genotyped
        # BCF's INFO field (see step 2 below). Negative for DEL records
        # (e.g. -304 for a 304 bp deletion). Column name kept as `SVLEN` rather
        # than `del_len` so it's directly traceable back to the source VCF/BCF
        # without any renaming step. For 3,827 clean DEL records, SVLEN equals
        # len(ALT) - len(REF). For ~18 complex DEL records where len(ALT) > 1,
        # SVLEN reports the biologically-meaningful deletion size (typically
        # -(len(REF) - 1)) while len(ALT) - len(REF) would under-report it.
        # SVAN's DEL_LEN field is anchor-INCLUDED (1 bp larger in absolute
        # value than SVLEN); we deliberately do not use it.
        svan[vid] = dict(
            variant_id=vid, family=fam.group(1) if fam else None,
            chrom=cols[0], pos=int(cols[1]),
            # SVLEN filled below from genotyped BCF
            strand=strand.group(1) if strand else '.',
            alt=cols[4],
            dtype_n='solo',
            filter_svan=cols[6],
            perc_resolved=float(pct.group(1)) if pct else None,
            vcf_source=args.filtered_vcf,
        )
print(f'SVAN solo TE DEL: {len(svan):,}')

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
    if ',' in cols[4]:
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

for vid in multi_allelic_ids: svan.pop(vid, None)

print(f'Variants joined with genotyped BCF (biallelic): {matched:,}')
print(f'Multi-allelic skipped: {len(multi_allelic_ids):,}')

for vid, rec in svan.items():
    for k, _ in NUM_FIELDS: rec.setdefault(k, None)
    rec.setdefault(f'{args.na_sample}_GT', '')
    rec.setdefault(f'{args.na_sample}_alt_dose', None)
    rec.setdefault('n_carriers', None)
    rec.setdefault('carriers_hom', '')
    rec.setdefault('carriers_het', '')
    # SVLEN already in rec from NUM_FIELDS extraction above; no renaming.

# Flag columns linking the full catalog to the Fig 3 eQTL analysis (DEL arm):
#   na12878_00            — NA12878 retains the TE (matches hg38; the eQTL candidate pool
#                           that mitigates AlphaGenome's reference-allele bias).
#   eqtl_testable_MAGE260 — variant entered the MAGE-260 cis-eQTL test, i.e. it appears in
#                           Supp Table 7 (NA12878=0/0 + MAF>0.01 + >=5 hom + >=5 ref in MAGE-260).
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

cols = ['variant_id','family','chrom','pos','SVLEN','strand','alt',
        'dtype_n','filter_svan','perc_resolved',
        'AC','AN','AF','MAF','NS','F_MISSING',
        'AC_Hom','AC_Het','AC_Hemi','HWE','ExcHet',
        f'{args.na_sample}_GT', f'{args.na_sample}_alt_dose',
        'na12878_00','eqtl_testable_MAGE260',
        'n_carriers','carriers_hom','carriers_het','vcf_source']
family_rank = {'LTR5_Hs': 1, 'SVA': 2, 'L1': 3, 'Alu': 4}
df = pd.DataFrame(list(svan.values()))[cols]
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
