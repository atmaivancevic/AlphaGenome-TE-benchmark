"""
Build Supp Table 2 (per-donor cohort flags) for the Schloissnig 908 phased
individuals. One row per donor: population/sex, Tier 1 RNA-seq flags (MAGE /
GEUVADIS), Tier 2 chromatin flags (ATAC / ChIP / DNase across 7 sources),
aggregated tier counts, and per-source accessions. Inputs are listed in the
README provenance block. Output: supptables/supp_table_2_individuals.tsv.

GM <-> NA: 1KGP IDs use HG/NA; Coriell uses GM for HapMap YRI/CEU lines
(GM18498 = NA18498). MAGE gives the explicit map; otherwise NA->GM substring sub.

Example usage:
python scripts/generate_supp_tables/build_supp_table2_individuals.py \
    --tsv supptables/supp_table_2_individuals.tsv
"""
import argparse, csv, gzip, re, sys
from pathlib import Path
import pandas as pd

p = argparse.ArgumentParser()
p.add_argument('--phased',   default='data/schloissnig_2025/phased_908_samples.txt')
p.add_argument('--popfile',  default='data/schloissnig_2025/bundle_extracted/sample.tsv')
p.add_argument('--mage',     default='data/rna_cohorts/mage_metadata.tsv')
p.add_argument('--mage-ena', default='data/rna_cohorts/mage_ena_fastq.tsv')
p.add_argument('--mage-md5', default='scripts/baffin_downloads/manifests_md5/MAGE.tsv',
               help='Baffin md5 manifest: sample\\tfilename\\tmd5')
p.add_argument('--geuvadis', default='data/rna_cohorts/geuvadis_samples.txt')
p.add_argument('--geuv-ena', default='data/rna_cohorts/geuvadis_ena_fastq.tsv')
p.add_argument('--geuv-md5', default='scripts/baffin_downloads/manifests_md5/GEUVADIS.tsv',
               help='Baffin md5 manifest: sample\\tfilename\\tmd5')
p.add_argument('--kumasaka', default='data/lcl_chromatin/kumasaka2019_ENA_PRJEB28318_analysis.tsv')
p.add_argument('--waszak',   default='data/lcl_chromatin/waszak2015_E-MTAB-3657.sdrf.txt')
p.add_argument('--grubert',  default='data/lcl_chromatin/grubert2015_GSE58852_series_matrix.txt.gz')
p.add_argument('--degner',   default='data/lcl_chromatin/degner2012_GSE31388_series_matrix.txt.gz')
p.add_argument('--afgr',     default='data/lcl_chromatin/afgr_atac_ENCODE.tsv')
p.add_argument('--encode',   default='data/lcl_chromatin/encode_NA12878_NA19238.tsv')
p.add_argument('--xlsx',     default='supptables/SuppTables.xlsx')
p.add_argument('--tsv',      default='supptables/supp_table_2_individuals.tsv')
p.add_argument('--sheet-name', default='Suppl Table 2 Individuals')
args = p.parse_args()

# 908 master roster — NA12878 (GM12878) pinned to row 1 as the flagship LCL,
# everything else in alphabetical order.
phased = sorted(Path(args.phased).read_text().split())
if 'NA12878' in phased:
    phased = ['NA12878'] + [s for s in phased if s != 'NA12878']
phased_set = set(phased)
print(f'Schloissnig phased: {len(phased)} (first row: {phased[0]})')

# Pop / sex
pop = {}
with open(args.popfile) as f:
    r = csv.DictReader(f, delimiter='\t')
    for row in r: pop[row['sample']] = (row['pop'], row['sex'])
SUPERPOP = {  # 1KGP super-pop
    'CEU':'EUR','TSI':'EUR','FIN':'EUR','GBR':'EUR','IBS':'EUR',
    'YRI':'AFR','LWK':'AFR','MSL':'AFR','ESN':'AFR','GWD':'AFR','ASW':'AFR','ACB':'AFR',
    'CHB':'EAS','JPT':'EAS','CHS':'EAS','CDX':'EAS','KHV':'EAS',
    'GIH':'SAS','PJL':'SAS','BEB':'SAS','STU':'SAS','ITU':'SAS',
    'MXL':'AMR','PUR':'AMR','CLM':'AMR','PEL':'AMR'}

# MAGE: explicit kgp_id ↔ coriell_id map + RNA tier 1
mage_kgp_to_coriell, mage_set = {}, set()
with open(args.mage) as f:
    r = csv.DictReader(f, delimiter='\t')
    for row in r:
        mage_kgp_to_coriell.setdefault(row['kgp_id'], row['coriell_id'])
        mage_set.add(row['kgp_id'])

geuv_set = set(Path(args.geuvadis).read_text().split())

# Per-donor FASTQ URLs (RNA tier 1)
def load_ena_fastq(path, id_extract):
    """Return {sample_id: [(run_accession, fastq_url1, fastq_url2_or_None), ...]}."""
    out = {}
    with open(path) as f:
        r = csv.DictReader(f, delimiter='\t')
        for row in r:
            sid = id_extract(row)
            if not sid: continue
            urls = row.get('fastq_ftp','').split(';')
            urls = ['https://' + u for u in urls if u]
            out.setdefault(sid, []).append((row['run_accession'], urls))
    return out

# MAGE: sample_alias is the kgp_id directly (the ENA browser shows the full
# "HG00100_batch04_rep1" form; the filereport export here uses the bare HG/NA)
mage_runs = load_ena_fastq(args.mage_ena, lambda row: row.get('sample_alias','').strip())
# GEUVADIS: sample_alias is "GEUV:NA12878"; sample_title is bare HG/NA
geuv_runs = load_ena_fastq(args.geuv_ena, lambda row: row.get('sample_title','').strip())

# Baffin md5 manifests: (sample, filename) → md5. Used to attach verified hashes
# to each URL in Supp Table 2 so the table is self-auditable.
def load_md5(path):
    out = {}
    with open(path) as f:
        r = csv.DictReader(f, delimiter='\t')
        for row in r:
            out[(row['sample'], row['filename'])] = row['md5']
    return out
mage_md5 = load_md5(args.mage_md5)
geuv_md5 = load_md5(args.geuv_md5)
print(f'md5 manifests: MAGE {len(mage_md5)} files, GEUVADIS {len(geuv_md5)} files')

def md5s_per_run(sid, runs, md5_map, missing, used):
    """Return ';'-joined md5 groups, one group per run, ','-joined within a run.
    e.g. for two PE runs: 'm1,m2;m3,m4'. Aligns 1:1 with the per-run URL list."""
    out = []
    for _, ftp_urls in runs:
        run_md5s = []
        for u in ftp_urls:
            fn = u.rsplit('/', 1)[-1]
            used.add((sid, fn))
            h = md5_map.get((sid, fn), '')
            if not h: missing.append((sid, fn))
            run_md5s.append(h)
        out.append(','.join(run_md5s))
    return ';'.join(out)

mage_md5_missing, geuv_md5_missing = [], []
mage_used, geuv_used = set(), set()  # (sample, filename) pairs cited in supp table

def coriell_for(sid):
    """Coriell catalog ID. MAGE-explicit if present; else NA→GM substitution."""
    if sid in mage_kgp_to_coriell: return mage_kgp_to_coriell[sid]
    return sid.replace('NA','GM') if sid.startswith('NA') else sid

# Kumasaka: 1KGP ID is in analysis_title column
kumasaka_acc, kumasaka_url = {}, {}  # sid → ENA analysis accession / CRAM URL
with open(args.kumasaka) as f:
    r = csv.DictReader(f, delimiter='\t')
    for row in r:
        sid = row['analysis_title'].strip()
        if re.match(r'^(HG|NA)\d+$', sid):
            kumasaka_acc.setdefault(sid, row['analysis_accession'])
            # submitted_ftp is "ftp.../X.cram;ftp.../X.cram.crai"
            urls = ['https://' + u for u in row.get('submitted_ftp','').split(';') if u]
            kumasaka_url.setdefault(sid, ';'.join(urls))

# Waszak: per-donor assays (need full H3K27ac+H3K4me1+H3K4me3 trio to count)
# + per-donor FASTQ URLs aggregated across all antibodies / runs
waszak_assays, waszak_fastq = {}, {}
with open(args.waszak) as f:
    r = csv.DictReader(f, delimiter='\t')
    for row in r:
        cor = row.get('Characteristics[coriell id]','').strip()
        ip  = row.get('Factor Value[immunoprecipitate]','').strip()
        uri = row.get('Comment[FASTQ_URI]','').strip()
        if cor and ip:
            waszak_assays.setdefault(cor, set()).add(ip)
        if cor and uri:
            waszak_fastq.setdefault(cor, []).append(uri.replace('ftp://','https://'))
TRIO = {'H3K27ac','H3K4me1','H3K4me3'}
waszak_trio = {s for s,a in waszak_assays.items() if TRIO.issubset(a)}

# Grubert (H3K27ac arm — GSE58852)
grubert_gms = set()
with gzip.open(args.grubert,'rt') as f:
    for line in f:
        if line.startswith('!Sample_title'):
            for t in re.findall(r'"([^"]+)"', line):
                m = re.match(r'(GM\d+)', t)
                if m: grubert_gms.add(m.group(1))
            break

# Degner DNase-seq (GSE31388)
degner_gms = set()
with gzip.open(args.degner,'rt') as f:
    for line in f:
        if line.startswith('!Sample_title'):
            for t in re.findall(r'"([^"]+)"', line):
                m = re.search(r'(GM\d+)', t)
                if m: degner_gms.add(m.group(1))
            break

# AFGR ATAC (ENCODE)
afgr_acc = {}
with open(args.afgr) as f:
    r = csv.DictReader(f, delimiter='\t')
    for row in r:
        afgr_acc.setdefault(row['sample'], row['encode_accession'])

# ENCODE NA12878 + NA19238 (per-assay counts + accession lists)
encode_counts, encode_accs = {}, {}  # (sample, assay) → count / accession list
with open(args.encode) as f:
    r = csv.DictReader(f, delimiter='\t')
    for row in r:
        k = (row['sample'], row['assay'])
        encode_counts[k] = encode_counts.get(k, 0) + 1
        encode_accs.setdefault(k, []).append(row['encode_accession'])

def encode_search_url(sample, assay):
    accs = encode_accs.get((sample, assay), [])
    return ';'.join(f'https://www.encodeproject.org/experiments/{a}/' for a in accs)

# Build per-donor rows
GEO_ACC = 'https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc='
rows = []
for sid in phased:
    pop_, sex = pop.get(sid, ('NA','NA'))
    cor = coriell_for(sid)

    rna_mage = sid in mage_set
    rna_geuv = sid in geuv_set
    in_tier1 = rna_mage or rna_geuv

    # One browser URL per run (no duplication for PE files); md5s grouped per
    # run with ',' within a run and ';' between runs so each URL ↔ md5-group
    # aligns 1:1.
    mage_runs_d  = mage_runs.get(sid, [])
    geuv_runs_d  = geuv_runs.get(sid, [])
    mage_run_ids   = ';'.join(r[0] for r in mage_runs_d)
    geuv_run_ids   = ';'.join(r[0] for r in geuv_runs_d)
    mage_fastq_url = ';'.join(f'https://www.ebi.ac.uk/ena/browser/view/{r[0]}' for r in mage_runs_d)
    geuv_fastq_url = ';'.join(f'https://www.ebi.ac.uk/ena/browser/view/{r[0]}' for r in geuv_runs_d)
    mage_fastq_md5 = md5s_per_run(sid, mage_runs_d, mage_md5, mage_md5_missing, mage_used)
    geuv_fastq_md5 = md5s_per_run(sid, geuv_runs_d, geuv_md5, geuv_md5_missing, geuv_used)

    # Chromatin per-source flags
    atac_NA12878  = encode_counts.get((sid, 'ATAC-seq'), 0) > 0
    atac_kumasaka = sid in kumasaka_acc
    atac_afgr     = sid in afgr_acc
    chip_NA12878  = (encode_counts.get((sid,'Histone ChIP-seq'),0)
                     + encode_counts.get((sid,'TF ChIP-seq'),0) > 0)
    chip_grubert  = cor in grubert_gms
    chip_waszak   = cor in waszak_trio
    dnase_NA12878 = encode_counts.get((sid,'DNase-seq'), 0) > 0 and sid == 'NA12878'
    dnase_NA19238 = encode_counts.get((sid,'DNase-seq'), 0) > 0 and sid == 'NA19238'
    dnase_degner  = cor in degner_gms

    has_atac  = atac_NA12878 or atac_kumasaka or atac_afgr
    has_chip  = chip_NA12878 or chip_grubert or chip_waszak
    has_dnase = dnase_NA12878 or dnase_NA19238 or dnase_degner
    n_chrom = int(has_atac) + int(has_chip) + int(has_dnase)
    in_tier2 = n_chrom > 0

    rows.append(dict(
        sample=sid, coriell_id=cor, pop=pop_, superpop=SUPERPOP.get(pop_,''), sex=sex,
        tier0_carrier=True,
        # RNA tier 1
        tier1_rna_mage=rna_mage, tier1_rna_geuvadis=rna_geuv, in_tier1=in_tier1,
        mage_run_accessions=mage_run_ids,
        mage_fastq_urls=mage_fastq_url,
        mage_fastq_md5s=mage_fastq_md5,
        geuvadis_run_accessions=geuv_run_ids,
        geuvadis_fastq_urls=geuv_fastq_url,
        geuvadis_fastq_md5s=geuv_fastq_md5,
        # Chromatin per-source flags
        tier2_atac_NA12878_ENCODE=atac_NA12878,
        tier2_atac_AFGR=atac_afgr,
        tier2_atac_Kumasaka=atac_kumasaka,
        tier2_chip_NA12878_ENCODE=chip_NA12878,
        tier2_chip_Grubert=chip_grubert,
        tier2_chip_Waszak_trio=chip_waszak,
        tier2_dnase_NA12878_ENCODE=dnase_NA12878,
        tier2_dnase_NA19238_ENCODE=dnase_NA19238,
        tier2_dnase_Degner=dnase_degner,
        has_atac=has_atac, has_chip=has_chip, has_dnase=has_dnase,
        n_chrom_assays=n_chrom, in_tier2=in_tier2,
        # Per-source download URLs / accessions
        kumasaka_ENA_analysis=kumasaka_acc.get(sid,''),
        kumasaka_cram_urls=kumasaka_url.get(sid,''),
        waszak_assays=','.join(sorted(waszak_assays.get(cor,set()))),
        waszak_fastq_urls=';'.join(waszak_fastq.get(cor,[])),
        grubert_GEO_url=(f'{GEO_ACC}GSE58852' if chip_grubert else ''),
        degner_GEO_url=(f'{GEO_ACC}GSE31388' if dnase_degner else ''),
        afgr_ENCODE_experiment=afgr_acc.get(sid,''),
        afgr_ENCODE_url=(f'https://www.encodeproject.org/experiments/{afgr_acc[sid]}/' if sid in afgr_acc else ''),
        encode_NA12878_ATAC_urls=encode_search_url(sid,'ATAC-seq')   if sid=='NA12878' else '',
        encode_NA12878_DNase_urls=encode_search_url(sid,'DNase-seq') if sid=='NA12878' else '',
        encode_NA12878_TF_ChIP_n=encode_counts.get((sid,'TF ChIP-seq'),0),
        encode_NA12878_histone_ChIP_n=encode_counts.get((sid,'Histone ChIP-seq'),0),
        encode_NA12878_TF_ChIP_search=(f'https://www.encodeproject.org/search/?type=Experiment&assay_title=TF+ChIP-seq&biosample_ontology.term_name=GM12878&status=released' if sid=='NA12878' else ''),
        encode_NA12878_histone_ChIP_search=(f'https://www.encodeproject.org/search/?type=Experiment&assay_title=Histone+ChIP-seq&biosample_ontology.term_name=GM12878&status=released' if sid=='NA12878' else ''),
        encode_NA19238_DNase_urls=encode_search_url(sid,'DNase-seq') if sid=='NA19238' else '',
    ))

df = pd.DataFrame(rows)
print(f'rows: {len(df)} / cols: {len(df.columns)}')
print('Tier 1:', int(df.in_tier1.sum()),
      ' Tier 2:', int(df.in_tier2.sum()),
      ' all-3-assays:', int((df.n_chrom_assays==3).sum()))
print('Per-source counts:')
for col in df.columns:
    if col.startswith(('tier1_','tier2_','has_')):
        print(f'  {col:38s} {int(df[col].sum())}')

# Cross-checks: Supp Table URLs vs Baffin md5 manifests must be 1:1
def cross_check(label, used, md5_map, missing):
    orphans = sorted(set(md5_map) - used)
    print(f'{label}: {len(used)} (sample,filename) cited / {len(md5_map)} on Baffin; '
          f'missing-md5: {len(missing)}, manifest-orphans: {len(orphans)}')
    if missing: print(f'  missing[:5]: {missing[:5]}')
    if orphans: print(f'  orphans[:5]: {orphans[:5]}')
    return len(missing) == 0 and len(orphans) == 0

ok_m = cross_check('MAGE',     mage_used, mage_md5, mage_md5_missing)
ok_g = cross_check('GEUVADIS', geuv_used, geuv_md5, geuv_md5_missing)
if not (ok_m and ok_g):
    sys.exit('Supp Table 2 ↔ Baffin manifest mismatch — fix before publishing')

Path(args.tsv).parent.mkdir(parents=True, exist_ok=True)
df.to_csv(args.tsv, sep='\t', index=False)
print(f'wrote {args.tsv}')

# Inject into xlsx, preserving existing sheets
xlsx = Path(args.xlsx)
if xlsx.exists():
    existing = pd.ExcelFile(xlsx)
    sheets = {s: pd.read_excel(xlsx, sheet_name=s) for s in existing.sheet_names}
else:
    sheets = {}
sheets[args.sheet_name] = df
with pd.ExcelWriter(xlsx, engine='openpyxl', mode='w') as w:
    for s, sdf in sheets.items():
        sdf.to_excel(w, sheet_name=s, index=False)
print(f'injected sheet "{args.sheet_name}" into {xlsx}')
