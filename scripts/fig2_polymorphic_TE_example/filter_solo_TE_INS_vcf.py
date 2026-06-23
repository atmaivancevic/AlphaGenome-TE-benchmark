"""
Filter the SVIM-asm hg38 SVAN-annotated BCF to specified SVAN classes on
canonical chromosomes (chr1..chr22, chrX, chrY). A record is kept if its
ITYPE_N (insertion class) or DTYPE_N (deletion class) matches one of the
requested values.

Values seen in the SVIM-asm hg38 SVAN BCF (release v1.1):
  ITYPE_N: solo, VNTR, DUP, DUP_INTERSPERSED, partnered, orphan, INV_DUP,
           PSD, NUMT, chimera, COMPLEX_DUP
  DTYPE_N: VNTR, solo, orphan, PSD, partnered, chimera

Default keeps only `ITYPE_N=solo` (canonical solo MEI insertions — cleanest
substrate for sequence-to-function benchmarks). Pass --dtype to also pull
deletion classes; pass an empty --itype to exclude insertions entirely.

Output is a gzipped VCF. The SVAN BCF ships with stale T2T contig metadata
in its header (variant positions are hg38). The output strips those and
substitutes the hg38 contig lines from the genotyped BCF, so downstream
tools see consistent GRCh38 metadata.

Usage:
    # Pass 1: solo MEI insertions only
    python scripts/fig2_polymorphic_TE_example/filter_solo_TE_INS_vcf.py \\
        --svan-bcf data/schloissnig_2025/vcf/svim.asm.hg38.noGt.SVAN_1.3.bcf \\
        --gt-bcf   data/schloissnig_2025/vcf/svim.asm.hg38.bcf \\
        --output   data/schloissnig_2025/vcf/solo_TE_INS.svim_asm.hg38.SVAN.vcf.gz

    # Everything except solo MEI insertions (other INS classes + all DEL classes)
    python scripts/fig2_polymorphic_TE_example/filter_solo_TE_INS_vcf.py \\
        --svan-bcf data/schloissnig_2025/vcf/svim.asm.hg38.noGt.SVAN_1.3.bcf \\
        --gt-bcf   data/schloissnig_2025/vcf/svim.asm.hg38.bcf \\
        --itype    VNTR,DUP,DUP_INTERSPERSED,partnered,orphan,INV_DUP,PSD,NUMT,chimera,COMPLEX_DUP \\
        --dtype    VNTR,solo,orphan,PSD,partnered,chimera \\
        --output   data/schloissnig_2025/vcf/nonsolo_INS_DEL.svim_asm.hg38.SVAN.vcf.gz
"""
import argparse, gzip, hashlib, re, subprocess
from collections import Counter

CANONICAL = {f'chr{i}' for i in range(1, 23)} | {'chrX', 'chrY'}

p = argparse.ArgumentParser()
p.add_argument('--svan-bcf', required=True, help='SVIM-asm hg38 SVAN-annotated BCF (no GTs)')
p.add_argument('--gt-bcf', required=True, help='SVIM-asm hg38 genotyped BCF (used for hg38 contig header)')
p.add_argument('--itype', default='solo', help='Comma-separated ITYPE_N values to keep (default: solo)')
p.add_argument('--dtype', default='', help='Comma-separated DTYPE_N values to keep (default: none)')
p.add_argument('--output', required=True, help='Output gzipped VCF')
args = p.parse_args()

itypes = set(t for t in args.itype.split(',') if t)
dtypes = set(t for t in args.dtype.split(',') if t)

# Pull canonical hg38 ##contig lines from the genotyped BCF
gt_header = subprocess.run(['bcftools', 'view', '-h', args.gt_bcf],
                           capture_output=True, text=True, check=True).stdout
hg38_contigs = []
for line in gt_header.splitlines():
    if line.startswith('##contig'):
        m = re.search(r'ID=([^,>]+)', line)
        if m and m.group(1) in CANONICAL:
            hg38_contigs.append(line)

header_other = []
n_total = 0
kept_by_type = Counter()
proc = subprocess.Popen(['bcftools', 'view', args.svan_bcf], stdout=subprocess.PIPE, text=True)
with gzip.open(args.output, 'wt') as out:
    for line in proc.stdout:
        if line.startswith('##contig'):
            continue
        if line.startswith('##'):
            header_other.append(line); continue
        if line.startswith('#CHROM'):
            for h in header_other: out.write(h)
            for c in hg38_contigs: out.write(c + '\n')
            out.write(line); continue
        n_total += 1
        cols = line.split('\t', 8)
        if cols[0] not in CANONICAL: continue
        info = cols[7] + ';'  # sentinel so ?TYPE_N=foo at end of INFO matches with ;
        mi = re.search(r'ITYPE_N=([^;]+);', info)
        md = re.search(r'DTYPE_N=([^;]+);', info)
        keep_label = None
        if mi and mi.group(1) in itypes:
            keep_label = ('INS', mi.group(1))
        elif md and md.group(1) in dtypes:
            keep_label = ('DEL', md.group(1))
        if keep_label is None: continue
        out.write(line)
        kept_by_type[keep_label] += 1
proc.wait()
if proc.returncode != 0:
    raise SystemExit(f'bcftools view failed (rc={proc.returncode})')

with open(args.output, 'rb') as f: md5 = hashlib.md5(f.read()).hexdigest()

n_kept = sum(kept_by_type.values())
print(f'Input records:      {n_total:,}')
print(f'Kept (canonical chroms + ITYPE_N/DTYPE_N match): {n_kept:,}')
if itypes:
    print(f'Insertion (ITYPE_N) breakdown:')
    for t in sorted(itypes):
        print(f'    {t:<20s} {kept_by_type[("INS", t)]:>7,}')
if dtypes:
    print(f'Deletion (DTYPE_N) breakdown:')
    for t in sorted(dtypes):
        print(f'    {t:<20s} {kept_by_type[("DEL", t)]:>7,}')
print(f'Output: {args.output}')
print(f'md5:    {md5}')
