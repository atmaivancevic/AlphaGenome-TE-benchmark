#!/usr/bin/env python3
"""
Fetch hg38 anchor + ALT sequence for one or more Schloissnig SVIM-asm variant
IDs, writing a tab file (ID, CHROM, POS, REF, ALT) for the AG variant-prediction
and scoring scripts. The literal ALT base string (anchor base + inserted seq) is only in
the BCF, so it's pulled here rather than from the supp tables.

Example usage (a few variants by ID):
python scripts/fig2_polymorphic_TE_example/fetch_variant_alt_seq.py \
    --ids SvimAsm00060017 SvimAsm00022857 SvimAsm00107233 SvimAsm00133580 \
    --out data/fig3_tier1_variants.tab

Example usage (bulk, IDs from a file):
python scripts/fig2_polymorphic_TE_example/fetch_variant_alt_seq.py \
    --ids $(cat /tmp/variant_ids_1322.txt | tr '\n' ' ') \
    --out data/fig3_all_testable_variants.tab
"""
import argparse, subprocess, sys, shutil

p = argparse.ArgumentParser()
p.add_argument('--ids', nargs='+', required=True, help='SvimAsm IDs to fetch')
p.add_argument('--bcf', default='data/schloissnig_2025/vcf/svim.asm.hg38.noGt.SVAN_1.3.bcf')
p.add_argument('--out', required=True, help='Output tab file (cols: ID CHROM POS REF ALT)')
args = p.parse_args()

if not shutil.which('bcftools'):
    sys.exit('bcftools not on PATH')

# bcftools -i 'ID=="X"||ID=="Y"' filter — single pass over the BCF.
expr = '||'.join(f'ID=="{i}"' for i in args.ids)
cmd = ['bcftools', 'query', '-i', expr, '-f', '%ID\t%CHROM\t%POS\t%REF\t%ALT\n', args.bcf]
out = subprocess.run(cmd, capture_output=True, text=True, check=True).stdout

found = {}
for line in out.strip().splitlines():
    vid, chrom, pos, ref, alt = line.split('\t')
    found[vid] = (chrom, pos, ref, alt)

with open(args.out, 'w') as fh:
    fh.write('ID\tCHROM\tPOS\tREF\tALT\n')
    for vid in args.ids:                       # preserve user-specified order
        if vid not in found:
            sys.stderr.write(f'WARN: {vid} not found in {args.bcf}\n'); continue
        chrom, pos, ref, alt = found[vid]
        fh.write(f'{vid}\t{chrom}\t{pos}\t{ref}\t{alt}\n')
        ins_len = len(alt) - len(ref)
        sys.stderr.write(f'{vid}\t{chrom}:{pos}\tINS {ins_len} bp\n')

sys.stderr.write(f'wrote {args.out}\n')
