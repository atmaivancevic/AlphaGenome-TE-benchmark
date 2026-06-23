#!/usr/bin/env python3
"""
Fetch hg38 anchor + ALT sequence for one or more Schloissnig SVIM-asm variant IDs.

Writes a tab-delimited file with the columns required by
`scripts/fig2_polymorphic_TE_example/plot_polymorphic_TE_insertion.py`: ID, CHROM, POS, REF, ALT.

Input  : Schloissnig hg38 SVIM-asm BCF (default: the no-GT SVAN_1.3 file —
         smaller, same coordinates + ALT seqs as the genotyped BCF).
Output : tab file (one row per variant) plus a stderr log line per ID with
         insert length so downstream issues are obvious.

Why a separate helper: the eQTL prototype's Supp Table 2 only stores
alt_length, not the actual inserted sequence — for AG variant calls we need
the literal ALT base string (anchor + insert), which lives only in the BCF.

Usage:
    python scripts/fig2_polymorphic_TE_example/fetch_variant_alt_seq.py \\
        --ids SvimAsm00060017 SvimAsm00042027 SvimAsm00107233 \\
        --out data/fig3_tier1_variants.tab
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
