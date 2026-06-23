"""
Annotate outlier candidates with peak width, nearest protein-coding gene +
distance, nearest-TSS distance, and a region class (promoter / proximal / distal).
Reads the TSV from find_panel_outliers.py, joins each peak against GENCODE v46
protein-coding genes, and writes an annotated TSV for use as a supp table. Also
prints a region_class x direction summary.

Region classes:
    promoter : |peak_mid - nearest TSS| < 2 kb
    proximal : 2 kb <= distance < 20 kb
    distal   : distance >= 20 kb

Example usage:
python scripts/fig1_baseline_chromatin/annotate_panel_outliers.py \
    --input  tasks/fig1_panel_BC_outlier_candidates.tsv \
    --gencode data/gencode.v46.annotation.feather \
    --output supptables/supp_table_fig1_BC_outliers.tsv
"""

import os, argparse, csv, bisect
from collections import defaultdict

import numpy as np
import pandas as pd

PROMOTER_MAX = 2_000     # bp
PROXIMAL_MAX = 20_000    # bp

DIRECTION_DISPLAY = {
    'B_over':  'AG over-predicts H3K27ac (red dots)',
    'C_under': 'AG under-predicts H3K27ac (blue dots)',
}


def load_pc_gene_index(gencode_path):
    """Build {chrom: (sorted_tss, names, starts, ends, mids)} for protein-coding genes."""
    g = pd.read_feather(gencode_path)
    g = g[(g['Feature'] == 'gene') & (g['gene_type'] == 'protein_coding')].copy()
    g['tss'] = np.where(g['Strand'] == '+', g['Start'], g['End'])
    g['mid'] = (g['Start'] + g['End']) // 2

    idx = {}
    for chrom, sub in g.groupby('Chromosome'):
        sub = sub.sort_values('tss').reset_index(drop=True)
        idx[chrom] = {
            'tss':    sub['tss'].to_numpy(),
            'starts': sub['Start'].to_numpy(),
            'ends':   sub['End'].to_numpy(),
            'mids':   sub['mid'].to_numpy(),
            'names':  sub['gene_name'].to_numpy(),
        }
    return idx


def nearest_tss(chrom, peak_mid, idx):
    """Return (gene_name, signed_dist_to_tss) for the closest TSS."""
    if chrom not in idx:
        return ('NA', None)
    tss_arr = idx[chrom]['tss']
    names = idx[chrom]['names']
    i = bisect.bisect_left(tss_arr, peak_mid)
    best = None  # (dist, name)
    for j in (i - 1, i):
        if 0 <= j < len(tss_arr):
            d = abs(int(tss_arr[j]) - peak_mid)
            if best is None or d < best[0]:
                best = (d, str(names[j]))
    return (best[1], best[0]) if best else ('NA', None)


def nearest_pc_gene(chrom, start, end, idx):
    """Return (gene_name, signed_dist) — 0 if peak overlaps the gene body."""
    if chrom not in idx:
        return ('NA', None)
    starts = idx[chrom]['starts']
    ends   = idx[chrom]['ends']
    mids   = idx[chrom]['mids']
    names  = idx[chrom]['names']

    overlap = np.where((starts <= end) & (ends >= start))[0]
    peak_mid = (start + end) // 2
    if len(overlap):
        i = overlap[np.argmin(np.abs(mids[overlap] - peak_mid))]
        return (str(names[i]), 0)

    i = bisect.bisect_left(list(mids), peak_mid)
    best = None
    for j in (i - 1, i):
        if 0 <= j < len(mids):
            d = abs(int(mids[j]) - peak_mid)
            if best is None or d < best[0]:
                best = (d, str(names[j]))
    return (best[1], best[0]) if best else ('NA', None)


def region_class(tss_dist):
    if tss_dist is None:
        return 'NA'
    if tss_dist < PROMOTER_MAX:
        return 'promoter'
    if tss_dist < PROXIMAL_MAX:
        return 'proximal'
    return 'distal'


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--input',  required=True, help='TSV from find_panel_outliers.py')
    ap.add_argument('--gencode', default='data/gencode.v46.annotation.feather')
    ap.add_argument('--output', required=True, help='Annotated TSV (supp table)')
    args = ap.parse_args()

    idx = load_pc_gene_index(args.gencode)

    out_rows = []
    with open(args.input) as f:
        for row in csv.DictReader(f, delimiter='\t'):
            chrom = row['chrom']
            s, e = int(row['start']), int(row['end'])
            width = e - s
            peak_mid = (s + e) // 2

            tss_gene, tss_dist = nearest_tss(chrom, peak_mid, idx)
            pc_gene,  pc_dist  = nearest_pc_gene(chrom, s, e, idx)
            cls = region_class(tss_dist)

            out_rows.append({
                'biosample': row['biosample'],
                'direction': DIRECTION_DISPLAY.get(row['direction'], row['direction']),
                'rank_in_biosample': row['rank'],
                'chrom': chrom,
                'start': s,
                'end': e,
                'width_bp': width,
                'exp_signal':  float(row['exp_signal']),
                'pred_signal': float(row['pred_signal']),
                'residual_pred_signal_minus_linear_fit': float(row['residual']),
                'nearest_gene': pc_gene,
                'distance_to_nearest_gene_bp': pc_dist if pc_dist is not None else '',
                'nearest_TSS': tss_gene,
                'distance_to_TSS_bp':  tss_dist if tss_dist is not None else '',
                'region_class': cls,
            })

    cols = ['biosample', 'direction', 'rank_in_biosample',
            'chrom', 'start', 'end', 'width_bp',
            'exp_signal', 'pred_signal',
            'residual_pred_signal_minus_linear_fit',
            'nearest_gene', 'distance_to_nearest_gene_bp',
            'nearest_TSS', 'distance_to_TSS_bp',
            'region_class']
    os.makedirs(os.path.dirname(os.path.abspath(args.output)), exist_ok=True)
    with open(args.output, 'w') as f:
        f.write('\t'.join(cols) + '\n')
        for r in out_rows:
            f.write('\t'.join(str(r[c]) for c in cols) + '\n')

    # ── Summary: region_class × direction (the headline observation) ──────
    counts = defaultdict(lambda: defaultdict(int))
    widths = defaultdict(list)
    for r in out_rows:
        counts[r['direction']][r['region_class']] += 1
        widths[r['direction']].append(r['width_bp'])

    print(f"\nWrote {len(out_rows)} annotated rows → {args.output}\n")
    print('Region-class breakdown (top 20 per biosample, 8 biosamples = 160 per direction):')
    short = {DIRECTION_DISPLAY['B_over']: 'over-pred (B)',
             DIRECTION_DISPLAY['C_under']: 'under-pred (C)'}
    print(f"  {'direction':14s}  {'promoter':>8s}  {'proximal':>8s}  {'distal':>8s}")
    for direction in DIRECTION_DISPLAY.values():
        c = counts[direction]
        total = sum(c.values())
        print(f"  {short[direction]:14s}  "
              f"{c['promoter']:>4d} ({100*c['promoter']/total:>4.1f}%)  "
              f"{c['proximal']:>4d} ({100*c['proximal']/total:>4.1f}%)  "
              f"{c['distal']:>4d} ({100*c['distal']/total:>4.1f}%)")
    print('\nPeak width (median / mean / max, bp):')
    for direction in DIRECTION_DISPLAY.values():
        w = np.array(widths[direction])
        print(f"  {short[direction]:14s}  median={int(np.median(w)):>7d}  "
              f"mean={int(w.mean()):>7d}  max={int(w.max()):>9d}")


if __name__ == '__main__':
    main()
