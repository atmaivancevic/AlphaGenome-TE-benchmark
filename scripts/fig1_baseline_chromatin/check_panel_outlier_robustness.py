"""
Sanity-check the Fig 1 Panel B/C claim that AG over-calls cluster at
promoters and under-calls at distal enhancers.

Runs four checks:
  1. Baseline region-class distribution of ALL matched peaks across the
     8 biosamples (not just outliers). Required to interpret the B/C
     percentages: comparing B/C to each other tells us they differ;
     comparing each to baseline tells us *which way* AG is biased.
  2. Binomial test of each direction's promoter/distal fraction vs the
     baseline proportion. Two-sided.
  3. Per-biosample stability of B promoter% and C distal%.
  4. Top-N sensitivity (re-rank with topn ∈ {5, 20, 50}). Checks that
     the headline isn't a top-20 artifact.

Usage:
    python scripts/fig1_baseline_chromatin/check_panel_outlier_robustness.py \
        --supp-table supptables/supp_table_fig1_BC_outliers.tsv
"""

import os, argparse, csv, gzip, bisect
from collections import defaultdict
import numpy as np
import pandas as pd
from scipy import stats

PROMOTER_MAX = 2_000     # bp
PROXIMAL_MAX = 20_000    # bp

# Direction labels — must match the human-readable values written by
# annotate_panel_outliers.py into the supp table.
DIR_B = 'AG over-predicts H3K27ac (red dots)'
DIR_C = 'AG under-predicts H3K27ac (blue dots)'

BIOSAMPLES = [
    ('HCT116',   'EFO_0002824_HCT116'),
    ('GM12878',  'EFO_0002784_GM12878'),
    ('K562',     'EFO_0002067_K562'),
    ('OCI-LY3',  'EFO_0006710_OCI-LY3'),
    ('H1',       'EFO_0003042_H1'),
    ('HUES6',    'EFO_0007086_HUES6'),
    ('placenta', 'UBERON_0001987_placenta'),
    ('PC-9',     'EFO_0002847_PC-9'),
]

SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
PROJECT_ROOT = os.path.abspath(os.path.join(SCRIPT_DIR, '..', '..'))
EXP_DIR  = os.path.join(PROJECT_ROOT, 'data/encode_h3k27ac/peaks/merged')
PRED_DIR = os.path.join(PROJECT_ROOT, 'results/AG_predicted_h3k27ac_batched')


def build_tss_index(gencode_path):
    g = pd.read_feather(gencode_path)
    g = g[(g['Feature']=='gene') & (g['gene_type']=='protein_coding')].copy()
    g['tss'] = np.where(g['Strand']=='+', g['Start'], g['End'])
    return {chrom: np.sort(sub['tss'].values) for chrom, sub in g.groupby('Chromosome')}


def region_class(chrom, start, end, tss_idx):
    if chrom not in tss_idx: return 'NA'
    arr = tss_idx[chrom]
    peak_mid = (start+end)//2
    i = bisect.bisect_left(arr, peak_mid)
    best = None
    for j in (i-1, i):
        if 0 <= j < len(arr):
            d = abs(int(arr[j]) - peak_mid)
            if best is None or d < best: best = d
    if best is None: return 'NA'
    if best < PROMOTER_MAX: return 'promoter'
    if best < PROXIMAL_MAX: return 'proximal'
    return 'distal'


def load_pairs(slug):
    out = {}
    for kind, path in (('exp',  os.path.join(EXP_DIR,  f'{slug}.narrowPeak.gz')),
                       ('pred', os.path.join(PRED_DIR, slug, f'{slug}.narrowPeak.gz'))):
        d = {}
        with gzip.open(path,'rt') as f:
            for line in f:
                if line.startswith('#') or not line.strip(): continue
                c = line.rstrip('\n').split('\t')
                d[(c[0],int(c[1]),int(c[2]))] = float(c[6])
        out[kind] = d
    return out['exp'], out['pred']


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--supp-table', required=True,
                    help='Annotated supp table from annotate_panel_outliers.py')
    ap.add_argument('--gencode', default='data/gencode.v46.annotation.feather')
    args = ap.parse_args()

    tss_idx = build_tss_index(args.gencode)

    # ── 1. Baseline ──────────────────────────────────────────────────────
    print("="*70)
    print("CHECK 1: Region-class distribution of ALL matched peaks (baseline)")
    print("="*70)
    baseline = defaultdict(int)
    per_bs_baseline = {}
    for label, slug in BIOSAMPLES:
        e, p = load_pairs(slug)
        shared = set(e) & set(p)
        bs = defaultdict(int)
        for c in shared:
            cls = region_class(*c, tss_idx)
            bs[cls] += 1
            baseline[cls] += 1
        per_bs_baseline[label] = bs
    total = sum(baseline.values())
    print(f"\nAll matched peaks across 8 biosamples (n = {total:,}):")
    for cls in ('promoter','proximal','distal','NA'):
        n = baseline[cls]
        print(f"  {cls:9s} {n:>7,} ({100*n/total:>5.1f}%)")
    p_prom  = baseline['promoter'] / (baseline['promoter']+baseline['proximal']+baseline['distal'])
    p_dist  = baseline['distal']   / (baseline['promoter']+baseline['proximal']+baseline['distal'])
    print(f"\nBaseline: P(promoter) = {p_prom:.3f}, P(distal) = {p_dist:.3f}")

    # ── 2. Binomial vs baseline ─────────────────────────────────────────
    print("\n" + "="*70)
    print("CHECK 2: B and C vs baseline (two-sided binomial)")
    print("="*70)
    rows = list(csv.DictReader(open(args.supp_table), delimiter='\t'))
    for direction in (DIR_B,DIR_C):
        c = defaultdict(int)
        for r in rows:
            if r['direction']==direction: c[r['region_class']] += 1
        n = c['promoter']+c['proximal']+c['distal']
        print(f"\n{direction} (n={n}):")
        for cls, base_p in (('promoter', p_prom), ('distal', p_dist)):
            obs = c[cls]
            exp = n * base_p
            fold = (obs/n) / base_p
            pval = stats.binomtest(obs, n, base_p, alternative='two-sided').pvalue
            arrow = '↑' if fold > 1 else '↓'
            print(f"  {cls:8s}: observed {obs}, expected {exp:.1f}  "
                  f"({100*obs/n:.1f}% vs baseline {100*base_p:.1f}%)  "
                  f"{arrow} {fold:.2f}×  p = {pval:.2e}")

    # ── 3. Per-biosample stability ──────────────────────────────────────
    print("\n" + "="*70)
    print("CHECK 3: Per-biosample stability (top 20)")
    print("="*70)
    print(f"{'biosample':9s}  {'B prom%':>7s} {'B dist%':>7s}  "
          f"{'C prom%':>7s} {'C dist%':>7s}  {'base prom%':>10s}")
    for label, _ in BIOSAMPLES:
        b = defaultdict(int); c = defaultdict(int)
        for r in rows:
            if r['biosample']!=label: continue
            (b if r['direction']==DIR_B else c)[r['region_class']] += 1
        bn = b['promoter']+b['proximal']+b['distal']
        cn = c['promoter']+c['proximal']+c['distal']
        bb = per_bs_baseline[label]
        bbn = bb['promoter']+bb['proximal']+bb['distal']
        print(f"{label:9s}  {100*b['promoter']/bn:>6.1f}% {100*b['distal']/bn:>6.1f}%  "
              f"{100*c['promoter']/cn:>6.1f}% {100*c['distal']/cn:>6.1f}%  "
              f"{100*bb['promoter']/bbn:>9.1f}%")

    # ── 4. Top-N sensitivity ────────────────────────────────────────────
    print("\n" + "="*70)
    print("CHECK 4: Top-N sensitivity (re-rank, region_class proportions)")
    print("="*70)
    for topn in (5, 20, 50):
        bc = defaultdict(lambda: defaultdict(int))
        for label, slug in BIOSAMPLES:
            e, p = load_pairs(slug)
            shared = sorted(set(e) & set(p))
            coords = np.array(shared, dtype=object)
            ev = np.array([e[c] for c in shared])
            pv = np.array([p[c] for c in shared])
            slope, intercept, *_ = stats.linregress(ev, pv)
            resid = pv - (slope*ev + intercept)
            order = np.argsort(resid)
            for direction, idxs in ((DIR_B, order[-topn:]), (DIR_C, order[:topn])):
                for i in idxs:
                    bc[direction][region_class(*coords[i], tss_idx)] += 1
        print(f"\ntop-{topn} (n per direction = {topn*8}):")
        for direction in (DIR_B,DIR_C):
            c = bc[direction]
            n = c['promoter']+c['proximal']+c['distal']
            print(f"  {direction:8s}  prom={100*c['promoter']/n:>4.1f}%  "
                  f"prox={100*c['proximal']/n:>4.1f}%  "
                  f"dist={100*c['distal']/n:>4.1f}%")


if __name__ == '__main__':
    main()
