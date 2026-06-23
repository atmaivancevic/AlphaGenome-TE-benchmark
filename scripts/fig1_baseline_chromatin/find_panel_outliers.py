"""
Find Fig 1 Panel B/C outlier candidates from the 8 AG-vs-ENCODE H3K27ac
scatter pairs.

For each biosample, peaks are matched by coordinate (same join used in
`plot_peak_correlation.py`), a linear regression is fit to (exp, pred)
on linear (raw signal) axes — matching the Panel A scatter — and peaks
are ranked by signed residual:

    residual = pred − (slope * exp + intercept)

  * Big positive residual → AG over-predicts (Panel B candidate).
  * Big negative residual → AG under-predicts (Panel C candidate).

Linear residuals naturally weight by signal magnitude, so the dots that
look most extreme on the scatter are the ones flagged. No floor filter
needed: low-signal noise has small residuals by construction.

Output: one TSV pooling the top-N candidates per direction per biosample.

Usage:
    python scripts/fig1_baseline_chromatin/find_panel_outliers.py \
        --topn 20 \
        --output tasks/fig1_panel_BC_outlier_candidates.tsv
"""

import os, argparse, gzip
import numpy as np
from scipy import stats

SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
PROJECT_ROOT = os.path.abspath(os.path.join(SCRIPT_DIR, '..', '..'))

# Panel A biosamples — keep in the same order as the figure grid.
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

EXP_DIR  = os.path.join(PROJECT_ROOT, 'data/encode_h3k27ac/peaks/merged')
PRED_DIR = os.path.join(PROJECT_ROOT, 'results/AG_predicted_h3k27ac_batched')


def load_narrowpeak(path):
    opener = gzip.open if path.endswith('.gz') else open
    coords, signals = [], []
    with opener(path, 'rt') as f:
        for line in f:
            if line.startswith('#') or not line.strip():
                continue
            cols = line.rstrip('\n').split('\t')
            coords.append((cols[0], int(cols[1]), int(cols[2])))
            signals.append(float(cols[6]))
    return coords, np.array(signals)


def biosample_outliers(label, slug, topn):
    exp_path  = os.path.join(EXP_DIR,  f'{slug}.narrowPeak.gz')
    pred_path = os.path.join(PRED_DIR, slug, f'{slug}.narrowPeak.gz')
    exp_coords, exp_sig   = load_narrowpeak(exp_path)
    pred_coords, pred_sig = load_narrowpeak(pred_path)

    exp_dict  = dict(zip(exp_coords, exp_sig))
    pred_dict = dict(zip(pred_coords, pred_sig))
    shared = sorted(set(exp_dict) & set(pred_dict))

    coords = np.array(shared, dtype=object)
    e = np.array([exp_dict[c]  for c in shared])
    p = np.array([pred_dict[c] for c in shared])

    # Linear regression on raw signal — matches the Panel A scatter axes.
    slope, intercept, r, _, _ = stats.linregress(e, p)
    resid = p - (slope * e + intercept)

    order  = np.argsort(resid)
    bottom = order[:topn]            # Panel C: under-predicted (resid << 0)
    top    = order[-topn:][::-1]     # Panel B: over-predicted  (resid >> 0)

    rows = []
    for direction, idxs in (('B_over', top), ('C_under', bottom)):
        for rank, i in enumerate(idxs, 1):
            chrom, start, end = coords[i]
            rows.append({
                'biosample': label,
                'biosample_slug': slug,
                'direction': direction,
                'rank': rank,
                'chrom': chrom,
                'start': start,
                'end': end,
                'exp_signal': float(e[i]),
                'pred_signal': float(p[i]),
                'residual': float(resid[i]),
                'r_biosample': float(r),
            })
    return rows


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--topn', type=int, default=20,
                    help='Top N candidates per direction per biosample.')
    ap.add_argument('--output', required=True, help='Output TSV path.')
    args = ap.parse_args()

    all_rows = []
    for label, slug in BIOSAMPLES:
        rows = biosample_outliers(label, slug, args.topn)
        all_rows.extend(rows)
        n_b = sum(1 for r in rows if r['direction'] == 'B_over')
        n_c = sum(1 for r in rows if r['direction'] == 'C_under')
        print(f'{label:8s}  r={rows[0]["r_biosample"]:.3f}  '
              f'B(top over) n={n_b}  C(top under) n={n_c}')

    cols = ['biosample', 'biosample_slug', 'direction', 'rank',
            'chrom', 'start', 'end',
            'exp_signal', 'pred_signal', 'residual', 'r_biosample']
    os.makedirs(os.path.dirname(os.path.abspath(args.output)), exist_ok=True)
    with open(args.output, 'w') as f:
        f.write('\t'.join(cols) + '\n')
        for r in all_rows:
            f.write('\t'.join(str(r[c]) for c in cols) + '\n')
    print(f'\nWrote {len(all_rows)} candidate rows → {args.output}')


if __name__ == '__main__':
    main()
