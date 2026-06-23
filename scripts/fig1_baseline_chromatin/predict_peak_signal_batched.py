"""
Predict AlphaGenome H3K27ac signal at experimental peak regions. For each peak,
centers a 1Mb window on the peak midpoint, caches the prediction, and reuses it
for any later peak whose midpoint falls in the window (much faster than one API
call per peak). Input must be sorted by chrom, start.

Example usage:
python scripts/fig1_baseline_chromatin/predict_peak_signal_batched.py \
    --regions data/encode_h3k27ac/peaks/merged/EFO_0002824_HCT116.narrowPeak.gz \
    --ontology EFO:0002824 \
    --label EFO_0002824_HCT116 \
    --outdir results/AG_predicted_h3k27ac_batched
"""

import os, time, argparse, gzip
import numpy as np

parser = argparse.ArgumentParser()
parser.add_argument('--regions', required=True,
                    help='BED or narrowPeak file (gzipped ok)')
parser.add_argument('--ontology', required=True,
                    help='Biosample ontology ID (e.g. EFO:0002824)')
parser.add_argument('--label', required=True,
                    help='Output label (e.g. EFO_0002824_HCT116)')
parser.add_argument('--outdir', required=True,
                    help='Output directory')
parser.add_argument('--test', action='store_true',
                    help='First 20 regions only')
args = parser.parse_args()

# ── 1. AlphaGenome setup ────────────────────────────────────────────────────

SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
# Script lives in scripts/fig1_baseline_chromatin/, so repo root is two levels up.
PROJECT_ROOT = os.path.dirname(os.path.dirname(SCRIPT_DIR))

from alphagenome.data import genome
from alphagenome.models import dna_client

with open(os.path.join(PROJECT_ROOT, 'scripts', 'my_api_key.txt')) as f:
    api_key = f.read().strip()

model = dna_client.create(api_key)
OUTPUTS = {dna_client.OutputType.CHIP_HISTONE}
WINDOW = dna_client.SEQUENCE_LENGTH_1MB  # 1,048,576

# ── 2. Load regions ─────────────────────────────────────────────────────────

opener = gzip.open if args.regions.endswith('.gz') else open
with opener(args.regions, 'rt') as f:
    regions = []
    for line in f:
        if line.startswith('#') or line.strip() == '':
            continue
        cols = line.strip().split('\t')
        chrom, start, end = cols[0], int(cols[1]), int(cols[2])
        name = cols[3] if len(cols) > 3 else f'{chrom}:{start}-{end}'
        regions.append((chrom, start, end, name))

if args.test:
    regions = regions[:20]

# Verify input is sorted by chrom, start. Unsorted input would silently collapse the cache hit rate.
last_chrom = None
last_start = -1
for chrom, start, _, _ in regions:
    if last_chrom is not None:
        if chrom < last_chrom:
            raise ValueError(f"Input not sorted by chromosome: {last_chrom} -> {chrom}")
        if chrom == last_chrom and start < last_start:
            raise ValueError(f"Input not sorted by start within {chrom}: "
                             f"{last_start} -> {start}")
    last_chrom = chrom
    last_start = start

print(f"Loaded {len(regions)} regions | {args.ontology} | {args.label}")

# ── 3. Signal extraction helper ─────────────────────────────────────────────

def extract_h3k27ac_max(h_meta, h_vals, pred_start, res, pk_start, pk_end):
    """Max H3K27ac signal across bins overlapping a peak.

    Asserts that AG returns exactly one H3K27ac track for the requested
    ontology (current AG schema invariant: every ontology with H3K27ac
    has exactly 1 track; verified across all 159 H3K27ac ontologies as
    of 2026-05-19). If AG ever ships a multi-replicate schema, this
    raises rather than silently picking the first.
    """
    h3k_idx = np.where(h_meta['histone_mark'].values == 'H3K27ac')[0]
    if len(h3k_idx) != 1:
        raise ValueError(
            f"Expected exactly 1 H3K27ac track in AG output; got {len(h3k_idx)}. "
            f"Aggregation semantics must be decided before proceeding.")
    vals = h_vals[:, h3k_idx[0]]
    s = max(0, int(np.floor((pk_start - pred_start) / res)))
    e = max(s + 1, int(np.ceil((pk_end - pred_start) / res)))
    e = min(e, vals.shape[0])
    if s >= e:
        return np.nan
    return float(vals[s:e].max())

# ── 4. Score each region (with window caching) ──────────────────────────────

OUT_DIR = os.path.join(args.outdir, args.label)
os.makedirs(OUT_DIR, exist_ok=True)

signals = np.full(len(regions), np.nan, dtype=np.float32)
n_api_calls = 0
errors = []  # (index, chrom, start, end, name, error_message)
t_start = time.time()

# Cache: list of (chrom, window_start, window_end, h_meta, h_vals, pred_start, res)
cache = []

for i, (chrom, start, end, name) in enumerate(regions):
    center = (start + end) // 2

    # Check if this peak falls entirely within a cached window
    hit = None
    for cached in cache:
        if (cached[0] == chrom and
                cached[1] <= start and end <= cached[2]):
            hit = cached
            break

    if hit is None:
        # New API call, centered on this peak
        interval = genome.Interval(chrom, center, center + 1).resize(WINDOW)
        try:
            out = model.predict_interval(
                interval, requested_outputs=OUTPUTS,
                ontology_terms=[args.ontology])
            n_api_calls += 1

            h_meta = out.chip_histone.metadata
            h_vals = out.chip_histone.values
            pred_start = out.chip_histone.interval.start
            res = out.chip_histone.resolution

            # Cache using model-returned bounds (may differ at chrom edges)
            pred_end = pred_start + h_vals.shape[0] * res
            cache.append((chrom, pred_start, pred_end,
                          h_meta, h_vals, pred_start, res))
            # Keep only the last 5 windows (input is position-sorted)
            cache = cache[-5:]

            signals[i] = extract_h3k27ac_max(
                h_meta, h_vals, pred_start, res, start, end)

        except (RuntimeError, ValueError, ConnectionError, TimeoutError, OSError) as e:
            errors.append((i, chrom, start, end, name, str(e)))
            n_api_calls += 1
            time.sleep(5)
    else:
        # Reuse cached prediction
        _, _, _, h_meta, h_vals, pred_start, res = hit
        signals[i] = extract_h3k27ac_max(
            h_meta, h_vals, pred_start, res, start, end)

    if (i + 1) % 500 == 0 or i == 0:
        eta = (time.time() - t_start) / (i + 1) * (len(regions) - i - 1)
        print(f"  [{i+1}/{len(regions)}] {chrom}:{start}-{end}  "
              f"sig={signals[i]:.1f}  calls={n_api_calls}  ETA {eta/60:.0f}min",
              flush=True)

n_failed = np.isnan(signals).sum()
if n_failed > 0:
    print(f"\n{n_failed} peaks failed. Retrying...")

# ── 5. Retry failed peaks ──────────────────────────────────────────────────

retry_errors = []
for idx, chrom, start, end, name, orig_err in errors:
    # Check cache first (a previous retry may cover this peak)
    hit = None
    for cached in cache:
        if (cached[0] == chrom and
                cached[1] <= start and end <= cached[2]):
            hit = cached
            break

    if hit is not None:
        _, _, _, h_meta, h_vals, pred_start, res = hit
        signals[idx] = extract_h3k27ac_max(
            h_meta, h_vals, pred_start, res, start, end)
        retry_errors.append((chrom, start, end, name, orig_err, 'resolved (from cache)'))
        continue

    center = (start + end) // 2
    interval = genome.Interval(chrom, center, center + 1).resize(WINDOW)
    try:
        out = model.predict_interval(
            interval, requested_outputs=OUTPUTS,
            ontology_terms=[args.ontology])
        n_api_calls += 1

        h_meta = out.chip_histone.metadata
        h_vals = out.chip_histone.values
        pred_start = out.chip_histone.interval.start
        res = out.chip_histone.resolution

        signals[idx] = extract_h3k27ac_max(
            h_meta, h_vals, pred_start, res, start, end)
        retry_errors.append((chrom, start, end, name, orig_err, 'resolved'))

        # Cache retry result for nearby failed peaks
        pred_end = pred_start + h_vals.shape[0] * res
        cache.append((chrom, pred_start, pred_end,
                      h_meta, h_vals, pred_start, res))
        cache = cache[-5:]

    except (RuntimeError, ValueError, ConnectionError, TimeoutError, OSError) as e:
        signals[idx] = np.nan
        retry_errors.append((chrom, start, end, name, orig_err, f'failed again: {e}'))
        n_api_calls += 1
        time.sleep(5)

# Write error log
if retry_errors:
    log_file = os.path.join(OUT_DIR, f'{args.label}_errors.tsv')
    with open(log_file, 'w') as f:
        f.write('chrom\tstart\tend\tname\toriginal_error\tretry_status\n')
        for row in retry_errors:
            f.write('\t'.join(str(x) for x in row) + '\n')
    n_resolved = sum(1 for r in retry_errors if r[5].startswith('resolved'))
    n_still_failed = len(retry_errors) - n_resolved
    print(f"  {n_resolved} resolved on retry, {n_still_failed} still failed")
    print(f"  Error log: {log_file}")

# ── 6. Write narrowPeak output (excluding unresolved peaks) ────────────────

n_input = len(regions)
out_file = os.path.join(OUT_DIR, f'{args.label}.narrowPeak.gz')
n_output = 0

with gzip.open(out_file, 'wt') as f:
    for i, (chrom, start, end, name) in enumerate(regions):
        if np.isnan(signals[i]):
            continue
        sig = float(signals[i])
        score = min(int(sig * 10), 1000)
        summit = (end - start) // 2
        f.write(f'{chrom}\t{start}\t{end}\t{name}\t{score}\t.\t{sig:.4f}\t-1\t-1\t{summit}\n')
        n_output += 1

# Write missed peaks to error file if counts don't match
n_missed = n_input - n_output
if n_missed > 0:
    missed_file = os.path.join(OUT_DIR, f'{args.label}_missed.bed')
    with open(missed_file, 'w') as f:
        for i, (chrom, start, end, name) in enumerate(regions):
            if np.isnan(signals[i]):
                f.write(f'{chrom}\t{start}\t{end}\t{name}\n')
    print(f"\nWARNING: {n_missed} peaks excluded (no prediction)")
    print(f"  Missed peaks: {missed_file}")

print(f"\nInput: {n_input} peaks | Output: {n_output} peaks | API calls: {n_api_calls}")
if n_input == n_output:
    print("All peaks scored successfully.")
print(f"Output: {out_file}")
print(f"Total time: {(time.time() - t_start)/60:.1f} min")
