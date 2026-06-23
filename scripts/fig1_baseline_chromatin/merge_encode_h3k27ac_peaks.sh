#!/usr/bin/env bash
# Merge narrowPeak BEDs for multi-experiment AG training tracks.
# For each H3K27ac biosample in the AG training manifest:
#   - if multiple ENCODE experiments map to one biosample (e.g. 38 of 159
#     in the full human set), concat + bedtools merge averaging score
#     columns
#   - if one experiment, pass through with a numeric (chrom,start) sort
#     applied (ENCODE narrowPeaks are typically lex-sorted on POS, which
#     would otherwise collapse the 1Mb-cache hit rate in
#     scripts/fig1_baseline_chromatin/predict_peak_signal_batched.py)
#
# Output: one 10-column narrowPeak per biosample.
# Manifest defaults to the full Avsec 2026 H3K27ac track list (159 human +
# 31 mouse) — pass arg 4 to override (e.g., the PASS/WARNING-filtered
# 133-track subset used in the original Apr 2026 run).
set -euo pipefail

ENRICHED="${1:-data/encode_h3k27ac/peak_resolution_human.tsv}"
PEAKDIR="${2:-data/encode_h3k27ac/peaks/human}"
OUTDIR="${3:-data/encode_h3k27ac/peaks/merged}"
MANIFEST="${4:-data/AG_training_H3K27ac_tracks.csv}"
mkdir -p "$OUTDIR"

SCRIPT_DIR="$(dirname "$0")"

python3 "$SCRIPT_DIR/group_peaks_by_biosample.py" \
    --enriched "$ENRICHED" \
    --manifest "$MANIFEST" \
    --peakdir  "$PEAKDIR" \
  | while IFS=$'\t' read -r label files; do
    IFS=$'\t' read -ra PEAKS <<< "$files"
    OUT="$OUTDIR/${label}.narrowPeak.gz"

    if [ ${#PEAKS[@]} -eq 0 ]; then
        echo "[skip] $label: no peak files"
        continue
    fi

    if [ ${#PEAKS[@]} -eq 1 ]; then
        gzcat "${PEAKS[0]}" | sort -k1,1 -k2,2n | gzip > "$OUT"
        echo "[copy+sort] $label"
        continue
    fi

    # Multi-experiment: concat → sort → bedtools merge (average scores) → 10-col narrowPeak
    gzcat "${PEAKS[@]}" \
      | sort -k1,1 -k2,2n \
      | bedtools merge -i - -c 5,7,8,9 -o mean,mean,mean,mean \
      | awk -v OFS='\t' '{
          printf "%s\t%s\t%s\tMergedPeak_%d\t%d\t.\t%.5f\t%.5f\t%.5f\t-1\n",
            $1, $2, $3, NR, int($4+0.5), $5, $6, $7
        }' \
      | gzip > "$OUT"
    echo "[merge] $label: ${#PEAKS[@]} files"
done

echo "[done] $(ls "$OUTDIR"/*.narrowPeak.gz 2>/dev/null | wc -l | tr -d ' ') merged peak files in $OUTDIR"
