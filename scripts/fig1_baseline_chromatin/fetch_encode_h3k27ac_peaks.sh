#!/bin/bash

## Download the ENCODE H3K27ac narrowPeaks resolved by
## resolve_encode_h3k27ac_peaks.sh, and verify md5sums.
##
## Example usage:
## bash scripts/fig1_baseline_chromatin/fetch_encode_h3k27ac_peaks.sh
set -euo pipefail

RES="${1:-data/encode_h3k27ac/peak_resolution_human.tsv}"
OUT="${2:-data/encode_h3k27ac/peaks/human}"
mkdir -p "$OUT"

N=$(awk -F'\t' 'NR>1 && $4!=""' "$RES" | wc -l | tr -d ' ')
echo "[info] $N peak files to fetch"

i=0; ok=0; fail=0
awk -F'\t' 'NR>1 && $4!="" {print $1"\t"$4"\t"$10}' "$RES" | while IFS=$'\t' read -r encsr pk md5; do
  i=$((i+1))
  FN="${encsr}_${pk}.bed.gz"
  DEST="$OUT/$FN"
  if [ -s "$DEST" ]; then
    have=$(md5 -q "$DEST" 2>/dev/null || md5sum "$DEST" | awk '{print $1}')
    if [ "$have" = "$md5" ]; then
      echo "[cache $i/$N] $FN"
      continue
    fi
    echo "[restale $i/$N] $FN md5 mismatch — redownloading"
  fi
  URL="https://www.encodeproject.org/files/${pk}/@@download/${pk}.bed.gz"
  if curl -sfL "$URL" -o "$DEST"; then
    have=$(md5 -q "$DEST" 2>/dev/null || md5sum "$DEST" | awk '{print $1}')
    if [ "$have" = "$md5" ]; then
      echo "[ok $i/$N] $FN"
    else
      echo "[md5-FAIL $i/$N] $FN  expected=$md5  got=$have" >&2
    fi
  else
    echo "[err $i/$N] $FN: curl failed" >&2
  fi
done

echo "[done] $(ls "$OUT"/*.bed.gz 2>/dev/null | wc -l | tr -d ' ') files in $OUT"
