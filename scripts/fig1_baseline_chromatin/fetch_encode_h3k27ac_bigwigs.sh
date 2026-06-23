#!/bin/bash

## Download ENCODE H3K27ac bigWigs by accession, verifying md5 from
## the ENCODE REST API. BigWigs are large, so only the biosamples we need are
## fetched. Output: data/encode_h3k27ac/bigwigs/human/<ENCSR>_<ENCFF>.bigWig
##
## Example usage (HCT116 + GM12878):
## bash scripts/fig1_baseline_chromatin/fetch_encode_h3k27ac_bigwigs.sh ENCFF169MCH ENCFF469WVA
set -euo pipefail

OUT="${OUT:-data/encode_h3k27ac/bigwigs/human}"
mkdir -p "$OUT"

if [ $# -eq 0 ]; then
  echo "usage: $0 <ENCFF...>" >&2
  exit 1
fi

i=0; n=$#
for acc in "$@"; do
  i=$((i+1))
  META="https://www.encodeproject.org/files/${acc}/?format=json"
  meta_json=$(curl -sfL -H "Accept: application/json" "$META")
  expected_md5=$(printf '%s' "$meta_json" | python3 -c 'import sys,json;print(json.load(sys.stdin)["md5sum"])')
  encsr=$(printf '%s' "$meta_json" | python3 -c 'import sys,json,os;d=json.load(sys.stdin);print((d.get("dataset") or "").rstrip("/").rsplit("/",1)[-1])')
  ftype=$(printf '%s' "$meta_json" | python3 -c 'import sys,json;d=json.load(sys.stdin);print(d.get("file_format",""))')
  if [ "$ftype" != "bigWig" ]; then
    echo "[skip $i/$n] $acc: file_format='$ftype' (expected bigWig)" >&2
    continue
  fi

  FN="${encsr}_${acc}.bigWig"
  DEST="$OUT/$FN"

  if [ -s "$DEST" ]; then
    have=$(md5 -q "$DEST" 2>/dev/null || md5sum "$DEST" | awk '{print $1}')
    if [ "$have" = "$expected_md5" ]; then
      echo "[cache $i/$n] $FN"
      continue
    fi
    echo "[restale $i/$n] $FN md5 mismatch — redownloading"
  fi

  URL="https://www.encodeproject.org/files/${acc}/@@download/${acc}.bigWig"
  echo "[get $i/$n] $FN ← $URL"
  curl -fL --progress-bar "$URL" -o "$DEST"
  have=$(md5 -q "$DEST" 2>/dev/null || md5sum "$DEST" | awk '{print $1}')
  if [ "$have" = "$expected_md5" ]; then
    size=$(du -h "$DEST" | awk '{print $1}')
    echo "[ok $i/$n] $FN ($size, md5=$have)"
  else
    echo "[md5-FAIL $i/$n] $FN  expected=$expected_md5  got=$have" >&2
    exit 1
  fi
done
