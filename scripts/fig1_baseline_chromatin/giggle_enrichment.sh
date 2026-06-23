#!/bin/bash

## GIGGLE enrichment of H3K27ac peaks against hg38 repeat families.
## Outputs per-biosample tables, a combined ranked table, and a filtered table
## of the strongest TE subfamilies. Run on a Linux HPC (needs giggle + bgzip on
## PATH); the `module load samtools` line below is HPC-specific.
##
## Example usage:
## bash scripts/fig1_baseline_chromatin/giggle_enrichment.sh \
##   data/encode_h3k27ac/peaks/merged \
##   data/encode_h3k27ac/giggle_results \
##   /path/to/giggle/hg38/repeats/indexed

module load samtools

PEAKDIR="${1:-data/encode_h3k27ac/peaks/merged}"
# Default output dir = where plot_giggle_bubbles.py reads its inputs.
OUTDIR="${2:-data/encode_h3k27ac/giggle_results}"
GIGGLE_INDEX="${3:-giggle/hg38/repeats/indexed}"
GENOME_SIZE=3209286105   # hg38 genome size (bp), for giggle -g

mkdir -p "$OUTDIR"

# Step 1: bgzip peak files for GIGGLE
for f in "$PEAKDIR"/*.narrowPeak.gz; do
  base=$(basename "$f" .narrowPeak.gz)
  zcat "$f" | cut -f1-3 | sort -k1,1 -k2,2n | bgzip > "$OUTDIR/${base}.bed.gz"
done

# Step 2: run GIGGLE search for each biosample
for f in "$OUTDIR"/*.bed.gz; do
  echo "$f"
  giggle search -q "$f" -i "$GIGGLE_INDEX" -s -g $GENOME_SIZE \
    | sed 's#sorted/##g' \
    | sed 's/.bed.gz//g' \
    | grep -v "#" \
    | sort -nrk8 \
    | sed '1i#file\tfile_size\toverlaps\todds_ratio\tfishers_two_tail\tfishers_left_tail\tfishers_right_tail\tcombo_score' \
    > "${f%.gz}.VSrepeats.tab"
done

# Step 3: combine all biosamples into one table
for f in "$OUTDIR"/*.VSrepeats.tab; do
  base=$(basename "$f" .bed.VSrepeats.tab)
  grep -v "#" "$f" \
    | awk -v b="$base" '{print b "\t" $0}'
done \
  | sort -nrk9 \
  | awk '{print $1 "\t" $2 "\t" $3 "\t" $4 "\t" $5 "\t" $9}' \
  | sed '1ibiosample\trepeat\tfilesize\toverlaps\toddsratio\tgigglescore' \
  > "$OUTDIR/h3k27ac_vs_TEs_rankedByScore.tab"

echo "[done] $(wc -l < "$OUTDIR/h3k27ac_vs_TEs_rankedByScore.tab") rows"

# Step 4: filtered table of the strongest TE subfamilies (odds ratio >5,
# GIGGLE score >=100, overlaps >40); plot_giggle_bubbles.py uses it for the default TE set.
awk -F'\t' 'NR==1 || ($4>40 && $5>5 && $6>=100)' \
  "$OUTDIR/h3k27ac_vs_TEs_rankedByScore.tab" \
  > "$OUTDIR/filtered_OR5_score100_overlaps40.tab"
echo "[done] filtered: $(($(wc -l < "$OUTDIR/filtered_OR5_score100_overlaps40.tab") - 1)) rows, $(tail -n +2 "$OUTDIR/filtered_OR5_score100_overlaps40.tab" | cut -f2 | sort -u | wc -l | tr -d ' ') TE subfamilies"
