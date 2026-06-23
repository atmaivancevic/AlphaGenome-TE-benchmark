#!/usr/bin/env bash
#
# scripts/fig3_polymorphic_TE_eQTLs/render_fig3_AG_panels.sh
#
# Reproducible batch renderer for Fig 3 candidate AlphaGenome panels.
# Calls scripts/fig2_polymorphic_TE_example/plot_polymorphic_TE_insertion.py three times for 9 candidates,
# one call per chromatin-track configuration:
#
#   1. rna_only/         — RNA-seq tracks only
#   2. active_marks/     — RNA + H3K27ac (active enhancer mark) + ATAC + DNase
#   3. full_chromatin/   — RNA + 5 histone marks (H3K27ac, H3K4me1, H3K4me3,
#                          H3K27me3, H3K9me3) + ATAC + DNase
#
# Each configuration renders 3 plot types (baseline, diff, overlay) per variant,
# so the full output is 9 variants × 3 layouts × 3 plots = 81 PDFs, organized
# as figures/fig3_AG_predictions/{layout}/{variant}_{label}_{plot}.pdf.
#
# AG predictions are made fresh on each run (the script does not cache); each
# layout call requires one AG predict_variant call per variant (~5-10s each).
# Full batch ~10 min on a typical session.
#
# Usage:
#   bash scripts/fig3_polymorphic_TE_eQTLs/render_fig3_AG_panels.sh
#   bash scripts/fig3_polymorphic_TE_eQTLs/render_fig3_AG_panels.sh data/fig3_9_candidates.tab
#
# The default variant tab is data/fig3_9_candidates.tab; override as positional
# arg if running on a different candidate set.

set -euo pipefail

VARIANTS="${1:-data/fig3_9_candidates.tab}"
OUT_BASE="figures/fig3_AG_predictions"
ONTOLOGY="EFO:0002784"   # GM12878 LCL
LABEL="GM12878"

if [[ ! -f "$VARIANTS" ]]; then
  echo "ERROR: variant tab not found: $VARIANTS" >&2; exit 1
fi

mkdir -p "$OUT_BASE"

# Common args for each layout
COMMON=(--variants "$VARIANTS" --ontology "$ONTOLOGY" --label "$LABEL"
        --plots baseline overlay diff)

echo "[1/3] Rendering RNA-only layout..."
python3 scripts/fig2_polymorphic_TE_example/plot_polymorphic_TE_insertion.py "${COMMON[@]}" \
    --assays RNA_SEQ \
    --out-dir "$OUT_BASE/rna_only"

echo "[2/3] Rendering active-marks layout (RNA + H3K27ac + ATAC + DNase)..."
python3 scripts/fig2_polymorphic_TE_example/plot_polymorphic_TE_insertion.py "${COMMON[@]}" \
    --assays RNA_SEQ CHIP_HISTONE ATAC DNASE \
    --chip-marks H3K27ac \
    --out-dir "$OUT_BASE/active_marks"

echo "[3/3] Rendering full-chromatin layout (RNA + 5 histone marks + ATAC + DNase)..."
python3 scripts/fig2_polymorphic_TE_example/plot_polymorphic_TE_insertion.py "${COMMON[@]}" \
    --assays RNA_SEQ CHIP_HISTONE ATAC DNASE \
    --chip-marks H3K27ac H3K4me1 H3K4me3 H3K27me3 H3K9me3 \
    --out-dir "$OUT_BASE/full_chromatin"

echo
echo "=== Done ==="
for layout in rna_only active_marks full_chromatin; do
    n=$(ls "$OUT_BASE/$layout"/*.pdf 2>/dev/null | wc -l | tr -d ' ')
    echo "  $OUT_BASE/$layout/ : $n PDFs"
done
