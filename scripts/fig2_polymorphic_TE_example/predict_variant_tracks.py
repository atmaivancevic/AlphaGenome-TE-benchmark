"""
Predict AlphaGenome track signals for both REF and ALT of each variant, and
extract the mean + max signal over a user-specified locus per track. Used
for the Fig 6 AP1-perturbation panels where we want to see the predicted
H3K27ac / FOSL1 / RNA peak height shrink as AP1 motifs are scrambled —
i.e. the absolute predicted signal per allele, not the variant log2FC.

Distinct from score_variant_*.py:
  - score_variant_*.py uses AG's score_variant() to get a single scalar
    log2FC per (variant, scorer, biosample) tuple. Useful for the "how big
    is the effect?" question.
  - This script uses AG's predict_variant() to get the FULL track arrays
    for both REF (baseline / WT-side) and ALT (perturbed side). Useful for
    the "how does the predicted peak look on each allele?" question.

For each variant row in --variants, the script runs predict_variant() once
and then extracts, per requested track:
    - signal_mean_REF, signal_max_REF       (track values averaged across
    - signal_mean_ALT, signal_max_ALT        bins overlapping the locus,
                                             for REF and ALT separately)
The locus is taken from the per-variant POS / REF length (the variant's own
genomic span), unless --locus is set to a fixed chrom:start-end window.

Outputs:
  --out-csv  one row per (variant, track) combo with mean/max signal for
             both REF and ALT alleles + variant metadata.
  Optionally --save-tracks-dir to dump per-variant full track arrays as
  numpy .npz files (large; useful for genome-browser-style visualisation).

Inputs:
  --variants  tab file (ID/CHROM/POS/REF/ALT). Supports VCF-biallelic,
              CRISPR-synthetic DEL (ALT='.'), and MNV (REF and ALT same
              length, e.g. Fig 6 AP1 perturbations).
  --ontology  cell-type ontology term, default EFO:0002824 (HCT116).
  --requested-outputs
              comma-separated list of AG OutputType names to request, e.g.
              'CHIP_HISTONE,CHIP_TF,RNA_SEQ'. Default
              'CHIP_HISTONE,CHIP_TF,RNA_SEQ'.
  --keep-tracks
              optional list of track name substrings to filter the output.
              e.g. 'H3K27ac FOSL1' to keep only those two ChIP tracks.
              Applied per output_type, case-insensitive.

Usage:
    python scripts/fig2_polymorphic_TE_example/predict_variant_tracks.py \\
        --variants  data/LTR10_ATG12_AP1_perturbations.tab \\
        --ontology  EFO:0002824 \\
        --cell-line HCT116 \\
        --keep-tracks H3K27ac FOSL1 total RNA-seq \\
        --out-csv   results/AG_predict_variant_tracks_LTR10_ATG12_AP1.csv
"""
import argparse
import os
from pathlib import Path
import numpy as np
import pandas as pd

import alphagenome
from alphagenome.data import genome
from alphagenome.models import dna_client
from alphagenome.models.dna_output import OutputType

p = argparse.ArgumentParser()
p.add_argument('--variants', required=True)
p.add_argument('--ontology', default='EFO:0002824',
               help='Cell type ontology term (default EFO:0002824 = HCT116)')
p.add_argument('--cell-line', default='HCT116',
               help='Used post-hoc to filter AG tracks by biosample_name')
p.add_argument('--requested-outputs', default='CHIP_HISTONE,CHIP_TF,RNA_SEQ',
               help='Comma-separated AG OutputType names to request')
p.add_argument('--keep-tracks', nargs='*', default=None,
               help='Optional list of track-name substrings to keep (e.g. H3K27ac FOSL1)')
p.add_argument('--locus', default=None,
               help='Optional fixed chrom:start-end window (1-based inclusive) to score signal over. '
                    "Default: each variant's own genomic span (POS .. POS+len(REF)-1).")
p.add_argument('--api-key', default='scripts/my_api_key.txt')
p.add_argument('--out-csv', required=True)
p.add_argument('--save-tracks-dir', default=None,
               help='If set, write per-variant track arrays as .npz files to this directory')
p.add_argument('--test', action='store_true', help='Run only the first variant')
args = p.parse_args()

with open(args.api_key) as f:
    api_key = f.read().strip()
model = dna_client.create(api_key)

# Parse requested outputs
requested_outputs = {
    OutputType[name.strip()] for name in args.requested_outputs.split(',') if name.strip()
}
print(f"AlphaGenome version: {alphagenome.__version__}")
print(f"Cell line: {args.cell_line}  |  Ontology: {args.ontology}")
print(f"Requested output types: {[ot.name for ot in requested_outputs]}")
if args.keep_tracks:
    print(f"Track-name filter: {args.keep_tracks}")

# Optional fixed locus override
fixed_locus = None
if args.locus:
    chrom, span = args.locus.split(':')
    start, end = (int(x) for x in span.split('-'))
    fixed_locus = (chrom, start, end)
    print(f"Locus override: {chrom}:{start:,}-{end:,}  (used for ALL variants)")

variant_df = pd.read_table(args.variants, sep='\t')
if args.test:
    variant_df = variant_df.head(1)
print(f"Loaded {len(variant_df)} variants from {args.variants}")

out_dir = None
if args.save_tracks_dir:
    out_dir = Path(args.save_tracks_dir); out_dir.mkdir(parents=True, exist_ok=True)

def signal_over_locus(track_arr: np.ndarray, interval_start: int, interval_end: int,
                      locus_start: int, locus_end: int):
    """Mean and max of track_arr over [locus_start, locus_end].

    track_arr is shape (num_bins, num_tracks). interval_start/end define the
    genomic span of the 1 Mb AG window; AG uses 128 bp resolution by default.
    We map locus coords to bin indices and take the mean/max per track.
    """
    n_bins = track_arr.shape[0]
    bin_size = (interval_end - interval_start) // n_bins
    lo_bin = max(0, (locus_start - interval_start) // bin_size)
    hi_bin = min(n_bins, (locus_end - interval_start) // bin_size + 1)
    sub = track_arr[lo_bin:hi_bin]
    return sub.mean(axis=0), sub.max(axis=0)

rows = []
for i, r in variant_df.iterrows():
    vid = r['ID']
    chrom, pos = r['CHROM'], int(r['POS'])
    var_ref = str(r['REF']).upper()
    var_alt_raw = '' if str(r['ALT']) == '.' else str(r['ALT']).upper()

    # AG Variant: handle MNV (same-length), INS, DEL, and CRISPR-synthetic-DEL (alt='').
    variant = genome.Variant(
        chromosome=chrom, position=pos,
        reference_bases=var_ref, alternate_bases=var_alt_raw,
        name=vid,
    )
    interval = variant.reference_interval.resize(dna_client.SEQUENCE_LENGTH_1MB)
    print(f"[{i+1}/{len(variant_df)}] {vid}  REF={len(var_ref)} bp  ALT={len(var_alt_raw)} bp")

    out = model.predict_variant(
        interval=interval,
        variant=variant,
        ontology_terms=[args.ontology],
        requested_outputs=requested_outputs,
        organism=dna_client.Organism.HOMO_SAPIENS,
    )
    # AG returns a VariantOutput with .reference and .alternate, each an
    # Output container with attributes per requested output type
    # (e.g. .chip_histone, .chip_tf, .rna_seq). Each is a TrackData-like
    # object with .values (np.ndarray, shape (bins, tracks)) and .metadata.
    int_start = interval.start
    int_end   = interval.end
    locus_start, locus_end = (
        fixed_locus[1], fixed_locus[2]
    ) if fixed_locus else (pos, pos + len(var_ref) - 1)

    # Attribute name on the Output container per OutputType
    OUT_ATTR = {
        OutputType.CHIP_HISTONE: 'chip_histone',
        OutputType.CHIP_TF:      'chip_tf',
        OutputType.ATAC:         'atac',
        OutputType.RNA_SEQ:      'rna_seq',
        OutputType.DNASE:        'dnase',
    }

    for ot in requested_outputs:
        attr = OUT_ATTR.get(ot, ot.name.lower())
        ref_track = getattr(out.reference, attr, None)
        alt_track = getattr(out.alternate, attr, None)
        if ref_track is None or alt_track is None:
            print(f"    (no {attr} output for {vid})"); continue

        ref_vals = ref_track.values  # (bins, tracks)
        alt_vals = alt_track.values
        track_meta = ref_track.metadata  # pandas DataFrame per AG SDK

        # Track-name filter
        if args.keep_tracks:
            mask = np.zeros(len(track_meta), dtype=bool)
            for needle in args.keep_tracks:
                mask |= track_meta['name'].astype(str).str.contains(
                    needle, case=False, na=False).values
            if not mask.any():
                continue
            ref_vals = ref_vals[:, mask]
            alt_vals = alt_vals[:, mask]
            track_meta = track_meta[mask].reset_index(drop=True)

        # Cell-line filter (if biosample_name in metadata)
        if 'biosample_name' in track_meta.columns:
            cl_mask = track_meta['biosample_name'].astype(str) == args.cell_line
            if cl_mask.any():
                ref_vals = ref_vals[:, cl_mask.values]
                alt_vals = alt_vals[:, cl_mask.values]
                track_meta = track_meta[cl_mask].reset_index(drop=True)

        # Compute signal mean/max over locus, per track
        ref_mean, ref_max = signal_over_locus(ref_vals, int_start, int_end, locus_start, locus_end)
        alt_mean, alt_max = signal_over_locus(alt_vals, int_start, int_end, locus_start, locus_end)

        for j, tr in track_meta.iterrows():
            rows.append({
                'variant_id':       vid,
                'chrom':            chrom,
                'pos':              pos,
                'ref_len':          len(var_ref),
                'alt_len':          len(var_alt_raw),
                'locus_start':      locus_start,
                'locus_end':        locus_end,
                'output_type':      ot.name,
                'track_name':       tr.get('name', ''),
                'biosample_name':   tr.get('biosample_name', ''),
                'signal_mean_REF':  float(ref_mean[j]),
                'signal_max_REF':   float(ref_max[j]),
                'signal_mean_ALT':  float(alt_mean[j]),
                'signal_max_ALT':   float(alt_max[j]),
            })

    if out_dir is not None:
        np.savez(out_dir / f'{vid}_tracks.npz',
                 ref_chip_histone=getattr(out.reference, 'chip_histone', None).values if hasattr(out.reference, 'chip_histone') and getattr(out.reference, 'chip_histone') is not None else np.empty(0),
                 alt_chip_histone=getattr(out.alternate, 'chip_histone', None).values if hasattr(out.alternate, 'chip_histone') and getattr(out.alternate, 'chip_histone') is not None else np.empty(0),
                 interval_start=int_start, interval_end=int_end)

# Write output
out_df = pd.DataFrame(rows)
Path(args.out_csv).parent.mkdir(parents=True, exist_ok=True)
out_df.to_csv(args.out_csv, index=False)
print(f"\nWrote {args.out_csv}  ({len(out_df)} (variant, track) rows)")
