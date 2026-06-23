"""
Predict + plot the per-track effect of a polymorphic TE insertion using
AlphaGenome. The companion of helper_functions.plot_difference (deletions),
adapted for insertion variants from Schloissnig 2025 Supp Table 2.

Default usage scores all variants in --variants against the GM12878 ontology
(EFO:0002784) for ChIP-histone, ChIP-TF, and RNA-seq tracks, and saves one
PDF per variant to figures/.

Insertion variant construction:
    reference_bases = single-base anchor at POS (REF column from VCF)
    alternate_bases = anchor + inserted sequence (ALT column from VCF)

Usage:
    python scripts/fig2_polymorphic_TE_example/plot_polymorphic_TE_insertion.py \\
        --variants data/polymorphicTE_variants.tab \\
        --ontology EFO:0002784 --label GM12878
"""
import argparse, os, sys
import matplotlib.pyplot as plt
import pandas as pd

import numpy as np
from alphagenome.data import gene_annotation, genome, track_data, transcript as transcript_utils
from alphagenome.models import dna_client
from alphagenome.visualization import plot_components

p = argparse.ArgumentParser()
p.add_argument('--variants', required=True, help='Variant tab file (ID, CHROM, POS, REF, ALT, ...)')
p.add_argument('--ontology', required=True, help='AG ontology term (e.g. EFO:0002784 for GM12878)')
p.add_argument('--label', required=True, help='Cell line label for figure titles + filenames (e.g. GM12878)')
p.add_argument('--assays', nargs='+', default=['CHIP_HISTONE', 'RNA_SEQ'],
               help='AG OutputType names to score. CHIP_TF is excluded by default '
                    '(GM12878 has 101 TF tracks, exceeds the AG plot cap of 50). '
                    'Pass an explicit subset of TF tracks (out of scope here) to include.')
p.add_argument('--flank-bp', type=int, default=5000, help='Flank around insertion site for visualization')
p.add_argument('--api-key', default='scripts/my_api_key.txt')
p.add_argument('--gtf-feather', default='data/gencode.v46.annotation.feather',
               help='Local cache for AG\'s pre-built GENCODE v46 transcript annotation (Feather format). '
                    'Auto-downloaded from `storage.googleapis.com/alphagenome/reference/gencode/hg38/...` if missing.')
p.add_argument('--out-dir', default='figures', help='Output directory for PDFs')
p.add_argument('--frame', choices=['alt', 'ref'], default='alt',
               help='Coordinate frame for plotting. "alt" (default) extends the x-axis by insert_len '
                    'at the variant site so the inserted sequence\'s own signal is visible (alt-genome frame); '
                    'ref shows zeros there. "ref" splices out the insert\'s rows and aligns downstream onto '
                    'reference coordinates (insert signal lost).')
p.add_argument('--chip-marks', nargs='+', default=None,
               help='If set, keep only ChIP_HISTONE tracks whose name contains one of these substrings '
                    '(e.g. --chip-marks H3K27ac H3K36me3 H3K4me1). Default: all marks.')
p.add_argument('--rna-track-filter', nargs='+', default=None,
               help='If set, keep only RNA_SEQ tracks whose name contains one of these substrings '
                    '(e.g. --rna-track-filter "polyA plus" or "total RNA-seq"). GM12878 has both polyA '
                    'and total RNA tracks; this flag picks one flavor for cleaner publication panels.')
p.add_argument('--plots', nargs='+', default=['baseline', 'diff', 'overlay'],
               choices=['baseline', 'diff', 'overlay'],
               help='Which plot types to generate (default: all three).')
args = p.parse_args()

# Setup
SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
# Script lives in scripts/fig2_polymorphic_TE_example/, so the repo root is two
# levels up (scripts/fig2_polymorphic_TE_example -> scripts -> repo root).
PROJECT_ROOT = os.path.dirname(os.path.dirname(SCRIPT_DIR))
out_dir = os.path.join(PROJECT_ROOT, args.out_dir)
os.makedirs(out_dir, exist_ok=True)

with open(os.path.join(PROJECT_ROOT, args.api_key)) as f:
    api_key = f.read().strip()
dna_model = dna_client.create(api_key)

# GTF for transcript annotation panel — AG-distributed Feather file (GENCODE v46)
import urllib.request
gtf_path = os.path.join(PROJECT_ROOT, args.gtf_feather)
if not os.path.isfile(gtf_path):
    os.makedirs(os.path.dirname(gtf_path), exist_ok=True)
    gtf_url = 'https://storage.googleapis.com/alphagenome/reference/gencode/hg38/gencode.v46.annotation.gtf.gz.feather'
    print(f'Downloading GENCODE v46 Feather from {gtf_url} ...')
    urllib.request.urlretrieve(gtf_url, gtf_path)
    print(f'  cached at {gtf_path}')

gtf = pd.read_feather(gtf_path)
# protein-coding + longest-transcript-per-gene only. Dropping TSL filter so that
# retrogenes / newer protein-coding loci (e.g. RPSA2 — all 13 transcripts at TSL>1
# or NA, would otherwise be filtered out) still show up in the transcript panel.
gtf_longest = gene_annotation.filter_to_longest_transcript(
    gene_annotation.filter_protein_coding(gtf)
)
transcript_extractor = transcript_utils.TranscriptExtractor(gtf_longest)
print(f'Loaded GTF (longest protein-coding transcript per gene): {len(gtf_longest):,} rows')

variant_df = pd.read_table(args.variants)
print(f'Loaded {len(variant_df)} variants from {args.variants}')

assay_attr_map = {
    'RNA_SEQ': 'rna_seq', 'CHIP_HISTONE': 'chip_histone', 'CHIP_TF': 'chip_tf',
    'DNASE': 'dnase', 'ATAC': 'atac', 'CAGE': 'cage',
    'SPLICE_SITES': 'splice_sites', 'SPLICE_SITE_USAGE': 'splice_site_usage',
    'SPLICE_JUNCTIONS': 'splice_junctions', 'CONTACT_MAPS': 'contact_maps', 'PROCAP': 'procap',
}
requested_outputs = {getattr(dna_client.OutputType, a) for a in args.assays}

def rasterize_fills():
    """Rasterize matplotlib fill collections so dense AG signal panels render
    cleanly in vector PDF / Illustrator (per helper_functions.plot_baseline)."""
    for ax in plt.gcf().get_axes():
        for collection in ax.collections:
            collection.set_rasterized(True)


_ACGT = set('ACGTN')  # N tolerated (unknown base; common in long-read assemblies)

for _, row in variant_df.iterrows():
    vid = row['ID']
    chrom, pos = row['CHROM'], int(row['POS'])
    var_ref, var_alt = row['REF'], row['ALT']
    # Reject anything that isn't a literal ACGT(N) sequence: NaN, empty, missing,
    # VCF symbolic alleles like <INS> / <DEL>, multi-allelic <ALT1,ALT2>, etc.
    # AG's score/predict APIs require a real sequence, not a symbolic placeholder.
    if not isinstance(var_alt, str) or not isinstance(var_ref, str):
        print(f'  skip {vid}: ALT/REF is not a string'); continue
    if var_alt in ('.', '', '*') or var_alt.startswith('<') or ',' in var_alt:
        print(f'  skip {vid}: ALT is symbolic / missing / multi-allelic ({var_alt!r})'); continue
    if not set(var_alt.upper()) <= _ACGT or not set(var_ref.upper()) <= _ACGT:
        print(f'  skip {vid}: ALT or REF contains non-ACGT bases'); continue
    # Normalise to uppercase once after validation. AG accepts either case,
    # but uppercase is the canonical genomic-sequence convention and keeps
    # downstream length / variant-construction calls consistent.
    var_ref = var_ref.upper()
    var_alt = var_alt.upper()
    insert_len = len(var_alt) - len(var_ref)
    if insert_len <= 0:
        print(f'  skip {vid}: ALT not longer than REF (not an insertion)'); continue
    # Per-variant visualization flanks. Two layers of overrides on the CLI
    # --flank-bp default:
    #   1. `FLANK_BP` column sets a symmetric flank for that row.
    #   2. `LEFT_FLANK_BP` / `RIGHT_FLANK_BP` columns enable asymmetric flanking
    #      when present (e.g. to crop the view so the variant + downstream
    #      target gene get more visual real estate without showing irrelevant
    #      upstream context). If only one of the two is present, the other
    #      falls back to FLANK_BP (or args.flank_bp).
    flank_bp = args.flank_bp
    if 'FLANK_BP' in row.index and pd.notna(row['FLANK_BP']) and str(row['FLANK_BP']).strip():
        flank_bp = int(row['FLANK_BP'])
    flank_left = flank_bp
    flank_right = flank_bp
    if 'LEFT_FLANK_BP' in row.index and pd.notna(row['LEFT_FLANK_BP']) and str(row['LEFT_FLANK_BP']).strip():
        flank_left = int(row['LEFT_FLANK_BP'])
    if 'RIGHT_FLANK_BP' in row.index and pd.notna(row['RIGHT_FLANK_BP']) and str(row['RIGHT_FLANK_BP']).strip():
        flank_right = int(row['RIGHT_FLANK_BP'])
    asymmetric = (flank_left != flank_right)
    flank_label = f'L={flank_left:,}/R={flank_right:,}' if asymmetric else f'{flank_bp:,}'
    print(f'\n{vid}: {chrom}:{pos} INS {insert_len} bp ({args.label}; flank={flank_label} bp)')

    # Build INSERTION variant
    variant = genome.Variant(
        chromosome=chrom, position=pos,
        reference_bases=var_ref, alternate_bases=var_alt,
        name=f'{vid}_Insertion',
    )
    # PREDICTION interval: ALWAYS symmetric around the insertion site, so when
    # AG resizes to 1 Mb it centers exactly on `pos` and the variant lies at
    # row N // 2 of the values arrays. The asymmetric LEFT/RIGHT flanks (and
    # the VIEW_START/VIEW_END override below) only affect VISUALIZATION, not
    # the prediction call (AG always returns 1 Mb centered on POS regardless
    # of the pre-resize interval length).
    predict_interval = genome.Interval(
        chromosome=chrom,
        start=max(0, pos - flank_bp),
        end=pos + flank_bp,
        strand='+',
    )

    # Optional VIEW_START / VIEW_END columns override the flank-based logic
    # entirely. Useful for views that aren't variant-centered (e.g. zooming in
    # on the target gene downstream of the insertion). Both must be present
    # together; values are absolute genome coordinates.
    use_view = False
    if ('VIEW_START' in row.index and 'VIEW_END' in row.index
        and pd.notna(row['VIEW_START']) and pd.notna(row['VIEW_END'])
        and str(row['VIEW_START']).strip() and str(row['VIEW_END']).strip()):
        use_view = True
        view_start = int(row['VIEW_START'])
        view_end = int(row['VIEW_END'])

    # VISUALIZATION interval: asymmetric flanks supported (left and right
    # independent). In alt frame, extend the right edge by insert_len so the
    # inserted sequence's own predicted signal fits in the displayed range.
    # If VIEW_START/VIEW_END are set, use them directly (no flank_extension
    # added — user-specified end is treated as a hard bound).
    flank_extension = insert_len if args.frame == 'alt' else 0
    if use_view:
        interval = genome.Interval(chromosome=chrom, start=view_start,
                                   end=view_end, strand='+')
    else:
        interval = genome.Interval(
            chromosome=chrom,
            start=max(0, pos - flank_left),
            end=pos + flank_right + flank_extension,
            strand='+',
        )

    out = dna_model.predict_variant(
        interval=predict_interval.resize(2**20),  # AG receptive field 1 Mb, centered on pos
        variant=variant,
        requested_outputs=requested_outputs,
        ontology_terms=[args.ontology],
    )

    # Per-assay (alt, ref) pairs.
    # AG centers its 1-Mb prediction window on the variant. For an INS of length L,
    # alt's UPSTREAM-of-insertion rows align with ref correctly; the next L rows
    # are AG's prediction *for the inserted sequence itself*; rows after that are
    # alt's downstream-of-insertion (which in ref coords are also aligned but
    # plotted L positions too far downstream).
    #
    # Two coordinate frames available:
    #   --frame alt (default): extended x-axis = [start, end + L bp). Both alt
    #     and ref get padded to length N + insert_rows. Alt keeps its insert
    #     signal at rows [k, k+L); ref gets L zero rows spliced in at row k
    #     (ref has no equivalent at the insertion site). Downstream of the
    #     insertion, both are correctly aligned. Variant rendered as a yellow
    #     highlight band [pos, pos+L) instead of a thin line.
    #   --frame ref: drop alt's insert rows [k, k+L), shift downstream up by L
    #     (insert signal lost; tighter ref-frame view).
    # Two parallel pair lists:
    #   pairs_orig — untouched alt/ref straight from AG. Used for the BASELINE
    #     plot (ref-only, "what's there before the insertion") and for ref-frame.
    #   pairs_var  — frame-transformed pair, used for the DIFF and OVERLAY plots
    #     where the variant effect is being visualized.
    pairs_orig = []
    pairs_var  = []
    INTERVAL_BP = 1 << 20  # AG receptive field after resize(2**20)
    for assay in args.assays:
        attr = assay_attr_map.get(assay)
        if attr is None: continue
        if not (hasattr(out.alternate, attr) and hasattr(out.reference, attr)):
            print(f'  skip {assay}: not in AG output'); continue
        # No try/except wrapping the filter + splice block: the hasattr check
        # above already handles "assay not returned by AG". Anything that throws
        # below this point (shape mismatch, metadata-column rename, etc.) is a
        # real bug that should surface, not get silently downgraded to "skip".
        alt_td = getattr(out.alternate, attr)
        ref_td = getattr(out.reference, attr)
        # Filter ChIP_HISTONE to a user-specified mark subset (substring match
        # against the track 'name' column in metadata). sys.exit (not raise)
        # on zero-match — see comment inline below for why.
        if assay == 'CHIP_HISTONE' and args.chip_marks:
            orig_n = alt_td.num_tracks
            track_names = alt_td.metadata['name'].astype(str)
            keep_idx = [i for i, n in enumerate(track_names)
                        if any(m in n for m in args.chip_marks)]
            if not keep_idx:
                # User-input mismatches are configuration errors and must
                # terminate the whole script with a clear listing of what
                # tracks WERE available.
                sys.exit(
                    f'ERROR: --chip-marks={args.chip_marks} matched zero of {orig_n} '
                    f'CHIP_HISTONE tracks for biosample. Available track names:\n  '
                    + '\n  '.join(sorted(set(track_names)))
                )
            alt_td = alt_td.select_tracks_by_index(keep_idx)
            ref_td = ref_td.select_tracks_by_index(keep_idx)
            print(f'  {assay}: kept {len(keep_idx)} of {orig_n} tracks ({args.chip_marks})')
        # Filter RNA_SEQ to a user-specified flavor subset (e.g. "polyA plus"
        # vs "total RNA-seq"). Same substring-match logic as --chip-marks.
        if assay == 'RNA_SEQ' and args.rna_track_filter:
            orig_n = alt_td.num_tracks
            track_names = alt_td.metadata['name'].astype(str)
            keep_idx = [i for i, n in enumerate(track_names)
                        if any(m in n for m in args.rna_track_filter)]
            if not keep_idx:
                sys.exit(
                    f'ERROR: --rna-track-filter={args.rna_track_filter} matched zero of {orig_n} '
                    f'RNA_SEQ tracks for biosample. Available track names:\n  '
                    + '\n  '.join(sorted(set(track_names)))
                )
            alt_td = alt_td.select_tracks_by_index(keep_idx)
            ref_td = ref_td.select_tracks_by_index(keep_idx)
            print(f'  {assay}: kept {len(keep_idx)} of {orig_n} tracks ({args.rna_track_filter})')

        # Drop unstranded RNA-seq tracks by default (Liu et al 2026 convention:
        # use only stranded polyA + stranded total tracks). If GENE_STRAND
        # column is set in the tab, further restrict to that strand only.
        if assay == 'RNA_SEQ':
            orig_n = alt_td.num_tracks
            strands = alt_td.metadata['strand'].astype(str)
            target_strands = ['+', '-']
            if 'GENE_STRAND' in row.index and pd.notna(row['GENE_STRAND']):
                gs = str(row['GENE_STRAND']).strip()
                if gs in ('+', '-'):
                    target_strands = [gs]
            keep_idx = [i for i, s in enumerate(strands) if s in target_strands]
            if not keep_idx:
                sys.exit(
                    f'ERROR: no RNA_SEQ tracks matched strand filter '
                    f'(target_strands={target_strands}). Available strands: '
                    f'{sorted(set(strands))}'
                )
            alt_td = alt_td.select_tracks_by_index(keep_idx)
            ref_td = ref_td.select_tracks_by_index(keep_idx)
            print(f'  {assay}: kept {len(keep_idx)} of {orig_n} stranded tracks '
                  f'(strands={target_strands})')
        pairs_orig.append((assay, alt_td, ref_td))

        if insert_len > 0:
            vals_a = alt_td.values
            vals_r = ref_td.values
            N = vals_a.shape[0]
            bp_per_row = INTERVAL_BP / N
            insert_rows = max(1, int(round(insert_len / bp_per_row)))
            # Derive the variant row from the returned TrackData interval and
            # assert it matches the expected center row (N//2 ± 1). Avoids the
            # hidden N//2 assumption — if a future AG API change shifted the
            # window centering or returned a cropped interval, the splice
            # would silently miscenter the yellow band; the assertion catches
            # that case loudly. Empirically matches N//2 exactly across
            # resolutions (RNA 1 bp, ChIP 128 bp, ATAC 1 bp) for AG SDK v0.6.1;
            # the ±1-row tolerance absorbs legitimate integer-rounding /
            # bin-boundary edge cases without crashing otherwise-valid runs.
            k = round((pos - alt_td.interval.start) / alt_td.resolution)
            assert abs(k - N // 2) <= 1, (
                f'AG window not centered on variant for {assay}: '
                f'k_derived={k} vs N//2={N//2} (drift>1 row) '
                f'(interval=[{alt_td.interval.start}, {alt_td.interval.end}], '
                f'pos={pos}, resolution={alt_td.resolution})'
            )

            if args.frame == 'alt':
                alt_ext = np.vstack([vals_a, np.zeros((insert_rows, vals_a.shape[1]))])
                ref_ext = np.vstack([
                    vals_r[:k, :],
                    np.zeros((insert_rows, vals_r.shape[1])),
                    vals_r[k:, :],
                ])
                extended_interval = genome.Interval(
                    chromosome=alt_td.interval.chromosome,
                    start=alt_td.interval.start,
                    end=alt_td.interval.end + insert_rows * alt_td.resolution,
                    strand=alt_td.interval.strand,
                )
                alt_td_var = track_data.TrackData(
                    alt_ext, alt_td.metadata,
                    resolution=alt_td.resolution, interval=extended_interval,
                )
                ref_td_var = track_data.TrackData(
                    ref_ext, ref_td.metadata,
                    resolution=ref_td.resolution, interval=extended_interval,
                )
            else:  # frame == 'ref'
                spliced = np.vstack([
                    vals_a[:k, :],
                    np.pad(vals_a[k + insert_rows:, :], ((0, insert_rows), (0, 0))),
                ])
                alt_td_var = track_data.TrackData(
                    spliced, alt_td.metadata,
                    resolution=alt_td.resolution, interval=alt_td.interval,
                )
                ref_td_var = ref_td
            pairs_var.append((assay, alt_td_var, ref_td_var))
        else:
            pairs_var.append((assay, alt_td, ref_td))

    # Output filename: include optional FAMILY and GENE columns from the variant
    # tab if present, so files are self-describing (e.g.
    # `LTR5_Hs_HLA-DQA2_SvimAsm00060017_GM12878_overlay.pdf`). Falls back to
    # `{vid}_{label}_{plot}.pdf` when those columns aren't in the tab.
    prefix_parts = []
    for col in ('FAMILY', 'GENE'):
        if col in row.index and pd.notna(row[col]) and str(row[col]).strip():
            prefix_parts.append(str(row[col]).strip())
    base = '_'.join(prefix_parts + [vid, args.label]) if prefix_parts else f'{vid}_{args.label}'

    # Original (un-extended) interval — for the baseline plot. Honours the
    # asymmetric LEFT/RIGHT flanks (same x-axis range as the alt-frame plots
    # minus the insert-length extension), or the VIEW_START/VIEW_END override.
    if use_view:
        interval_orig = genome.Interval(chromosome=chrom, start=view_start,
                                        end=view_end, strand='+')
    else:
        interval_orig = genome.Interval(
            chromosome=chrom,
            start=max(0, pos - flank_left),
            end=pos + flank_right,
            strand='+',
        )

    def render_and_save(tracks, suffix, title_extra, plot_interval, draw_band):
        plot_components.plot(
            tracks,
            annotations=([plot_components.VariantAnnotation([variant], alpha=0.8)]
                         if not draw_band else []),
            interval=plot_interval,
            title=f'{vid} ({chrom}:{pos}, {insert_len} bp INS) in {args.label} — {title_extra}',
        )
        # Highlight the inserted sequence as a yellow band on signal axes when
        # in alt frame (skip the transcript-annotation axis at top). The band
        # marks the TRUE genome-coordinate insertion region [pos, pos+insert_len)
        # plus a half-ChIP-bin pad on each side to absorb the line-interpolation
        # wedge artifact at band edges (line plots connecting row k-1 with
        # diff≈0 to row k with diff=S_sva create a triangular fill straddling
        # the band edge, ~half a bin wide for ChIP tracks at 128 bp resolution;
        # RNA tracks at 1 bp resolution barely show the artifact).
        #
        # Bin-quantization caveat (multi-assay figures): at coarse resolutions
        # (ChIP 128 bp), the rendered data extends only to pos + insert_rows*128,
        # which can be slightly less than pos + insert_len (e.g. 256 vs 312 for
        # a 312 bp insert on a 128 bp bin grid). The band marks the true
        # insertion in genome coordinates, not the bin-quantized data extent —
        # the small right-edge gap on ChIP panels is honest bin sampling, not
        # a misalignment. RNA / ATAC panels (1 bp resolution) match the band
        # exactly.
        if draw_band:
            BAND_VISUAL_PAD_BP = 64  # = 128 bp / 2 (half a ChIP bin width)
            # Skip the transcript-annotation axis (the first one) when present,
            # so the yellow variant band only paints onto signal panels.
            all_axes = plt.gcf().get_axes()
            sig_axes = all_axes[1:] if transcript_extractor is not None else all_axes
            for ax in sig_axes:
                ax.axvspan(pos - BAND_VISUAL_PAD_BP,
                           pos + insert_len + BAND_VISUAL_PAD_BP,
                           color='gold', alpha=0.25, zorder=0, lw=0)
        rasterize_fills()
        plt.gcf().set_size_inches(18, 14)
        path = os.path.join(out_dir, f'{base}_{suffix}.pdf')
        plt.savefig(path, dpi=200, bbox_inches='tight')
        plt.close()
        print(f'  saved {path}')

    def transcript_panel_for(iv):
        return ([plot_components.TranscriptAnnotation(transcript_extractor.extract(iv))]
                if transcript_extractor is not None else [])

    def _ylabel_for(assay, td):
        """Show strand info for RNA-seq / CAGE / PROCAP tracks (where +/- vs .
        carries biological meaning); plain biosample+name for chromatin tracks."""
        has_meta = hasattr(td, 'metadata') and 'biosample_name' in td.metadata.columns
        if assay in ('RNA_SEQ', 'CAGE', 'PROCAP'):
            return '{biosample_name}\n{name}\nstrand: {strand}' if has_meta else '{name}\nstrand: {strand}'
        return '{biosample_name}\n{name}' if has_meta else '{name}'

    # 1) BASELINE — reference signal at original interval, no extension, no band
    if 'baseline' in args.plots:
        baseline_tracks = list(transcript_panel_for(interval_orig))
        for assay, _alt, ref_td in pairs_orig:
            baseline_tracks.append(plot_components.Tracks(
                tdata=ref_td, ylabel_template=_ylabel_for(assay, ref_td), filled=True
            ))
        render_and_save(baseline_tracks, 'baseline', 'baseline / reference (pre-insertion)',
                        plot_interval=interval_orig, draw_band=False)

    # 2) DIFFERENCE — alt − ref in the chosen frame
    if 'diff' in args.plots:
        diff_tracks = list(transcript_panel_for(interval))
        for assay, alt_td, ref_td in pairs_var:
            diff = alt_td - ref_td
            diff_tracks.append(plot_components.Tracks(
                tdata=diff, ylabel_template=_ylabel_for(assay, diff), filled=True
            ))
        render_and_save(diff_tracks, 'diff', 'difference (alt − ref; up=gain, down=loss)',
                        plot_interval=interval, draw_band=(args.frame == 'alt'))

    # 3) OVERLAY — gray ref + red alt
    if 'overlay' in args.plots:
        overlay_tracks = list(transcript_panel_for(interval))
        for assay, alt_td, ref_td in pairs_var:
            overlay_tracks.append(plot_components.OverlaidTracks(
                {'Reference': ref_td, 'Insertion': alt_td},
                colors={'Reference': 'dimgrey', 'Insertion': 'red'},
            ))
        render_and_save(overlay_tracks, 'overlay', 'reference (gray) vs insertion (red)',
                        plot_interval=interval, draw_band=(args.frame == 'alt'))

print(f'\nDone — figures in {out_dir}/')
