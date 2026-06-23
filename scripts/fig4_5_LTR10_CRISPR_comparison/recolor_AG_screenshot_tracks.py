"""
Recolour the AlphaGenome-SDK genome-browser screenshot tracks in-place by
editing each track image's indexed palette.

The AG SDK renders each track as a 2-colour indexed-palette image (white +
one fill colour). Each image lives as a separate XObject in the page
resources. Changing each image's palette swaps the track colour without
touching geometry, axes, labels, or anything else.

Two PDF formats are supported:

  Old (matplotlib ≤ ~3.8): XObjects named /Im0../Im6, palette stored as a
  shared indirect /Stream (Im0+Im1+Im3 may all reference one stream).
  Mapping for the 7-track HCT116-only browser shot used in the CSHL poster.

  New (matplotlib ≥ ~3.9): XObjects named /I1../I9, palette stored as an
  inline pikepdf.String per image (no sharing). Mapping for the 9-track
  HCT116 + transverse-colon + sigmoid-colon browser shots used in Fig 5.

Usage:
    python scripts/fig4_5_LTR10_CRISPR_comparison/recolor_AG_screenshot_tracks.py \\
        --input  figures/FIG5_FINAL/browser_shots/raw/LTR10.ATG12_browser.pdf \\
        --output figures/FIG5_FINAL/browser_shots/LTR10.ATG12_browser.pdf
"""
import argparse
from pathlib import Path
import pikepdf


# Locked assay palette (cross-figure). Hex values match
# reference_assay_color_palette.md.
NAVY      = (0x16, 0x33, 0x59)  # TF FOSL1
STEELBLUE = (0x5B, 0x7F, 0xB5)  # TF JUND
DARKRED   = (0x8B, 0x1A, 0x1F)  # H3K27ac
ORANGE    = (0xF0, 0x8C, 0x2D)  # H3K4me1
DARKGREEN = (0x1A, 0x5F, 0x2A)  # ATAC
PURPLE    = (0x7A, 0x33, 0x70)  # RNA (both strands)

# 9-track stack (Fig 4/Fig 5). All HCT116 tracks first; "Normal colon
# (sigmoid)" group at the bottom. Locked palette preserved; sigmoid
# tracks reuse the same assay colour as their HCT116 counterpart so the
# cancer-vs-healthy contrast reads at a glance.
TRACK_COLOURS_NEW = {
    'I1': STEELBLUE,  # HCT116 JUND
    'I2': NAVY,       # HCT116 FOSL1
    'I3': DARKRED,    # HCT116 H3K27ac
    'I4': ORANGE,     # HCT116 H3K4me1
    'I5': DARKGREEN,  # HCT116 ATAC
    'I6': PURPLE,     # HCT116 total RNA, + strand
    'I7': PURPLE,     # HCT116 total RNA, - strand
    'I8': DARKRED,    # sigmoid colon H3K27ac
    'I9': ORANGE,     # sigmoid colon H3K4me1
}

# 7-track stack (legacy — CSHL poster ATG12 panel, HCT116-only).
TRACK_COLOURS_OLD = {
    'Im5': NAVY,       # FOSL1
    'Im6': STEELBLUE,  # JUND
    'Im0': DARKRED,    # H3K27ac
    'Im2': ORANGE,     # H3K4me1
    'Im1': DARKGREEN,  # ATAC
    'Im3': PURPLE,     # RNA+
    'Im4': PURPLE,     # RNA-
}


def recolour_indexed_image(pdf, img_obj, new_rgb):
    """Replace the second palette entry (after white) with new_rgb.

    Handles both palette storage formats:
    - Indirect /Stream (old PDFs): may be shared across multiple images,
      so we clone into a fresh indirect stream before rewriting.
    - Inline pikepdf.String (new PDFs): no sharing — rewrite in place.
    """
    cs = img_obj['/ColorSpace']
    if not isinstance(cs, pikepdf.Array) or len(cs) < 4 or str(cs[0]) != '/Indexed':
        raise ValueError(f'Not an Indexed colour-space image: {cs}')
    palette = cs[3]
    if isinstance(palette, pikepdf.String):
        raw = bytes(palette)
        if len(raw) != 6:
            raise ValueError(f'Expected 2-entry palette (6 bytes), got {len(raw)}: {raw!r}')
        cs[3] = pikepdf.String(raw[:3] + bytes(new_rgb))
    else:
        raw = bytes(palette.read_bytes())
        if len(raw) != 6:
            raise ValueError(f'Expected 2-entry palette (6 bytes), got {len(raw)}: {raw!r}')
        fresh = pdf.make_stream(raw[:3] + bytes(new_rgb))
        img_obj['/ColorSpace'] = pikepdf.Array([cs[0], cs[1], cs[2], fresh])


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--input',  required=True)
    ap.add_argument('--output', required=True)
    args = ap.parse_args()

    pdf = pikepdf.open(args.input)
    page = pdf.pages[0]
    xobjs = page.Resources.XObject

    # Pick the right map by detecting which keys are present.
    if any(f'/I{i}' in xobjs for i in range(1, 10)):
        track_map, fmt = TRACK_COLOURS_NEW, '8-track (new format)'
    else:
        track_map, fmt = TRACK_COLOURS_OLD, '7-track (old format)'
    print(f'Detected layout: {fmt}')

    for img_name, rgb in track_map.items():
        key = '/' + img_name
        if key not in xobjs:
            print(f'  skip {img_name}: not on page')
            continue
        recolour_indexed_image(pdf, xobjs[key], rgb)
        print(f'  {img_name} -> rgb({rgb[0]},{rgb[1]},{rgb[2]}) = #{rgb[0]:02X}{rgb[1]:02X}{rgb[2]:02X}')

    # Strip Illustrator's private artwork cache (PieceInfo et al). If the
    # PDF was saved from Illustrator, it embeds the editable AI-native
    # artwork as private data inside /PieceInfo. Illustrator reads that
    # cache and ignores our public PDF colour edits — Preview shows the
    # new colours, Illustrator shows the old. Removing the cache forces
    # Illustrator to reopen the PDF using the public page content.
    illustrator_cache_keys = ['/PieceInfo', '/LastModified', '/Thumb']
    for key in illustrator_cache_keys:
        if key in page:
            del page[key]

    Path(args.output).parent.mkdir(parents=True, exist_ok=True)
    pdf.save(args.output)
    print(f'\nWrote {args.output}')


if __name__ == '__main__':
    main()
