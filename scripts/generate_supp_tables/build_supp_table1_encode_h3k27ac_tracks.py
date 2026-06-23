"""
Build Supp Table 1 — ENCODE H3K27ac peak files matched to AlphaGenome training tracks.

Joins the Avsec et al. 2026 Suppl Table 2 (full track metadata, sheet 2 of
the MOESM3 xlsx) with our resolved ENCODE narrowPeak manifest (one row per
AG-training bigWig with peak_ENCFF / pipeline_version / md5sum).

Output:
    supptables/supp_table_1_encode_h3k27ac_tracks.tsv

One row per AG training biosample (159 human H3K27ac × encode tracks).
Multi-experiment biosamples list their ENCSR / bigWig ENCFF / peak ENCFF /
pipeline_version / md5sum fields as comma-separated parallel arrays
(aligned to ordering of Experiment accession).

Usage:
    python scripts/generate_supp_tables/build_supp_table1_encode_h3k27ac_tracks.py
"""

import csv
import openpyxl
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
AVSEC_XLSX = ROOT / 'Avsec_Nature2026' / '41586_2025_10014_MOESM3_ESM.xlsx'
RES_TSV = ROOT / 'data' / 'encode_h3k27ac' / 'peak_resolution_human.tsv'
OUT_TSV = ROOT / 'supptables' / 'supp_table_1_encode_h3k27ac_tracks.tsv'

# ENCODE files are served via a stable URL pattern derived from the accession.
# We expose both the bigWig and narrowPeak URLs in the supp table so a reader
# can click straight through without needing to know the @@download convention.
ENCODE_URL_FMT = 'https://www.encodeproject.org/files/{accession}/@@download/{accession}.{ext}'

def make_urls(accession_csv: str, ext: str) -> str:
    """Comma-joined accessions → comma-joined ENCODE download URLs."""
    accs = [a.strip() for a in (accession_csv or '').split(',') if a.strip()]
    return ','.join(ENCODE_URL_FMT.format(accession=a, ext=ext) for a in accs)

# ── Load resolution TSV: (ENCSR, bigwig_ENCFF) -> per-peak metadata ────────
res = {}
with open(RES_TSV) as f:
    for r in csv.DictReader(f, delimiter='\t'):
        res[(r['ENCSR'], r['AG_bigwig_ENCFF'])] = {
            'peak_ENCFF':       r.get('peak_ENCFF', ''),
            'peak_output_type': r.get('peak_output_type', ''),
            'ENCAN':            r.get('ENCAN', ''),
            'pipeline_version': r.get('pipeline_version', ''),
            'md5sum':           r.get('md5sum', ''),
            'preferred_default': r.get('preferred_default', ''),
            'resolution_note':  r.get('note', ''),
        }

# ── Read Avsec Supp 2 (full track metadata) ────────────────────────────────
wb = openpyxl.load_workbook(AVSEC_XLSX, read_only=True, data_only=True)
sheet = wb['Suppl Table 2 Track metadata (f']
iters = sheet.iter_rows(values_only=True)
header = next(iters)

# Avsec Supp 2 columns to carry forward verbatim (same order as the
# published Avsec xlsx). The 4dn_* columns (enzyme/dataset/condition) are
# dropped — they're 4DN-specific and always blank for ENCODE H3K27ac.
CARRY = [
    'organism', 'output_type', 'name', 'strand', 'track_index',
    'Assay title', 'File assembly', 'data_source', 'Target label',
    'ontology_curie', 'biosample_name', 'biosample_type',
    'gtex_tissue', 'gtex_tissue_group',
    'is_treated', 'is_fractional', 'genetically_modified', 'genetic_modification',
    'ontology_type', 'audit_flag_priority', 'endedness', 'biosample_life_stage',
    'Experiment accession', 'File accession',
    'frip', 'nonzero_mean',
]

# Output columns: 29 Avsec verbatim + 3 ENCODE peak-resolution columns
# (matching the names used in our existing SuppTables.xlsx Supp Table 1
# tab) + 4 internal-traceability columns we add for reproducibility +
# 2 self-describing ENCODE download-URL columns (`bigwig_URL` next to the
# bigWig accession `File accession`, `peak_URL` next to `peak_ENCFF`).
OUT_COLS = (
    list(CARRY[:CARRY.index('File accession') + 1])
    + ['bigwig_URL']
    + list(CARRY[CARRY.index('File accession') + 1:])
    + [
        'peak_ENCFF', 'peak_URL',
        'ENCODE4 GRCh38 pipeline_version',
        'md5sum',
        'ENCAN', 'peak_output_type', 'preferred_default', 'resolution_note',
    ]
)

rows_out = []
n_h3k27ac_human = 0
for r in iters:
    d = dict(zip(header, r))
    # Filter to H3K27ac × encode × human (matches the manuscript Methods scope)
    if d.get('Target label') != 'H3K27ac': continue
    if d.get('data_source') != 'encode':  continue
    if d.get('organism')    != 'human':   continue
    n_h3k27ac_human += 1

    # Avsec's File accession field is sometimes "ENCFF000XXX," (trailing comma)
    encsrs  = [x.strip() for x in str(d.get('Experiment accession') or '').split(',') if x.strip()]
    bigwigs = [x.strip() for x in str(d.get('File accession')       or '').split(',') if x.strip()]

    # Build the per-peak fields by looking up each (ENCSR, bigwig).
    # Internal keys (`pipeline_version`) come from our resolution TSV; the
    # output column rename to `ENCODE4 GRCh38 pipeline_version` happens
    # below to match the published Supp Table 1 schema.
    per_keys = ('peak_ENCFF', 'peak_output_type', 'ENCAN',
                'pipeline_version', 'md5sum', 'preferred_default', 'resolution_note')
    per = {k: [] for k in per_keys}
    for encsr, bw in zip(encsrs, bigwigs):
        hit = res.get((encsr, bw), {})
        for k in per_keys:
            per[k].append(hit.get(k, ''))

    out = {c: ('' if d.get(c) is None else str(d.get(c)).rstrip(',')) for c in CARRY}
    # Re-canonicalize the Experiment / File accession fields (strip trailing commas)
    out['Experiment accession'] = ','.join(encsrs)
    out['File accession']       = ','.join(bigwigs)
    out['bigwig_URL']           = make_urls(out['File accession'], 'bigWig')
    for k, vals in per.items():
        col = 'ENCODE4 GRCh38 pipeline_version' if k == 'pipeline_version' else k
        out[col] = ','.join(vals)
    out['peak_URL'] = make_urls(out.get('peak_ENCFF', ''), 'bed.gz')
    rows_out.append(out)

# ── Write TSV ──────────────────────────────────────────────────────────────
OUT_TSV.parent.mkdir(parents=True, exist_ok=True)
with open(OUT_TSV, 'w', newline='') as f:
    w = csv.DictWriter(f, fieldnames=OUT_COLS, delimiter='\t', extrasaction='ignore')
    w.writeheader()
    for row in rows_out:
        w.writerow(row)

print(f"Avsec Supp 2 H3K27ac × encode × human rows: {n_h3k27ac_human}")
print(f"Rows written: {len(rows_out)}")
print(f"Output: {OUT_TSV.relative_to(ROOT)}")

# ── Sanity: per-row peak join completeness ─────────────────────────────────
missing_peak = sum(1 for r in rows_out if not r['peak_ENCFF'])
n_audit = {}
for r in rows_out:
    af = r.get('audit_flag_priority', '') or '(unknown)'
    n_audit[af] = n_audit.get(af, 0) + 1
print(f"\nRows missing peak_ENCFF (any element): {missing_peak}")
print(f"audit_flag_priority breakdown:")
for k, v in sorted(n_audit.items(), key=lambda kv: -kv[1]):
    print(f"  {k:<16} {v}")
