"""
Build Supp Table 10 — AlphaGenome predictions for the six CRISPRi-validated LTR10
enhancers (plus one matched control deletion) joined with the Ivancevic 2024
CRISPRi DESeq2 results (Tables S15-S20). One row per (variant, gene).

Three categories of rows per variant (sorted within each variant by
|crispr_log2FoldChange| descending; AG-only rows with no CRISPR value are
appended at the bottom of each block sorted by |AG_RNA_raw_score| descending):

  1. Both:        gene scored by AG (within ±500 kb of variant, protein_coding)
                  AND listed in the CRISPRi DESeq2 table. AG + CRISPR cols
                  populated; notes blank.
  2. CRISPR-only: protein_coding gene with padj < 0.05 AND
                  log2FoldChange < 0 (down-regulated) in the DESeq2 table,
                  same chromosome as the variant, within ±10 Mb of the variant
                  midpoint, but outside AG's ±500 kb scoring window. AG cols
                  blank; notes flag the architectural exclusion. Up-regulated
                  CRISPR-only hits are excluded — they are typically secondary
                  / trans effects and don't represent the canonical "enhancer
                  positively regulates target gene" relationship AG should
                  have caught.
  3. AG-only:     gene scored by AG but absent from the DESeq2 table (likely
                  not expressed in HCT116 / filtered out by DESeq2 independent
                  filtering). CRISPR cols blank; notes flag the absence.

The ±10 Mb cis scope matches the coord-scatter plot window
(`scripts/fig4_5_LTR10_CRISPR_comparison/plot_variant_coord_scatter.py --window-mb 10`), so every CRISPR-sig
"escaped" gene that's relevant to Fig 4 supp panels has a row here.

Scope is restricted to protein_coding genes throughout, with an explicit
allowlist override for legacy / unsymbolized loci where the protein-coding
filter would discard a real AG↔CRISPR match. Currently the allowlist contains
just ENSG00000248112 (legacy symbol AC108174.1, a lncRNA ~120 kb upstream of
LTR10.XRCC4 that AG predicts at raw_score=-1.60 and CRISPR validates at
log2FC=-3.09). The Ivancevic 2024 CRISPRi DESeq2 tables use the legacy
symbol; an alias map (AC108174.1 -> ENSG00000248112) normalises the join.

The control deletion (`control_deletion_near_ATG12`) was not CRISPR-tested in
Ivancevic 2024, so its rows are all AG-only by construction.

Inputs:
  --variants        data/LTR10_variants.tab (ID/CHROM/POS/REF/ALT/...)
  --rna-dir         results/AG_LFC_LTR10_CRISPR/  (per-variant AG RNA CSVs)
  --chromatin-dir   results/AG_chromatin_LTR10_CRISPR/  (per-variant chromatin)
  --crispr-dir      Ivancevic_SciAdv2024/CRISPRi_results/ (Tables S15-S20 xlsx)
  --gencode         data/gencode.v46.annotation.feather (TSS lookup)
  --cell-line       HCT116 (file suffix; retained as a column)

Per-row columns:
  variant_id, chrom, pos, SVLEN, type, target_gene, study_id, cell_line,
  gene_name, gene_id, gene_type, gene_strand, gene_TSS,
  crispr_log2FoldChange, crispr_padj, crispr_distance_kb, notes,
  AG_RNA_raw_score, AG_RNA_quantile_score,
  AG_ATAC_raw_score, AG_ATAC_quantile_score,
  AG_H3K27ac_raw_score, AG_H3K27ac_quantile_score,
  AG_H3K4me1_raw_score, AG_H3K4me1_quantile_score

RNA convention: mean of stranded RNA-seq tracks only (track_strand in {+,-});
unstranded polyA dropped. Matches Supp Tables 4/5/7/8.

Per-variant identity columns (variant_id, chrom, pos, SVLEN, type, target_gene,
study_id, cell_line) and the four chromatin columns hold their literal values
only on the first row of each variant block; subsequent rows leave those cells
blank.

Usage:
    python scripts/generate_supp_tables/build_supp_table10_LTR10_CRISPR_AG_predictions.py \\
        --variants      data/LTR10_variants.tab \\
        --rna-dir       results/AG_LFC_LTR10_CRISPR \\
        --chromatin-dir results/AG_chromatin_LTR10_CRISPR \\
        --crispr-dir    Ivancevic_SciAdv2024/CRISPRi_results \\
        --gencode       data/gencode.v46.annotation.feather \\
        --cell-line     HCT116 \\
        --tsv           supptables/supp_table_10_LTR10_CRISPR_AG_predictions.tsv
"""
import argparse
from pathlib import Path
import numpy as np
import pandas as pd

p = argparse.ArgumentParser()
p.add_argument('--variants',      required=True)
p.add_argument('--rna-dir',       required=True)
p.add_argument('--chromatin-dir', required=True)
p.add_argument('--crispr-dir',    required=True, help='Ivancevic 2024 CRISPRi DESeq2 xlsx directory')
p.add_argument('--gencode',       required=True, help='gencode.v46.annotation.feather (TSS lookup)')
p.add_argument('--cell-line',     default='HCT116')
p.add_argument('--cis-window-mb', type=float, default=10.0,
               help='Same-chrom cis-window (Mb) for including CRISPR-sig genes outside AG (default 10 Mb)')
p.add_argument('--crispr-padj',   type=float, default=0.05,
               help='padj cutoff for "CRISPR-significant" rows (default 0.05)')
p.add_argument('--tsv',           required=True)
args = p.parse_args()

rna_dir    = Path(args.rna_dir)
chr_dir    = Path(args.chromatin_dir)
crispr_dir = Path(args.crispr_dir)

# Variant display order (enhancers first, control last). ATG12+XRCC4 lead
# because their CRISPRi-validated targets sit furthest from the deletion.
VARIANT_ORDER = [
    'LTR10.ATG12', 'LTR10.XRCC4', 'LTR10.KDM6A',
    'LTR10.FGF2',  'LTR10.MCPH1', 'LTR10.MEF2D',
    'control_deletion_near_ATG12',
]

# Variant → (CRISPR xlsx filename, sheet name). Note S20 uses CRISPR-KO
# (physical enhancer deletion) rather than CRISPRi silencing — different
# perturbation modality, but log2FC is still ALT-REF and commensurable.
# Explicit allowlist of non-protein-coding genes whose AG and/or CRISPR rows
# should appear in supp table 10 despite the protein_coding biotype filter.
# Add an entry here only when there's a clear AG↔CRISPR concordance worth
# documenting (Atma's approval required, per 2026-05-12 decision).
INCLUDE_NON_CODING = {
    'ENSG00000248112',  # legacy symbol AC108174.1 — lncRNA ~120 kb upstream
                        # of LTR10.XRCC4; AG raw=-1.60, CRISPR log2FC=-3.09
                        # (within AG window; appears as a Both row)
    'ENSG00000270069',  # MIR222HG — lncRNA ~845 kb downstream of LTR10.KDM6A;
                        # CRISPR log2FC=-0.81, padj=0.024. Outside AG window
                        # → enters as a CRISPR-only row.
}

# Alias map: CRISPR table gene name -> canonical gene_id used by AG / GENCODE
# v46. The Ivancevic 2024 DESeq2 tables predate Ensembl's drop of legacy
# clone-based symbols, so we need an explicit rename to join AC108174.1 to
# its current Ensembl ID. Add entries here only when needed to recover a
# specific allowlisted gene's CRISPR row.
CRISPR_NAME_TO_ENSG = {
    'AC108174.1': 'ENSG00000248112',
}

# Display-name overrides applied at the very end of the build, after all
# joins. Used to surface legacy / informative gene symbols that GENCODE v46
# has dropped (Ensembl reports the bare ENSG ID instead). Reader-facing
# benefit: AC108174.1 is the symbol Ivancevic 2024 published, so the table
# should match the paper rather than the v46 internal name.
GENE_NAME_OVERRIDES = {
    'ENSG00000248112': 'AC108174.1',
}

CRISPR_FILES = {
    'LTR10.ATG12':  ('Table_S15.xlsx', 'DEseq2 LTR10.ATG12 CRISPRi'),
    'LTR10.XRCC4':  ('Table_S16.xlsx', 'DEseq2 LTR10.XRCC4 CRISPRi'),
    'LTR10.MEF2D':  ('Table_S17.xlsx', 'DEseq2 LTR10.MEF2D CRISPRi'),
    'LTR10.FGF2':   ('Table_S18.xlsx', 'DEseq2 LTR10.FGF2 CRISPRi'),
    'LTR10.MCPH1':  ('Table_S19.xlsx', 'DEseq2 LTR10.MCPH1 CRISPRi'),
    'LTR10.KDM6A':  ('Table_S20.xlsx', 'DEseq2 LTR10.KDM6A CRISPR-KO'),
    'control_deletion_near_ATG12': None,  # not CRISPR-tested
}

variants = pd.read_table(args.variants, sep='\t')
variants = variants.set_index('ID').reindex(VARIANT_ORDER).reset_index().rename(columns={'index':'ID'})
missing = variants['ID'][variants['CHROM'].isna()].tolist()
if missing:
    raise SystemExit(f"Missing rows in --variants for: {missing}")

# GENCODE TSS lookup (all biotypes — the biotype filter is applied per-row
# below, with the INCLUDE_NON_CODING allowlist as override). Loading all
# biotypes lets us look up lncRNAs like MIR222HG when they're allowlisted.
print(f"Loading GENCODE...")
gtf = pd.read_feather(args.gencode)
gtf = gtf[gtf['Feature'] == 'gene'].copy()
gtf['TSS'] = gtf.apply(lambda r: r['Start'] if r['Strand'] == '+' else r['End'], axis=1)
gtf['gene_id_base'] = gtf['gene_id'].str.split('.').str[0]
gencode_by_name = gtf.drop_duplicates('gene_name', keep='first').set_index('gene_name')

def target_gene(vid: str) -> str:
    if vid.startswith('LTR10.'): return vid.split('.', 1)[1]
    if vid == 'control_deletion_near_ATG12': return 'ATG12'
    raise ValueError(f"unrecognized variant id: {vid}")

def chromatin_per_variant(chr_csv: Path) -> dict:
    df = pd.read_csv(chr_csv)
    out = {}
    for mark, mask in [
        ('ATAC',    df['output_type'] == 'ATAC'),
        ('H3K27ac', (df['output_type'] == 'CHIP_HISTONE') &
                     df['track_name'].str.contains('H3K27ac', na=False)),
        ('H3K4me1', (df['output_type'] == 'CHIP_HISTONE') &
                     df['track_name'].str.contains('H3K4me1', na=False)),
    ]:
        sub = df[mask]
        if len(sub) == 0:
            out[f'AG_{mark}_raw_score']      = pd.NA
            out[f'AG_{mark}_quantile_score'] = pd.NA
        else:
            out[f'AG_{mark}_raw_score']      = float(sub['raw_score'].mean())
            out[f'AG_{mark}_quantile_score'] = float(sub['quantile_score'].mean())
    return out

def load_crispr(vid: str) -> pd.DataFrame | None:
    entry = CRISPR_FILES.get(vid)
    if entry is None: return None
    fn, sh = entry
    df = pd.read_excel(crispr_dir / fn, sheet_name=sh)
    # Apply legacy-symbol aliases so the merge joins legacy CRISPR names to
    # current AG / GENCODE v46 gene IDs (e.g. AC108174.1 -> ENSG00000248112).
    df = df.copy()
    df['gene'] = df['gene'].replace(CRISPR_NAME_TO_ENSG)
    return df[['gene','log2FoldChange','padj']].rename(
        columns={'gene':'gene_name','log2FoldChange':'crispr_log2FoldChange','padj':'crispr_padj'})

CIS_WIN_BP = int(args.cis_window_mb * 1_000_000)
BLANK = ''

all_rows = []
for _, v in variants.iterrows():
    vid   = v['ID']
    chrom = v['CHROM']
    pos   = int(v['POS'])
    alt   = v['ALT'] if isinstance(v['ALT'], str) else ''
    svlen = len(alt.replace('.', '')) - len(v['REF'])
    var_mid = pos + max(abs(svlen), 1) // 2  # signed-SVLEN: mid of REF span
    vtype = 'control' if vid.startswith('control_') else 'enhancer'
    tgt = target_gene(vid)
    identity = {
        'variant_id': vid, 'chrom': chrom, 'pos': pos, 'SVLEN': svlen,
        'type': vtype, 'target_gene': tgt,
        'study_id': v.get('Study_ID', ''), 'cell_line': args.cell_line,
    }
    chr_cols = chromatin_per_variant(chr_dir / f'{vid}_{args.cell_line}.csv')

    # --- (1) AG protein-coding genes ---
    rna = pd.read_csv(rna_dir / f'{vid}_{args.cell_line}.csv')
    rna = rna[rna['track_strand'].isin(['+', '-'])]
    # Biotype filter with INCLUDE_NON_CODING allowlist override: keep
    # protein_coding rows, plus any row whose gene_id is explicitly approved
    # for inclusion (e.g. AC108174.1 / ENSG00000248112).
    rna = rna[(rna['gene_type'] == 'protein_coding') | (rna['gene_id'].isin(INCLUDE_NON_CODING))]
    ag = (rna.groupby(['gene_id','gene_name','gene_type','gene_strand'], as_index=False)
            .agg(AG_RNA_raw_score=('raw_score','mean'),
                 AG_RNA_quantile_score=('quantile_score','mean')))

    # --- (2) CRISPR DESeq2 table ---
    crispr = load_crispr(vid)

    # --- merge AG with CRISPR by gene_name ---
    if crispr is not None:
        merged = ag.merge(crispr, on='gene_name', how='outer', indicator=True)
    else:
        merged = ag.copy()
        merged['crispr_log2FoldChange'] = pd.NA
        merged['crispr_padj'] = pd.NA
        merged['_merge'] = 'left_only'

    # --- Restrict CRISPR-only outer rows to (a) CRISPR-sig, (b) same chrom,
    #     (c) within cis-window of variant midpoint, (d) protein_coding OR
    #     explicitly allowlisted via INCLUDE_NON_CODING (e.g. MIR222HG) ---
    def gene_tss(name: str):
        return gencode_by_name['TSS'].get(name) if name in gencode_by_name.index else None
    def gene_chrom(name: str):
        return gencode_by_name['Chromosome'].get(name) if name in gencode_by_name.index else None
    def gene_strand(name: str):
        return gencode_by_name['Strand'].get(name) if name in gencode_by_name.index else None
    def gene_id_base(name: str):
        gid = gencode_by_name['gene_id'].get(name) if name in gencode_by_name.index else None
        if gid is None: return None
        return str(gid).split('.')[0]  # strip GENCODE version suffix
    def gene_biotype(name: str):
        return gencode_by_name['gene_type'].get(name) if name in gencode_by_name.index else None

    crispr_only = merged[merged['_merge'] == 'right_only'].copy()
    crispr_only['_chrom']   = crispr_only['gene_name'].map(gene_chrom)
    crispr_only['_tss']     = crispr_only['gene_name'].map(gene_tss)
    crispr_only['_biotype'] = crispr_only['gene_name'].map(gene_biotype)
    crispr_only['_gid']     = crispr_only['gene_name'].map(gene_id_base)
    crispr_only = crispr_only[crispr_only['_chrom'] == chrom]
    crispr_only = crispr_only[crispr_only['_tss'].notna()]
    crispr_only = crispr_only[(crispr_only['_tss'] - var_mid).abs() <= CIS_WIN_BP]
    crispr_only = crispr_only[crispr_only['crispr_padj'] < args.crispr_padj]
    # Down-regulated only — up-regulated CRISPR-sig hits in this category
    # are typically trans / secondary effects and don't represent the
    # canonical enhancer-positively-regulates-target relationship AG should
    # have caught. The down-only filter keeps the "AG missed these targets"
    # narrative honest.
    crispr_only = crispr_only[crispr_only['crispr_log2FoldChange'] < 0]
    # Biotype filter with INCLUDE_NON_CODING allowlist override.
    crispr_only = crispr_only[(crispr_only['_biotype'] == 'protein_coding') |
                              (crispr_only['_gid'].isin(INCLUDE_NON_CODING))]
    crispr_only['gene_id']     = crispr_only['_gid']
    crispr_only['gene_type']   = crispr_only['_biotype']
    crispr_only['gene_strand'] = crispr_only['gene_name'].map(gene_strand)

    # AG-only and Both rows retain their original gene_id; lookup TSS for them
    inner_or_left = merged[merged['_merge'].isin(['both','left_only'])].copy()
    inner_or_left['_tss'] = inner_or_left['gene_name'].map(gene_tss)

    block = pd.concat([inner_or_left, crispr_only], ignore_index=True)
    # gene_TSS + crispr_distance_kb (signed = TSS - variant_mid; negative = upstream)
    block['gene_TSS']     = block['_tss']
    block['crispr_distance_kb']  = (block['_tss'] - var_mid) / 1000.0
    block['has_crispr']   = block['crispr_log2FoldChange'].notna()
    block['has_ag']       = block['AG_RNA_raw_score'].notna()

    # Notes column
    def make_note(r):
        if r['has_ag'] and r['has_crispr']: return ''
        if r['has_ag'] and not r['has_crispr']:
            if vid == 'control_deletion_near_ATG12':
                return 'control deletion; not measured in CRISPRi assay'
            return 'not listed in CRISPRi results'
        if r['has_crispr'] and not r['has_ag']:
            d = r['crispr_distance_kb']
            side = 'downstream' if d > 0 else 'upstream'
            return f'outside AG \xb1500 kb window (TSS {abs(d):.0f} kb {side} of variant midpoint)'
        return ''
    block['notes'] = block.apply(make_note, axis=1)

    # Sort within variant: rows with CRISPR data by |log2FC| desc, then AG-only rows
    block_with_crispr    = block[block['has_crispr']].copy()
    block_with_crispr    = block_with_crispr.reindex(
        block_with_crispr['crispr_log2FoldChange'].abs().sort_values(ascending=False).index)
    block_ag_only        = block[~block['has_crispr']].copy()
    block_ag_only        = block_ag_only.reindex(
        block_ag_only['AG_RNA_raw_score'].abs().sort_values(ascending=False).index)
    block_sorted = pd.concat([block_with_crispr, block_ag_only], ignore_index=True)

    AS_ABOVE = 'as above'
    for idx, g in block_sorted.iterrows():
        # Identity columns: literal values on row 0, blank on subsequent rows.
        # Chromatin columns: literal values on row 0, 'as above' on subsequent
        # rows — these are variant-level scores so 'as above' makes clear the
        # value isn't missing, it's the same one shown at the top of the block.
        if idx == 0:
            ident_row = dict(identity); chr_row = dict(chr_cols)
        else:
            ident_row = {k: BLANK for k in identity}
            chr_row   = {k: AS_ABOVE for k in chr_cols}
        row = dict(ident_row)
        # Apply legacy-symbol override so the published table matches Ivancevic
        # 2024's gene names rather than v46's bare ENSG fallbacks.
        gid_str = str(g['gene_id']) if pd.notna(g['gene_id']) else ''
        gid_base = gid_str.split('.')[0]
        display_name = GENE_NAME_OVERRIDES.get(gid_base, g['gene_name'])
        row.update({
            'gene_name':   display_name,
            'gene_id':     g['gene_id'] if pd.notna(g['gene_id']) else BLANK,
            'gene_type':   g['gene_type'] if pd.notna(g['gene_type']) else BLANK,
            'gene_strand': g['gene_strand'] if pd.notna(g['gene_strand']) else BLANK,
            'gene_TSS':    int(g['gene_TSS']) if pd.notna(g['gene_TSS']) else BLANK,
            'crispr_distance_kb': f"{g['crispr_distance_kb']:+.1f}" if pd.notna(g['crispr_distance_kb']) else BLANK,
            'AG_RNA_raw_score':      g['AG_RNA_raw_score']      if pd.notna(g['AG_RNA_raw_score']) else BLANK,
            'AG_RNA_quantile_score': g['AG_RNA_quantile_score'] if pd.notna(g['AG_RNA_quantile_score']) else BLANK,
            'crispr_log2FoldChange': g['crispr_log2FoldChange'] if pd.notna(g['crispr_log2FoldChange']) else BLANK,
            'crispr_padj':           g['crispr_padj']           if pd.notna(g['crispr_padj']) else BLANK,
        })
        row.update(chr_row)
        row['notes'] = g['notes']
        all_rows.append(row)

out = pd.DataFrame(all_rows)
cols = ['variant_id','chrom','pos','SVLEN','type','target_gene','study_id','cell_line',
        'gene_name','gene_id','gene_type','gene_strand','gene_TSS',
        'crispr_log2FoldChange','crispr_padj','crispr_distance_kb','notes',
        'AG_RNA_raw_score','AG_RNA_quantile_score',
        'AG_ATAC_raw_score','AG_ATAC_quantile_score',
        'AG_H3K27ac_raw_score','AG_H3K27ac_quantile_score',
        'AG_H3K4me1_raw_score','AG_H3K4me1_quantile_score']
out = out[cols]
Path(args.tsv).parent.mkdir(parents=True, exist_ok=True)
out.to_csv(args.tsv, sep='\t', index=False)
print(f"Wrote {args.tsv}  ({len(out)} rows × {out.shape[1]} cols)")

# Summary diagnostics
print(f"\nPer-variant row counts:")
out_dbg = out.copy()
out_dbg['variant'] = out_dbg['variant_id'].replace('', pd.NA).ffill()
for v, g in out_dbg.groupby('variant', sort=False):
    has_ag = g['AG_RNA_raw_score'].astype(str) != ''
    has_cr = g['crispr_log2FoldChange'].astype(str) != ''
    n_both    = (has_ag & has_cr).sum()
    n_ag_only = (has_ag & ~has_cr).sum()
    n_cr_only = (~has_ag & has_cr).sum()
    print(f"  {v:35s}  total={len(g):3d}  both={n_both:3d}  AG_only={n_ag_only:3d}  CRISPR_only={n_cr_only:3d}")
