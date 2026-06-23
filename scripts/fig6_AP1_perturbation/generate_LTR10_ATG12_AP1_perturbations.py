"""
Generate AP1-motif perturbation alleles for LTR10.ATG12 (Fig 6 input).

WT LTR10.ATG12: 2,358 bp at chr5:115,928,579-115,930,936.
Contains 20 AP1 motifs (FIMO + JASPAR FOSL1, per Ivancevic 2024 Table S25):
  - 12 in VNTR_1 (chr5:115,928,828-115,929,147)
  - 8  in VNTR_2 (chr5:115,930,258-115,930,520)
Reproducible with regex TGA[CGA]TCA = TGA[CG]TCA canonical (17) + TGAATCA
relaxed (3) = 20 hits, matching Ivancevic's count exactly.

Two perturbation classes (8 alleles total):

  TITRATION (4 alleles): scramble N of 20 motifs (N = 5, 10, 15, 20). Each
  scrambled motif is replaced with a *per-site composition-matched random
  permutation of its own 7 nucleotides* — guaranteed not to match the AP1
  regex (forward or reverse-complement strand). Per-site permutations are
  drawn with a fixed RNG seed for reproducibility. Lower-N alleles are
  strict subsets of higher-N alleles so the titration is monotonic
  (scram05's 5 sites ⊂ scram10's 10 ⊂ scram15's 15 ⊂ scram20).

  TF SUBSTITUTION (4 alleles): replace all 20 AP1 motifs with the canonical
  7-bp consensus core of TP53, CTCF, GATA1, or HNF1A. Same substitution
  applied to every motif site; flanking sequence left intact. Used to test
  whether AG's predicted regulatory output is AP1-specific or whether any
  other strong TF binder gives equivalent activity.

Why per-site composition-matched scrambles (not Ivancevic's single-string
scramble): Ivancevic 2024 used a 10-bp scramble AGCCCCGTTA to replace the
10-bp CTGAGTCACC window in their LTR10A consensus — that 10-bp window has
matched composition (4C:2T:2G:2A on both sides). For our 7-bp motif-only
substitution, a single fixed 7-mer can't match composition across both
canonical TGAGTCA (2T:2G:2A:1C) and relaxed TGAATCA (2T:1G:3A:1C) sites.
Permuting each motif's own nucleotides eliminates the composition confound
completely.

Outputs:
  --out         data/LTR10_ATG12_AP1_perturbations.tab
                MNV-format variant tab: 8 alleles, REF + ALT both 2,358 bp.
                One row per perturbation. Compatible with score_variant_lfc.py
                and score_variant_chromatin.py with --allow-mnv.
  --motif-log   data/LTR10_ATG12_AP1_perturbations.md
                Human-readable companion documenting:
                  * the 20 AP1 motif positions + sequences
                  * the titration random-pick order (seed-derived)
                  * each allele's exact motif-by-motif substitution map

Usage:
    python scripts/fig6_AP1_perturbation/generate_LTR10_ATG12_AP1_perturbations.py \\
        --variants  data/LTR10_variants.tab \\
        --out       data/LTR10_ATG12_AP1_perturbations.tab \\
        --motif-log data/LTR10_ATG12_AP1_perturbations.md
"""
import argparse
import random
import re
from itertools import permutations
from pathlib import Path
import pandas as pd

p = argparse.ArgumentParser()
p.add_argument('--variants', default='data/LTR10_variants.tab')
p.add_argument('--out',      required=True)
p.add_argument('--motif-log', required=True,
               help='Companion .md documenting the per-allele motif substitutions')
p.add_argument('--seed',     type=int, default=20260512,
               help='RNG seed for titration ordering + per-motif scramble picks')
args = p.parse_args()

# ── Load WT sequence ────────────────────────────────────────────────────────
df = pd.read_table(args.variants, sep='\t')
wt = df[df['ID'] == 'LTR10.ATG12'].iloc[0]
chrom = wt['CHROM']; pos = int(wt['POS'])
seq = wt['REF'].upper()
assert len(seq) == 2358, f"Unexpected WT length: {len(seq)}"

# ── Identify the 20 AP1 motif positions ────────────────────────────────────
# Same set Ivancevic counts: canonical TGA[CG]TCA + relaxed TGAATCA.
ap1_re = re.compile(r'(?=(TGA[CGA]TCA))')

def revcomp(s: str) -> str:
    return s.translate(str.maketrans('ACGT', 'TGCA'))[::-1]

def is_ap1_or_rc(s: str) -> bool:
    """True if s contains AP1 (TGA[CGA]TCA) on forward OR reverse strand."""
    return bool(ap1_re.search(s) or ap1_re.search(revcomp(s)))

motif_hits = []
for m in ap1_re.finditer(seq):
    s = m.start(1)
    motif_hits.append((s, s + 7, m.group(1)))
motif_hits.sort()
assert len(motif_hits) == 20, f"Expected 20 AP1 hits, got {len(motif_hits)}"

# ── Per-site composition-matched scramble (Approach 1) ─────────────────────
# For each motif, enumerate all unique permutations of its 7 letters,
# discard any that match AP1 (forward or RC), then pick one at random.
rng = random.Random(args.seed)
scrambles_per_motif: dict[int, str] = {}
for i, (start, end, original) in enumerate(motif_hits, start=1):
    # Unique permutations of this motif's nucleotides
    all_perms = sorted({''.join(p) for p in permutations(original)})
    # Drop the original (we want a CHANGE) and any AP1-containing permutation
    candidates = [p for p in all_perms
                  if p != original and not is_ap1_or_rc(p)]
    if not candidates:
        raise RuntimeError(f"No non-AP1 permutation found for motif {i} ({original})")
    scrambles_per_motif[i] = rng.choice(candidates)

# Titration order: pre-shuffle motif numbers 1..20 once; allele scramN uses
# the first N from this list. Monotonic by construction.
titration_order = list(range(1, 21))
rng.shuffle(titration_order)

# ── TF substitution sequences (all 7-bp; AP1 motifs are 7-bp) ──────────────
# Canonical JASPAR consensus cores. Same 7-mer applied to all 20 motifs.
TF_SUBS = {
    'TP53':  'GGGCATG',  # JASPAR MA0106 p53 half-site core
    'CTCF':  'CCGCGGG',  # JASPAR MA0139 C-box core
    'GATA1': 'TGATAAC',  # JASPAR MA0035 GATA half-site (extended)
    'HNF1A': 'TAATCAT',  # JASPAR MA0046 half-site core
}
for tf, s in TF_SUBS.items():
    assert len(s) == 7, tf
    assert not is_ap1_or_rc(s), f"TF sub {tf} ({s}) contains AP1; pick a different consensus"

def substitute_at(seq: str, edits: list[tuple[int, str]]) -> str:
    """Apply each (start, new_7mer) edit to seq; assumes edits don't overlap."""
    s = list(seq)
    for start, new in edits:
        assert len(new) == 7
        s[start:start + 7] = list(new)
    return ''.join(s)

# ── Build perturbation alleles ──────────────────────────────────────────────
alleles = []  # (suffix, alt_seq, list_of_motif_indices_changed, sub_label, description)

# Titration: per-motif scramble (each gets its own composition-matched perm)
for n in [5, 10, 15, 20]:
    picked_motif_idx = sorted(titration_order[:n])  # 1-based motif IDs
    edits = [(motif_hits[i - 1][0], scrambles_per_motif[i]) for i in picked_motif_idx]
    alt = substitute_at(seq, edits)
    alleles.append((
        f'scram{n:02d}', alt, picked_motif_idx,
        'per-site composition-matched scramble',
        f'{n}/20 AP1 motifs scrambled (each with its own random non-AP1 permutation of its own nucleotides)'
    ))

# TF substitutions: all 20 motifs replaced with one TF consensus
for tf, sub in TF_SUBS.items():
    picked_motif_idx = list(range(1, 21))
    edits = [(motif_hits[i - 1][0], sub) for i in picked_motif_idx]
    alt = substitute_at(seq, edits)
    alleles.append((
        f'allTF_{tf}', alt, picked_motif_idx,
        f'{tf}={sub}',
        f'all 20 AP1 motifs → {tf} core ({sub})'
    ))

# ── Sanity check ─────────────────────────────────────────────────────────────
for sfx, alt, motif_idx, _, _ in alleles:
    assert len(alt) == len(seq), sfx
    diff_positions = sum(1 for a, b in zip(seq, alt) if a != b)
    assert diff_positions <= 7 * len(motif_idx), sfx

# ── Write variant tab ───────────────────────────────────────────────────────
rows = []
for sfx, alt, motif_idx, sub_label, desc in alleles:
    rows.append({
        'ID':                       f'LTR10.ATG12_{sfx}',
        'CHROM':                    chrom,
        'POS':                      pos,
        'REF':                      seq,
        'ALT':                      alt,
        'Output':                   1,
        'Study_ID':                 'Fig6_AP1_perturbation',
        'Study_Variant_ID':         f'LTR10.ATG12_{sfx}',
        'scrambled_motif_indices':  ','.join(str(i) for i in motif_idx),
        'substitute_sequence':      sub_label,
        'description':              desc,
    })
out = pd.DataFrame(rows)
Path(args.out).parent.mkdir(parents=True, exist_ok=True)
out.to_csv(args.out, sep='\t', index=False)
print(f"Wrote {args.out}  ({len(out)} alleles)")

# ── Write motif-log markdown ───────────────────────────────────────────────
log = []
log.append('# LTR10.ATG12 AP1 perturbation alleles\n')
log.append(f'**WT element**: {chrom}:{pos:,}-{pos + len(seq) - 1:,}  ({len(seq)} bp)\n')
log.append(f'**AP1 motif detection**: regex `TGA[CGA]TCA` (matches Ivancevic 2024 Table S25 FIMO+JASPAR_FOSL1 count of 20)\n')
log.append(f'**RNG seed**: {args.seed}  (controls both the titration ordering and each motif\'s scramble permutation pick)\n')
log.append('\n## The 20 AP1 motifs in LTR10.ATG12\n')
log.append('| # | element pos | genomic pos (hg38) | motif | VNTR |\n|---|---|---|---|---|\n')
VNTR1 = (250, 568)  # element-relative VNTR_1 span (0-indexed half-open)
VNTR2 = (1680, 1941)
for i, (s, e, m) in enumerate(motif_hits, start=1):
    in_vntr = 'VNTR_1' if VNTR1[0] <= s < VNTR1[1] else ('VNTR_2' if VNTR2[0] <= s < VNTR2[1] else '-')
    log.append(f'| {i} | {s}-{e} | chr5:{pos+s:,}-{pos+e-1:,} | `{m}` | {in_vntr} |\n')

log.append('\n## Per-site scramble assignments (composition-matched permutations)\n')
log.append('| # | original | scramble | composition |\n|---|---|---|---|\n')
def comp(s: str) -> str:
    return f'{s.count("A")}A:{s.count("C")}C:{s.count("G")}G:{s.count("T")}T'
for i, (s, e, m) in enumerate(motif_hits, start=1):
    scr = scrambles_per_motif[i]
    log.append(f'| {i} | `{m}` ({comp(m)}) | `{scr}` ({comp(scr)}) | {"✓ matched" if comp(m) == comp(scr) else "✗ MISMATCH"} |\n')

log.append('\n## TF substitution consensus 7-mers\n')
log.append('| TF | substitute | JASPAR source |\n|---|---|---|\n')
log.append('| TP53  | `GGGCATG` | MA0106 p53 half-site core |\n')
log.append('| CTCF  | `CCGCGGG` | MA0139 C-box core |\n')
log.append('| GATA1 | `TGATAAC` | MA0035 GATA half-site (extended) |\n')
log.append('| HNF1A | `TAATCAT` | MA0046 half-site core |\n')

log.append('\n## Titration motif pick order (seed-derived)\n')
log.append(f'`titration_order` = {titration_order}\n\n')
log.append(f'- `scram05` scrambles motifs `{sorted(titration_order[:5])}`\n')
log.append(f'- `scram10` scrambles motifs `{sorted(titration_order[:10])}` (superset of scram05)\n')
log.append(f'- `scram15` scrambles motifs `{sorted(titration_order[:15])}` (superset of scram10)\n')
log.append(f'- `scram20` scrambles motifs `{sorted(titration_order[:20])}` (all 20)\n')

log.append('\n## Per-allele substitution map\n')
for sfx, alt, motif_idx, sub_label, desc in alleles:
    log.append(f'\n### `LTR10.ATG12_{sfx}`\n')
    log.append(f'*{desc}*\n\n')
    log.append('| motif # | original | substituted |\n|---|---|---|\n')
    for i in motif_idx:
        original_motif = motif_hits[i - 1][2]
        if sfx.startswith('scram'):
            new_seq = scrambles_per_motif[i]
        else:
            tf = sfx.replace('allTF_', '')
            new_seq = TF_SUBS[tf]
        log.append(f'| {i} | `{original_motif}` | `{new_seq}` |\n')

with open(args.motif_log, 'w') as f:
    f.writelines(log)
print(f"Wrote {args.motif_log}")

# ── Stdout summary ──────────────────────────────────────────────────────────
print('\nAlleles produced:')
for sfx, alt, motif_idx, sub_label, desc in alleles:
    print(f'  LTR10.ATG12_{sfx:14s}  motifs_changed={len(motif_idx):2d}  sub={sub_label}')
