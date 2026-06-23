"""
Generate AP1-motif ADDITION alleles for LTR10.ATG12 (Fig 6 Panel A, past WT).
Companion to generate_LTR10_ATG12_AP1_perturbations.py. Overwrites non-AP1 7-bp
windows of the WT element with TGAGTCA in place — length-preserving substitution
variants (REF + ALT both 2,358 bp), no coordinate shifts. Adds 5..220 new motifs
(totals 25..240); new positions are non-overlapping, seeded, lower-N subsets of higher-N.

Example usage:
python scripts/fig6_AP1_perturbation/generate_LTR10_ATG12_AP1_insertions.py \
    --variants  data/LTR10_variants.tab \
    --out       data/LTR10_ATG12_AP1_insertions.tab \
    --motif-log data/LTR10_ATG12_AP1_insertions.md
"""
import argparse
import random
import re
from pathlib import Path
import pandas as pd

p = argparse.ArgumentParser()
p.add_argument('--variants',  default='data/LTR10_variants.tab')
p.add_argument('--out',       required=True)
p.add_argument('--motif-log', required=True)
p.add_argument('--seed',      type=int, default=20260513,
               help='RNG seed for new-motif position picking')
args = p.parse_args()

# ── Load WT element ─────────────────────────────────────────────────────────
df = pd.read_table(args.variants, sep='\t')
wt = df[df['ID'] == 'LTR10.ATG12'].iloc[0]
chrom = wt['CHROM']; pos = int(wt['POS'])
seq = wt['REF'].upper()
assert len(seq) == 2358

# ── Helpers: AP1 motif detection (forward + reverse complement) ────────────
ap1_re = re.compile(r'(?=(TGA[CGA]TCA))')

def revcomp(s: str) -> str:
    return s.translate(str.maketrans('ACGT', 'TGCA'))[::-1]

def is_ap1_or_rc(s: str) -> bool:
    return bool(ap1_re.search(s) or ap1_re.search(revcomp(s)))

# Existing AP1 motif intervals in WT
existing = [(m.start(1), m.start(1) + 7) for m in ap1_re.finditer(seq)]
assert len(existing) == 20, f"Expected 20 AP1 motifs in WT, got {len(existing)}"
existing_set = set()
for s, e in existing:
    existing_set.update(range(s, e))

# ── Build candidate positions for new motifs ───────────────────────────────
# A candidate position p is valid if:
#   (a) the 7-bp window [p, p+7) does not overlap any existing AP1 motif
#   (b) the window does not already match AP1 forward/RC by chance
candidates = []
for p_start in range(len(seq) - 7 + 1):
    win_idxs = set(range(p_start, p_start + 7))
    if win_idxs & existing_set:
        continue
    if is_ap1_or_rc(seq[p_start:p_start + 7]):
        continue
    candidates.append(p_start)
print(f"WT has {len(existing)} AP1 motifs at positions: {[s for s,e in existing]}")
print(f"{len(candidates)} non-AP1 7-bp positions available for new motifs")

# ── Pick non-overlapping positions: shuffle then greedy ────────────────────
# Random shuffle (seeded) then greedy accept. Cap at 220 picks; lower-N
# alleles are strict prefixes. Random shuffle spreads new motifs across
# the full element length instead of packing them at the left edge.
rng = random.Random(args.seed)
shuffled = list(candidates); rng.shuffle(shuffled)
picked: list[int] = []
picked_idx_set: set[int] = set()
NEW_MOTIF = 'TGAGTCA'
for p_start in shuffled:
    win_idxs = set(range(p_start, p_start + 7))
    if win_idxs & picked_idx_set:
        continue
    picked.append(p_start)
    picked_idx_set |= win_idxs
    if len(picked) >= 220:
        break

if len(picked) < 220:
    raise RuntimeError(f"Couldn't pick 220 non-overlapping new motif positions "
                       f"(only got {len(picked)}). Reduce add220 or try a different seed.")

print(f"Picked {len(picked)} non-overlapping positions for new motifs")

def apply_new_motifs(seq: str, positions: list[int]) -> str:
    s = list(seq)
    for p_start in positions:
        s[p_start:p_start + 7] = list(NEW_MOTIF)
    return ''.join(s)

# ── Generate alleles ────────────────────────────────────────────────────────
ADD = [
    ('add05',     5),  # off-grid (25 total) — pre-existing fine-grained point
    ('add10',    10),  # 30  total
    ('add20',    20),  # 40
    ('add30',    30),  # 50
    ('add40',    40),  # 60
    ('add50',    50),  # 70
    ('add60',    60),  # 80
    ('add70',    70),  # 90
    ('add80',    80),  # 100
    ('add90',    90),  # 110
    ('add100',  100),  # 120
    ('add110',  110),  # 130
    ('add120',  120),  # 140
    ('add130',  130),  # 150
    ('add140',  140),  # 160
    ('add150',  150),  # 170
    ('add160',  160),  # 180
    ('add170',  170),  # 190
    ('add180',  180),  # 200
    ('add190',  190),  # 210
    ('add200',  200),  # 220
    ('add210',  210),  # 230
    ('add220',  220),  # 240
]
rows = []
for suffix, n_new in ADD:
    new_positions = picked[:n_new]  # subset (monotonic)
    alt = apply_new_motifs(seq, new_positions)
    total_motifs = 20 + n_new
    # Verify the result actually has the expected motif count
    actual = len(list(ap1_re.finditer(alt)))
    assert actual == total_motifs, \
        f"{suffix}: expected {total_motifs} AP1 motifs in ALT, got {actual}"
    rows.append({
        'ID':                  f'LTR10.ATG12_{suffix}',
        'CHROM':               chrom,
        'POS':                 pos,
        'REF':                 seq,
        'ALT':                 alt,
        'Output':              1,
        'Study_ID':            'Fig6_AP1_addition',
        'Study_Variant_ID':    f'LTR10.ATG12_{suffix}',
        'n_new_motifs':        n_new,
        'total_motifs':        total_motifs,
        'new_motif_positions': ','.join(str(p) for p in sorted(new_positions)),
        'description':         f'{n_new} new TGAGTCA motifs written into non-AP1 7-bp '
                               f'windows (substitution variant, same length as WT) → {total_motifs} motifs total'
    })

out = pd.DataFrame(rows)
Path(args.out).parent.mkdir(parents=True, exist_ok=True)
out.to_csv(args.out, sep='\t', index=False)
print(f"\nWrote {args.out}  ({len(out)} alleles)")
print(out[['ID','POS','n_new_motifs','total_motifs']].to_string(index=False))

# ── Companion .md log ──────────────────────────────────────────────────────
log = []
log.append('# LTR10.ATG12 AP1 ADDITION alleles\n')
log.append(f'**WT element**: {chrom}:{pos:,}-{pos+len(seq)-1:,}  ({len(seq)} bp; 20 AP1 motifs in two VNTR arrays)\n')
log.append(f'**Strategy**: write new TGAGTCA motifs in place over non-AP1 7-bp windows (substitution variant — same length as WT)\n')
log.append(f'**RNG seed**: {args.seed}  (controls which non-AP1 windows are picked)\n')

log.append(f'\n## Candidate pool\n')
log.append(f'- 7-bp windows in WT that overlap no existing AP1 motif AND do not '
           f'themselves contain TGA[CGA]TCA by chance: **{len(candidates)}** windows\n')
log.append(f'- After non-overlapping greedy selection from the shuffled pool: '
           f'**{len(picked)}** positions reserved for new motifs.\n')

log.append(f'\n## Pick order (seed-derived; lower-N alleles are subsets of higher-N)\n')
log.append(f'`picked_positions` (first 60, monotonic) = `{picked[:60]}`\n')

log.append(f'\n## ADDITION alleles\n')
log.append('| Allele | New motifs | Total motifs | New motif positions in element |\n|---|---:|---:|---|\n')
for r in rows:
    log.append(f'| `{r["ID"]}` | {r["n_new_motifs"]} | {r["total_motifs"]} | `{r["new_motif_positions"]}` |\n')

log.append(f'\n## File format\n')
log.append('REF and ALT are both 2,358 bp (substitution variant — same length, different content at the new-motif positions).\n')
log.append('Compatible with score_variant_lfc.py / score_variant_chromatin.py / predict_variant_tracks.py\n')
log.append('via the existing `--allow-sub` flag (same as the scrambling perturbations).\n')

with open(args.motif_log, 'w') as f:
    f.writelines(log)
print(f"Wrote {args.motif_log}")
