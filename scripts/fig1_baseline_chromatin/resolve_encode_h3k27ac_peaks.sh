#!/usr/bin/env bash
# For each AG H3K27ac bigWig (human), find the matching default bed narrowPeak
# file from the SAME ENCODE analysis. Writes a TSV; downloads nothing.
# Match rule: same ENCSR, shares an analysis (ENCAN) with the AG bigWig,
# file_format=bed narrowPeak, preferred_default=true.
set -euo pipefail

MANIFEST="${1:-data/AG_training_H3K27ac_tracks_passwarning.csv}"
LIMIT="${2:-0}"
export LIMIT
OUT="data/encode_h3k27ac/peak_resolution_human.tsv"
mkdir -p "$(dirname "$OUT")"

echo -e "ENCSR\tAG_bigwig_ENCFF\tENCAN\tpeak_ENCFF\tpeak_output_type\tpeak_assembly\tpreferred_default\tnote" > "$OUT"

python3 - "$MANIFEST" "$OUT" <<'PY'
import csv, json, sys, os, urllib.request, urllib.parse, time
LIMIT = int(os.environ.get("LIMIT","0"))

manifest, out = sys.argv[1], sys.argv[2]
BASE = "https://www.encodeproject.org"

def fetch(path):
    req = urllib.request.Request(BASE + path, headers={"Accept": "application/json"})
    with urllib.request.urlopen(req, timeout=30) as r:
        return json.load(r)

rows = []
with open(manifest) as f:
    for r in csv.DictReader(f):
        if r["organism"] != "human":
            continue
        encsrs = [x.strip() for x in r["Experiment accession"].split(",") if x.strip()]
        encffs = [x.strip() for x in r["File accession"].split(",") if x.strip()]
        for e, bw in zip(encsrs, encffs):
            rows.append((e, bw))

if LIMIT: rows = rows[:LIMIT]
print(f"[info] {len(rows)} (ENCSR, bigWig) pairs", file=sys.stderr)

with open(out, "a") as fo:
    for i, (encsr, bw) in enumerate(rows, 1):
        try:
            bwj = fetch(f"/files/{bw}/?format=json")
            def ext(a):
                if isinstance(a, dict):
                    return a.get("@id","").rstrip("/").split("/")[-1]
                return a.rstrip("/").split("/")[-1]
            encans = [x for x in (ext(a) for a in bwj.get("analyses", [])) if x]
            if not encans:
                fo.write(f"{encsr}\t{bw}\t\t\t\t\t\tno_analysis_on_bigwig\n")
                print(f"[skip {i}/{len(rows)}] {encsr} {bw}: no analysis", file=sys.stderr)
                continue
            expj = fetch(f"/experiments/{encsr}/?format=json")
            def passes_common(f):
                if f.get("file_format") != "bed": return False
                if f.get("file_format_type") != "narrowPeak": return False
                if f.get("status") != "released": return False
                if f.get("assembly") not in ("GRCh38","mm10"): return False
                dc = f.get("date_created","")
                if dc[:4] < "2020": return False
                # exclude files flagged for missing chromosome data
                audits = f.get("audit", {}) or {}
                for level in ("ERROR","NOT_COMPLIANT","WARNING","INTERNAL_ACTION"):
                    for a in audits.get(level, []) or []:
                        cat = (a.get("category","") + " " + a.get("detail","")).lower()
                        if "missing chromosom" in cat or "missing chr " in cat:
                            return False
                return True
            exp_peaks = [f for f in expj.get("files", []) if passes_common(f)]
            cand = []
            for f in exp_peaks:
                f_encans = [x for x in (ext(a) for a in f.get("analyses", [])) if x]
                if not any(a in encans for a in f_encans): continue
                cand.append(f)
            note_prefix = ""
            if not cand:
                # fallback: any preferred_default peak on the experiment (diff analysis)
                cand = [f for f in exp_peaks if f.get("preferred_default")]
                note_prefix = "diff_analysis_" if cand else ""
                if not cand:
                    fo.write(f"{encsr}\t{bw}\t{','.join(encans)}\t\t\t\t\tno_peak_in_analysis\n")
                    print(f"[skip {i}/{len(rows)}] {encsr}: no matching peak", file=sys.stderr)
                    continue
            # Prefer preferred_default=true; else replicated > pseudoreplicated > peaks
            pref = [f for f in cand if f.get("preferred_default")]
            pick = None
            if pref:
                pick = pref[0]
                note = note_prefix + "default"
            else:
                order = {"replicated peaks":0,"pseudoreplicated peaks":1,"peaks":2}
                cand.sort(key=lambda f: order.get(f.get("output_type",""), 9))
                pick = cand[0]
                note = note_prefix + "no_default_fallback"
            pacc = pick["accession"]
            otype = pick.get("output_type","")
            asm = pick.get("assembly","")
            pd = "true" if pick.get("preferred_default") else "false"
            fo.write(f"{encsr}\t{bw}\t{','.join(encans)}\t{pacc}\t{otype}\t{asm}\t{pd}\t{note}\n")
            fo.flush()
            print(f"[ok {i}/{len(rows)}] {encsr} bw={bw} -> peak={pacc} ({otype}, default={pd})", file=sys.stderr)
        except Exception as ex:
            fo.write(f"{encsr}\t{bw}\t\t\t\t\t\terror:{ex}\n")
            print(f"[err {i}/{len(rows)}] {encsr}: {ex}", file=sys.stderr)
        time.sleep(0.05)
PY

echo "[done] wrote $OUT"
