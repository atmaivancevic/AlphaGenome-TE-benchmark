"""
Group ENCODE H3K27ac peak files by AG-training biosample. Reads the resolution
TSV + AG track manifest and prints one line per human biosample to stdout. Normally called by
merge_encode_h3k27ac_peaks.sh, not run directly.

Example usage:
python scripts/fig1_baseline_chromatin/group_peaks_by_biosample.py \
    --enriched data/encode_h3k27ac/peak_resolution_human.tsv \
    --manifest data/AG_training_H3K27ac_tracks.csv \
    --peakdir  data/encode_h3k27ac/peaks/human
"""

import argparse
import csv

parser = argparse.ArgumentParser()
parser.add_argument('--enriched', required=True,
                    help='Resolution TSV (one row per AG-training bigWig with peak_ENCFF)')
parser.add_argument('--manifest', required=True,
                    help='AG training-track manifest CSV (one row per biosample)')
parser.add_argument('--peakdir', required=True,
                    help='Directory containing per-ENCSR peak files named {ENCSR}_{peak_ENCFF}.bed.gz')
args = parser.parse_args()

# Build (ENCSR, bigwig) -> peak_ENCFF lookup from the enriched resolution.
pk = {}
with open(args.enriched) as f:
    for r in csv.DictReader(f, delimiter='\t'):
        if r['peak_ENCFF']:
            pk[(r['ENCSR'], r['AG_bigwig_ENCFF'])] = r['peak_ENCFF']

# Read each biosample row from the manifest. For biosamples with multiple
# bigWigs, the Experiment/File accession fields hold comma-separated lists.
with open(args.manifest) as f:
    for r in csv.DictReader(f):
        if r['organism'] != 'human':
            continue
        tag = r['ontology_curie'].replace(':', '_')
        name = (r['biosample_name']
                .replace(' ', '_')
                .replace(',', '')
                .replace("'", '')
                .replace('.', '_'))
        encsrs = [x.strip() for x in r['Experiment accession'].split(',') if x.strip()]
        bigwigs = [x.strip() for x in r['File accession'].split(',') if x.strip()]
        files = []
        for encsr, bigwig in zip(encsrs, bigwigs):
            peak = pk.get((encsr, bigwig))
            if peak:
                files.append(f'{args.peakdir}/{encsr}_{peak}.bed.gz')
        print(f'{tag}_{name}\t' + '\t'.join(files))
