#!/usr/bin/env python3
import argparse
import json
import sys
from pathlib import Path

parser = argparse.ArgumentParser()
parser.add_argument('--busco-dirs', nargs='+', required=True)
parser.add_argument('--warning-threshold', type=float, required=True)
parser.add_argument('--output', required=True)
args = parser.parse_args()

rows = []
for directory_name in args.busco_dirs:
    directory = Path(directory_name)
    summaries = sorted(directory.rglob('short_summary*.json'))
    if not summaries:
        raise SystemExit(f'No BUSCO JSON summary found in {directory}')
    with summaries[-1].open() as handle:
        summary = json.load(handle)
    results = summary.get('results', {})
    complete = results.get('Complete percentage')
    single = results.get('Single copy percentage')
    duplicated = results.get('Multi copy percentage')
    if complete is None or single is None or duplicated is None:
        raise SystemExit(f'BUSCO summary in {summaries[-1]} lacks expected completeness metrics.')
    species_id = directory.name
    rows.append((species_id, float(complete), float(single), float(duplicated), str(summaries[-1])))

with open(args.output, 'w') as handle:
    handle.write('species_id\tcomplete_percent\tsingle_copy_percent\tduplicated_percent\tsummary_json\n')
    for species_id, complete, single, duplicated, summary_path in sorted(rows):
        handle.write(f'{species_id}\t{complete:.2f}\t{single:.2f}\t{duplicated:.2f}\t{summary_path}\n')
        if complete < args.warning_threshold:
            print(
                f'WARNING: BUSCO completeness for {species_id} is {complete:.2f}%, below '
                f'the configured {args.warning_threshold:.2f}% threshold.',
                file=sys.stderr,
            )
