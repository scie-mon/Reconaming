#!/usr/bin/env python3
import argparse
import csv

parser = argparse.ArgumentParser()
parser.add_argument('--species-inputs', required=True)
parser.add_argument('--output', required=True)
args = parser.parse_args()

with open(args.species_inputs, newline='') as handle:
    reader = csv.DictReader(handle, delimiter='\t')
    required = {'species_id', 'input_file', 'input_type', 'role'}
    if reader.fieldnames is None or not required.issubset(reader.fieldnames):
        raise SystemExit('species_inputs must contain: species_id, input_file, input_type, role.')
    rows = list(reader)

focal = []
seen = set()
for row in rows:
    if row['role'].strip() != 'focal':
        continue
    species_id = row['species_id'].strip()
    input_file = row['input_file'].strip()
    input_type = row['input_type'].strip().lower()
    if not species_id or not input_file:
        raise SystemExit('A focal row has an empty species_id or input_file.')
    if species_id in seen:
        raise SystemExit(f'Duplicate focal species_id: {species_id}')
    if input_type != 'genome':
        raise SystemExit(
            f"Annotation mode requires input_type=genome for focal species '{species_id}', not '{input_type}'."
        )
    seen.add(species_id)
    focal.append({'species_id': species_id, 'genome_file': input_file, 'role': 'focal'})

if not focal:
    raise SystemExit('species_inputs contains no focal genome rows for annotation mode.')

with open(args.output, 'w', newline='') as handle:
    writer = csv.DictWriter(handle, fieldnames=['species_id', 'genome_file', 'role'], delimiter='\t')
    writer.writeheader()
    writer.writerows(focal)
