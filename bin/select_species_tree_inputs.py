#!/usr/bin/env python3
import argparse
import csv

parser = argparse.ArgumentParser()
parser.add_argument('--representative-manifest', required=True)
parser.add_argument('--species-tree-inputs', required=True)
parser.add_argument('--busco-mode', choices=('auto', 'genome', 'proteome'), default='auto')
parser.add_argument('--output', required=True)
args = parser.parse_args()

with open(args.representative_manifest, newline='') as handle:
    manifest_rows = list(csv.DictReader(handle, delimiter='\t'))
if not manifest_rows or 'species_id' not in manifest_rows[0]:
    raise SystemExit('Representative manifest must be a non-empty TSV with a species_id column.')
gene_species = {row['species_id'].strip() for row in manifest_rows if row.get('species_id', '').strip()}
if not gene_species:
    raise SystemExit('Representative manifest contains no species_id values.')

with open(args.species_tree_inputs, newline='') as handle:
    reader = csv.DictReader(handle, delimiter='\t')
    required = {'species_id', 'input_file', 'input_type', 'role'}
    if reader.fieldnames is None or not required.issubset(reader.fieldnames):
        raise SystemExit('species_inputs must contain: species_id, input_file, input_type, role.')
    rows = list(reader)

seen = set()
focal = []
outgroups = []
for row in rows:
    species_id = row['species_id'].strip()
    input_file = row['input_file'].strip()
    input_type = row['input_type'].strip().lower()
    role = row['role'].strip()
    if not species_id or not input_file:
        raise SystemExit('species_inputs contains an empty species_id or input_file.')
    if species_id in seen:
        raise SystemExit(f'Duplicate species_id in species_inputs: {species_id}')
    seen.add(species_id)
    if input_type not in {'genome', 'proteome'}:
        raise SystemExit(f'Unsupported input_type for {species_id}: {input_type}')
    if args.busco_mode != 'auto':
        input_type = args.busco_mode
    row = {'species_id': species_id, 'input_file': input_file, 'input_type': input_type, 'role': role}
    if role == 'focal':
        focal.append(row)
    elif role == 'outgroup':
        outgroups.append(row)
    else:
        raise SystemExit(f'Unsupported role for {species_id}: {role}')

if len(outgroups) != 1:
    raise SystemExit('species_inputs must contain exactly one outgroup row.')
focal_species = {row['species_id'] for row in focal}
missing = gene_species - focal_species
extra = focal_species - gene_species
if missing:
    raise SystemExit('Missing focal BUSCO input(s): ' + ', '.join(sorted(missing)))
if extra:
    raise SystemExit('Unused focal BUSCO input(s): ' + ', '.join(sorted(extra)))
if outgroups[0]['species_id'] in gene_species:
    raise SystemExit('The auxiliary outgroup must not be a gene-family species.')

selected = sorted(focal, key=lambda row: row['species_id']) + outgroups
with open(args.output, 'w', newline='') as handle:
    writer = csv.DictWriter(handle, fieldnames=['species_id', 'input_file', 'input_type', 'role'], delimiter='\t')
    writer.writeheader()
    writer.writerows(selected)
