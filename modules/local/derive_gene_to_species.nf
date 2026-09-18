process DERIVE_GENE_TO_SPECIES {
    tag 'derive gene-to-species mapping'

    input:
    path representative_manifest

    output:
    path 'gene_to_species.tsv', emit: gene_to_species

    script:
    """
python3 - '${representative_manifest}' <<'PY'
import csv
import sys

TAB = chr(9)
NEWLINE = chr(10)

with open(sys.argv[1], newline='') as handle:
    rows = list(csv.DictReader(handle, delimiter=TAB))

if not rows:
    raise SystemExit('Representative manifest is empty.')

required = {'sequence_id', 'species_id'}
missing = required - set(rows[0])
if missing:
    raise SystemExit(f'Representative manifest lacks required column(s): {", ".join(sorted(missing))}')

seen_sequences = set()
with open('gene_to_species.tsv', 'w', newline='') as handle:
    writer = csv.writer(handle, delimiter=TAB, lineterminator=NEWLINE)
    writer.writerow(['sequence_id', 'species_id'])
    for row in rows:
        sequence_id = row['sequence_id'].strip()
        species_id = row['species_id'].strip()
        if not sequence_id or not species_id:
            raise SystemExit('Representative manifest contains an empty sequence_id or species_id.')
        if sequence_id in seen_sequences:
            raise SystemExit(f'Representative manifest contains duplicate sequence_id: {sequence_id}')
        seen_sequences.add(sequence_id)
        writer.writerow([sequence_id, species_id])
PY
    """
}
