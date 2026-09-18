#!/usr/bin/env python3
import argparse
import csv
import statistics
import sys


def fasta(path):
    sequences = {}
    header = None
    chunks = []
    with open(path) as handle:
        for line in handle:
            if line.startswith('>'):
                if header is not None:
                    if header in sequences:
                        sys.exit(f'Duplicate FASTA identifier: {header}')
                    sequences[header] = ''.join(chunks)
                header = line[1:].strip().split()[0]
                chunks = []
            else:
                chunks.append(line.strip())
    if header is not None:
        if header in sequences:
            sys.exit(f'Duplicate FASTA identifier: {header}')
        sequences[header] = ''.join(chunks)
    return sequences


def main(args):
    with open(args.manifest, newline='') as handle:
        rows = list(csv.DictReader(handle, delimiter='\t'))

    required = {'sequence_id', 'gene_key', 'protein_length'}
    if not rows or not required.issubset(rows[0]):
        sys.exit('Manifest lacks required columns: sequence_id, gene_key, protein_length')

    proteins = fasta(args.proteins)
    by_gene_key = {}
    sequence_lookup = {}
    for row in rows:
        sequence_id = row['sequence_id'].strip()
        gene_key = row['gene_key'].strip()
        if not sequence_id or not gene_key:
            sys.exit('Manifest contains an empty sequence_id or gene_key.')
        if sequence_id in sequence_lookup:
            sys.exit(f'Manifest contains duplicate sequence_id: {sequence_id}')
        if sequence_id not in proteins:
            sys.exit(f'Manifest sequence absent from FASTA: {sequence_id}')
        try:
            if int(row['protein_length']) < 1:
                raise ValueError
        except ValueError:
            sys.exit(f'Invalid protein_length for sequence_id: {sequence_id}')
        sequence_lookup[sequence_id] = row
        by_gene_key.setdefault(gene_key, []).append(row)

    extra_fasta_ids = set(proteins) - set(sequence_lookup)
    if extra_fasta_ids:
        sys.exit(f'FASTA contains IDs absent from manifest: {", ".join(sorted(extra_fasta_ids)[:10])}')

    chosen = {}
    mode = 'median'
    if args.selection_table:
        with open(args.selection_table, newline='') as handle:
            selected = list(csv.DictReader(handle, delimiter='\t'))
        if not selected or 'sequence_id' not in selected[0]:
            sys.exit('Selection table requires a sequence_id column.')
        for selected_row in selected:
            sequence_id = selected_row['sequence_id'].strip()
            if sequence_id not in sequence_lookup:
                sys.exit(f'Unknown selected sequence_id: {sequence_id}')
            row = sequence_lookup[sequence_id]
            gene_key = row['gene_key']
            if gene_key in chosen:
                sys.exit(f'Multiple selected isoforms for gene group: {gene_key}')
            chosen[gene_key] = row
        missing = set(by_gene_key) - set(chosen)
        if missing:
            sys.exit('Selection table does not select every gene group.')
        mode = 'explicit'
    else:
        median_length = statistics.median(int(row['protein_length']) for row in rows)
        for gene_key, group_rows in by_gene_key.items():
            chosen[gene_key] = sorted(
                group_rows,
                key=lambda row: (
                    abs(int(row['protein_length']) - median_length),
                    -int(row['protein_length']),
                    row['sequence_id'],
                ),
            )[0]

    with open(args.representatives, 'w') as handle:
        for gene_key, row in sorted(chosen.items()):
            sequence_id = row['sequence_id']
            handle.write(f'>{sequence_id}\n{proteins[sequence_id]}\n')

    with open(args.output_manifest, 'w', newline='') as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=list(rows[0]) + ['selection_mode'],
            delimiter='\t',
            lineterminator='\n',
        )
        writer.writeheader()
        for row in sorted(chosen.values(), key=lambda item: (item['gene_key'], item['sequence_id'])):
            output_row = dict(row)
            output_row['selection_mode'] = mode
            writer.writerow(output_row)

    with open(args.report, 'w', newline='') as handle:
        writer = csv.writer(handle, delimiter='\t', lineterminator='\n')
        writer.writerow(['gene_key', 'selected_sequence_id', 'selection_mode'])
        for gene_key, row in sorted(chosen.items()):
            writer.writerow([gene_key, row['sequence_id'], mode])


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--proteins', required=True)
    parser.add_argument('--manifest', required=True)
    parser.add_argument('--selection-table')
    parser.add_argument('--representatives', required=True)
    parser.add_argument('--output-manifest', required=True)
    parser.add_argument('--report', required=True)
    main(parser.parse_args())
