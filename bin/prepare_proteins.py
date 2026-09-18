#!/usr/bin/env python3
import argparse
import base64
import csv
from collections import defaultdict
from urllib.parse import unquote

COMP = str.maketrans('ACGTNacgtn', 'TGCANtgcan')
CODONS = {
    'TTT':'F','TTC':'F','TTA':'L','TTG':'L','TCT':'S','TCC':'S','TCA':'S','TCG':'S',
    'TAT':'Y','TAC':'Y','TAA':'*','TAG':'*','TGT':'C','TGC':'C','TGA':'*','TGG':'W',
    'CTT':'L','CTC':'L','CTA':'L','CTG':'L','CCT':'P','CCC':'P','CCA':'P','CCG':'P',
    'CAT':'H','CAC':'H','CAA':'Q','CAG':'Q','CGT':'R','CGC':'R','CGA':'R','CGG':'R',
    'ATT':'I','ATC':'I','ATA':'I','ATG':'M','ACT':'T','ACC':'T','ACA':'T','ACG':'T',
    'AAT':'N','AAC':'N','AAA':'K','AAG':'K','AGT':'S','AGC':'S','AGA':'R','AGG':'R',
    'GTT':'V','GTC':'V','GTA':'V','GTG':'V','GCT':'A','GCC':'A','GCA':'A','GCG':'A',
    'GAT':'D','GAC':'D','GAA':'E','GAG':'E','GGT':'G','GGC':'G','GGA':'G','GGG':'G'
}


def parse_attributes(text):
    result = {}
    for item in text.split(';'):
        if '=' in item:
            key, value = item.split('=', 1)
            result[unquote(key)] = value
    return result


def parent_ids(attributes):
    return [unquote(value) for value in attributes.get('Parent', '').split(',') if value]


def fasta(path):
    sequences, header, chunks = {}, None, []
    with open(path) as handle:
        for line in handle:
            if line.startswith('>'):
                if header:
                    sequences[header] = ''.join(chunks).upper()
                header, chunks = line[1:].split()[0], []
            else:
                chunks.append(line.strip())
    if header:
        sequences[header] = ''.join(chunks).upper()
    return sequences


def reverse_complement(sequence):
    return sequence.translate(COMP)[::-1]


def translate(sequence):
    return ''.join(CODONS.get(sequence[index:index + 3], 'X')
                   for index in range(0, len(sequence), 3))


def sequence_id(species_id, source_feature_id):
    payload = f'{species_id}\x1f{source_feature_id}'.encode('utf-8')
    return 'seq_' + base64.urlsafe_b64encode(payload).decode('ascii').rstrip('=')


def register_feature(features, feature_id, feature_type, parents, attributes_text, line_number):
    if not feature_id:
        return
    prior = features.get(feature_id)
    if prior is None:
        features[feature_id] = {
            'type': feature_type,
            'parents': parents,
            'attributes': attributes_text,
            'line_number': line_number,
        }
    elif prior['type'] != feature_type or prior['parents'] != parents:
        raise SystemExit(
            f'Inconsistent repeated GFF3 ID {feature_id!r} at line {line_number}; '
            f'first seen at line {prior["line_number"]}.'
        )


def gene_ancestor(feature_id, features):
    visited, current = set(), feature_id
    while current:
        if current in visited:
            raise SystemExit(f'Cycle in GFF3 Parent graph while resolving {feature_id}.')
        visited.add(current)
        feature = features.get(current)
        if feature is None:
            raise SystemExit(f'GFF3 feature {feature_id!r} references missing Parent/ID {current!r}.')
        if feature['type'].lower() == 'gene':
            return current
        if not feature['parents']:
            return ''
        if len(feature['parents']) != 1:
            raise SystemExit(f'Cannot resolve one gene ancestor for {feature_id!r}: multiple Parents.')
        current = feature['parents'][0]
    return ''


parser = argparse.ArgumentParser()
parser.add_argument('--genome', required=True)
parser.add_argument('--gff', required=True)
parser.add_argument('--species-id', required=True)
parser.add_argument('--proteins', required=True)
parser.add_argument('--manifest', required=True)
parser.add_argument('--report', required=True)
args = parser.parse_args()

genome = fasta(args.genome)
features = {}
cds_by_unit = defaultdict(list)

with open(args.gff) as handle:
    for line_number, line in enumerate(handle, start=1):
        if not line.strip() or line.startswith('#'):
            continue
        fields = line.rstrip().split('\t')
        if len(fields) != 9:
            continue
        contig, _, feature_type, start, end, _, strand, _, attributes_text = fields
        attributes = parse_attributes(attributes_text)
        feature_id = unquote(attributes['ID']) if 'ID' in attributes else None
        parents = parent_ids(attributes)
        register_feature(features, feature_id, feature_type, parents, attributes_text, line_number)
        if feature_type != 'CDS' or not parents:
            continue
        try:
            start, end = int(start), int(end)
        except ValueError:
            start = end = None
        for unit_id in parents:
            cds_by_unit[unit_id].append((start, end, contig, strand))

accepted_rows = []
report_rows = []
seen_ids = set()
with open(args.proteins, 'w') as protein_handle:
    for source_feature_id in sorted(cds_by_unit):
        feature = features.get(source_feature_id)
        if feature is None:
            raise SystemExit(f'CDS Parent {source_feature_id!r} has no corresponding GFF3 feature with that ID.')
        gene_id = gene_ancestor(source_feature_id, features)
        gene_key = gene_id or source_feature_id
        canonical_id = sequence_id(args.species_id, source_feature_id)
        if canonical_id in seen_ids:
            raise SystemExit(f'Generated duplicate sequence_id: {canonical_id}')
        seen_ids.add(canonical_id)

        intervals = cds_by_unit[source_feature_id]
        reason = ''
        contigs = {interval[2] for interval in intervals}
        strands = {interval[3] for interval in intervals}
        if len(contigs) != 1 or len(strands) != 1:
            reason = 'inconsistent_CDS_contig_or_strand'
        else:
            contig, strand = next(iter(contigs)), next(iter(strands))
            if contig not in genome or strand not in ('+', '-'):
                reason = 'missing_CDS_or_contig'
            else:
                intervals.sort(key=lambda interval: interval[0] if interval[0] is not None else -1,
                               reverse=(strand == '-'))
                fragments = []
                for start, end, cds_contig, cds_strand in intervals:
                    if (start is None or end is None or cds_contig != contig or cds_strand != strand or
                            start < 1 or start > end or end > len(genome[contig])):
                        reason = 'invalid_CDS_interval'
                        break
                    fragment = genome[contig][start - 1:end]
                    fragments.append(reverse_complement(fragment) if strand == '-' else fragment)
                coding_sequence = ''.join(fragments)
                if not reason and len(coding_sequence) % 3:
                    reason = 'incomplete_CDS'
                protein = translate(coding_sequence) if not reason else ''
                if not reason and '*' in protein[:-1]:
                    reason = 'internal_stop_codon'

        if reason:
            report_rows.append([canonical_id, source_feature_id, gene_id, 'excluded', reason])
            continue
        protein = protein.rstrip('*')
        protein_handle.write(f'>{canonical_id}\n{protein}\n')
        accepted_rows.append([
            canonical_id, gene_key, source_feature_id, ','.join(feature['parents']), gene_id,
            feature['type'], feature['attributes'], args.species_id, len(protein), len(intervals),
        ])
        report_rows.append([canonical_id, source_feature_id, gene_id, 'accepted', ''])

with open(args.manifest, 'w', newline='') as manifest_handle:
    writer = csv.writer(manifest_handle, delimiter='\t')
    writer.writerow([
        'sequence_id', 'gene_key', 'source_feature_id', 'source_parent_id', 'gene_id',
        'source_feature_type', 'source_attributes', 'species_id', 'protein_length', 'cds_count',
    ])
    writer.writerows(accepted_rows)

with open(args.report, 'w', newline='') as report_handle:
    writer = csv.writer(report_handle, delimiter='\t')
    writer.writerow(['sequence_id', 'source_feature_id', 'gene_id', 'status', 'reason'])
    writer.writerows(report_rows)
