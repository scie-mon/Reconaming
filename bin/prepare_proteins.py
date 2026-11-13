#!/usr/bin/env python3
import argparse
import base64
import csv
import sys
from collections import defaultdict
from urllib.parse import unquote

COMP = str.maketrans('ACGTNacgtn', 'TGCANtgcan')
CODONS = {
    'TTT':'F','TTC':'F','TTA':'L','TTG':'L','TCT':'S','TCC':'S','TCA':'S','TCG':'S',
    'TAT':'Y','TAC':'Y','TAA':'*','TAG':'*','TGT':'C','TGC':'C','TGA':'*','TGG':'W',
    'CTT':'L','CTC':'L','CTA':'L','CTG':'L','CCT':'P','CCC':'P','CCA':'P','CCG':'P',
    'CAT':'H','CAC':'H','CAA':'Q','CAG':'Q','CGC':'R','CGG':'R','CGA':'R','CGT':'R',
    'ATT':'I','ATC':'I','ATA':'I','ATG':'M','ACT':'T','ACC':'T','ACA':'T','ACG':'T',
    'AAT':'N','AAC':'N','AAA':'K','AGT':'S','AGC':'S','AGA':'R','AGG':'R',
    'GTT':'V','GTC':'V','GTA':'V','GTG':'V','GCT':'A','GCC':'A','GCA':'A','GCG':'A',
    'GAT':'D','GAC':'D','GAA':'E','GAG':'E','GGT':'G','GGC':'G','GGA':'G','GGG':'G',
}


def parse_attributes(text):
    attributes = {}
    for item in text.split(';'):
        if '=' not in item:
            continue
        key, value = item.split('=', 1)
        attributes[unquote(key)] = unquote(value)
    return attributes


def parent_ids(attributes):
    return [value.strip() for value in attributes.get('Parent', '').split(',') if value.strip()]


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
    return ''.join(CODONS.get(sequence[index:index + 3], 'X') for index in range(0, len(sequence), 3))


def sequence_id(species_id, source_feature_id):
    payload = f'{species_id}\x1f{source_feature_id}'.encode('utf-8')
    return 'seq_' + base64.urlsafe_b64encode(payload).decode('ascii').rstrip('=')


def register_feature(features, feature_id, feature_type, parents, attributes, attributes_text, line_number):
    if not feature_id:
        return
    prior = features.get(feature_id)
    if prior is None:
        features[feature_id] = {
            'type': feature_type,
            'parents': parents,
            'attributes': attributes,
            'attributes_text': attributes_text,
            'line_number': line_number,
        }
    elif (prior['type'] != feature_type or prior['parents'] != parents
          or (feature_type != 'CDS' and prior['attributes'] != attributes)):
        raise SystemExit(
            f'Inconsistent repeated GFF3 ID {feature_id!r} at line {line_number}; '
            f'first seen at line {prior["line_number"]}.'
        )


def gene_ancestor(feature_id, features):
    visited, current, gene_id = set(), feature_id, ''
    while current:
        if current in visited:
            raise SystemExit(f'Cycle in GFF3 Parent graph while resolving {feature_id!r}.')
        visited.add(current)
        feature = features.get(current)
        if feature is None:
            raise SystemExit(f'GFF3 feature {current!r} referenced while resolving {feature_id!r} is missing.')
        if feature['type'].lower() == 'gene' and not gene_id:
            gene_id = current
        if len(feature['parents']) > 1:
            raise SystemExit(
                f'Multiple Parents on GFF3 feature {current!r} at line '
                f'{feature["line_number"]} while resolving {feature_id!r}.'
            )
        current = feature['parents'][0] if feature['parents'] else ''
    return gene_id


def required_attribute(feature, feature_id, attribute, context):
    if not attribute:
        return ''
    value = feature['attributes'].get(attribute, '').strip()
    if not value:
        raise SystemExit(
            f'Requested {context} attribute {attribute!r} is absent or empty on '
            f'{feature["type"]} feature {feature_id!r}.'
        )
    if ',' in value:
        raise SystemExit(
            f'Requested {context} attribute {attribute!r} is multi-valued on '
            f'{feature["type"]} feature {feature_id!r}: {value!r}.'
        )
    return value



def corrected_translation(source_feature_id, cds_rows):
    """Resolve complete corrected proteins; never concatenate per-CDS values."""
    sequences = defaultdict(list)
    for row in cds_rows:
        value = row['attributes'].get('translation', '')
        if not value.strip():
            continue
        protein = ''.join(value.split()).upper()
        if protein.endswith('*'):
            protein = protein[:-1]
        location = f"line {row['line_number']}"
        if not protein or '*' in protein or any(
            residue not in 'ACDEFGHIKLMNPQRSTVWYBXZJUO' for residue in protein
        ):
            raise SystemExit(
                f'Invalid translation attribute for {source_feature_id!r} at {location}; '
                'expected a nonempty complete protein with no internal stops. '
                'Genomic translation was not used as fallback.'
            )
        sequences[protein].append(location)
    if len(sequences) > 1:
        locations = '; '.join(', '.join(lines) for lines in sequences.values())
        raise SystemExit(
            f'Conflicting translation attributes for {source_feature_id!r} at {locations}.'
        )
    return next(iter(sequences), None)


parser = argparse.ArgumentParser()
parser.add_argument('--genome', required=True)
parser.add_argument('--gff', required=True)
parser.add_argument('--species-id', required=True)
parser.add_argument('--proteins', required=True)
parser.add_argument('--manifest', required=True)
parser.add_argument('--report', required=True)
parser.add_argument('--gene-id-attribute', default='ID')
parser.add_argument('--transcript-id-attribute', default='ID')
parser.add_argument('--protein-id-attribute', default='')
parser.add_argument(
    '--no-parent-cds', choices=('ignore', 'fail'), default='fail',
    help='Policy for CDS without a nonempty Parent: fail (default) or ignore with warning.',
)
args = parser.parse_args()

for option in ('gene_id_attribute', 'transcript_id_attribute', 'protein_id_attribute'):
    value = getattr(args, option).strip()
    if ';' in value or '=' in value:
        raise SystemExit(f'Invalid attribute name for --{option.replace("_", "-")}: {value!r}')
    setattr(args, option, value)

if not args.transcript_id_attribute:
    raise SystemExit('--transcript-id-attribute must name a CDS attribute.')

genome = fasta(args.genome)
features = {}
cds_by_unit = defaultdict(list)
parentless_cds = []
missing_cds_ids = []
cds_unit_by_feature = {}

with open(args.gff) as handle:
    for line_number, line in enumerate(handle, start=1):
        if not line.strip() or line.startswith('#'):
            continue
        fields = line.rstrip().split('\t')
        if len(fields) != 9:
            continue
        contig, _, feature_type, start, end, _, strand, _, attributes_text = fields
        attributes = parse_attributes(attributes_text)
        feature_id = attributes.get('ID')
        parents = parent_ids(attributes)
        if feature_type == 'CDS' and not (feature_id or '').strip():
            missing_cds_ids.append(line_number)
            continue
        if feature_type == 'CDS' and not any(parent.strip() for parent in parents):
            parentless_cds.append((feature_id or '', line_number))
            continue
        register_feature(features, feature_id, feature_type, parents, attributes, attributes_text, line_number)
        if feature_type != 'CDS':
            continue
        try:
            start, end = int(start), int(end)
        except ValueError:
            start = end = None
        unit_id = required_attribute(
            {'type': 'CDS', 'attributes': attributes},
            f'{feature_id} (CDS line {line_number})',
            args.transcript_id_attribute, 'CDS candidate ID',
        )
        prior_unit = cds_unit_by_feature.setdefault(feature_id, unit_id)
        if prior_unit != unit_id:
            raise SystemExit(
                f'Inconsistent {args.transcript_id_attribute!r} values for repeated '
                f'CDS ID {feature_id!r} at line {line_number}.'
            )
        if cds_by_unit[unit_id]:
            prior = cds_by_unit[unit_id][0]
            if (prior['contig'] != contig or prior['strand'] != strand
                    or prior['parents'] != parents):
                raise SystemExit(
                    f'Inconsistent chromosome, strand or Parent for CDS candidate {unit_id!r} '
                    f'at line {line_number}; first seen at line {prior["line_number"]}.'
                )
        cds_by_unit[unit_id].append({
            'start': start, 'end': end, 'contig': contig, 'strand': strand,
            'attributes': attributes, 'line_number': line_number,
            'feature_id': feature_id, 'parents': parents,
        })

if missing_cds_ids:
    raise SystemExit(
        'CDS records missing a nonempty ID at GFF line(s): '
        + ', '.join(str(number) for number in missing_cds_ids)
        + '. CDS IDs are required, including with --no-parent-cds ignore.'
    )

accepted_rows = []
report_rows = []
seen_ids = set()

if parentless_cds:
    status = 'error' if args.no_parent_cds == 'fail' else 'excluded'
    for feature_id, line_number in parentless_cds:
        report_rows.append([
            '', feature_id, '', status,
            f'parentless_CDS;line={line_number};policy={args.no_parent_cds}', '',
        ])
    details = '\n'.join(
        f'  line {line_number}: ID={feature_id or "<missing>"}'
        for feature_id, line_number in parentless_cds
    )
    message = (
        f'{args.gff}: {len(parentless_cds)} CDS record(s) without a nonempty Parent:\n'
        f'{details}'
    )
    if args.no_parent_cds == 'fail':
        with open(args.report, 'w', newline='') as report_handle:
            writer = csv.writer(report_handle, delimiter='\t')
            writer.writerow(['sequence_id', 'source_feature_id', 'source_gene_id', 'status', 'reason', 'translation_source'])
            writer.writerows(report_rows)
        raise SystemExit(
            f'ERROR: {message}\n'
            'Fix the annotation or explicitly use --no-parent-cds ignore '
            '(pipeline: --no_parent_cds ignore) to skip these records.'
        )
    print(f'WARNING: {message}\nIgnored because --no-parent-cds ignore was supplied.', file=sys.stderr)

with open(args.proteins, 'w') as protein_handle:
    for source_feature_id in sorted(cds_by_unit):
        cds_rows = cds_by_unit[source_feature_id]
        feature_ids = sorted({row['feature_id'] for row in cds_rows})
        source_feature = features[feature_ids[0]]
        source_ids = ','.join(feature_ids)
        gene_ids = {gene_ancestor(feature_id, features) for feature_id in feature_ids}
        if len(gene_ids) != 1:
            raise SystemExit(f'Inconsistent gene membership for CDS candidate {source_feature_id!r}.')
        source_gene_id = next(iter(gene_ids))
        gene_feature = features.get(source_gene_id) if source_gene_id else None
        if args.gene_id_attribute and gene_feature is None:
            raise SystemExit(
                f'Cannot obtain requested gene identifier for {source_feature_id!r}: '
                'no ancestral gene feature was found.'
            )

        gene_key = source_gene_id or source_feature_id
        canonical_id = sequence_id(args.species_id, source_feature_id)
        if canonical_id in seen_ids:
            raise SystemExit(f'Generated duplicate sequence_id: {canonical_id}')
        seen_ids.add(canonical_id)

        gene_id = required_attribute(gene_feature, source_gene_id, args.gene_id_attribute, 'gene ID') if gene_feature else ''
        transcript_id = source_feature_id  # Configured CDS attribute, not a parent-feature identifier.
        protein_values = {
            required_attribute({'type': 'CDS', 'attributes': row['attributes']}, f'{source_feature_id} (CDS line {row["line_number"]})', args.protein_id_attribute, 'protein ID')
            for row in cds_rows
        } if args.protein_id_attribute else {''}
        if len(protein_values) != 1:
            raise SystemExit(
                f'Requested protein ID attribute {args.protein_id_attribute!r} is inconsistent among '
                f'CDS records for {source_feature_id!r}: {sorted(protein_values)!r}.'
            )
        protein_id = next(iter(protein_values))

        reason = ''
        corrected = corrected_translation(source_feature_id, cds_rows)
        translation_source = 'translation_attribute' if corrected is not None else 'genomic_CDS'
        contigs = {row['contig'] for row in cds_rows}
        strands = {row['strand'] for row in cds_rows}
        if len(contigs) != 1 or len(strands) != 1:
            reason = 'inconsistent_CDS_contig_or_strand'
        else:
            contig, strand = next(iter(contigs)), next(iter(strands))
            if contig not in genome or strand not in ('+', '-'):
                reason = 'missing_CDS_or_contig'
            else:
                cds_rows.sort(key=lambda row: row['start'] if row['start'] is not None else -1, reverse=(strand == '-'))
                fragments = []
                for row in cds_rows:
                    start, end = row['start'], row['end']
                    if start is None or end is None or start < 1 or start > end or end > len(genome[contig]):
                        reason = 'invalid_CDS_interval'
                        break
                    if corrected is None:
                        fragment = genome[contig][start - 1:end]
                        fragments.append(reverse_complement(fragment) if strand == '-' else fragment)
                if corrected is not None:
                    protein = corrected
                else:
                    coding_sequence = ''.join(fragments)
                    if not reason and len(coding_sequence) % 3:
                        reason = 'incomplete_CDS'
                    protein = translate(coding_sequence) if not reason else ''
                    if not reason and '*' in protein[:-1]:
                        reason = 'internal_stop_codon'
    
        if reason:
            report_rows.append([canonical_id, source_ids, source_gene_id, 'excluded', reason, translation_source])
            continue
        protein = protein.rstrip('*')
        protein_handle.write(f'>{canonical_id}\n{protein}\n')
        accepted_rows.append([
            canonical_id, gene_key, source_ids, ','.join(source_feature['parents']), source_gene_id,
            gene_id, transcript_id, protein_id, source_feature['type'], source_feature['attributes_text'],
            args.species_id, len(protein), len(cds_rows),
        ])
        report_rows.append([canonical_id, source_ids, source_gene_id, 'accepted', '', translation_source])

with open(args.manifest, 'w', newline='') as manifest_handle:
    writer = csv.writer(manifest_handle, delimiter='\t')
    writer.writerow([
        'sequence_id', 'gene_key', 'source_feature_id', 'source_parent_id', 'source_gene_id',
        'gene_id', 'transcript_id', 'protein_id', 'source_feature_type', 'source_attributes',
        'species_id', 'protein_length', 'cds_count',
    ])
    writer.writerows(accepted_rows)

with open(args.report, 'w', newline='') as report_handle:
    writer = csv.writer(report_handle, delimiter='\t')
    writer.writerow(['sequence_id', 'source_feature_id', 'source_gene_id', 'status', 'reason', 'translation_source'])
    writer.writerows(report_rows)
