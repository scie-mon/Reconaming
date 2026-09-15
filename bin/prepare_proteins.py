#!/usr/bin/env python3
import argparse
import csv
from collections import defaultdict

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


def attr(text):
    return {item.split('=', 1)[0]: item.split('=', 1)[1]
            for item in text.split(';') if '=' in item}


def fasta(path):
    sequences = {}
    header = None
    chunks = []
    with open(path) as handle:
        for line in handle:
            if line.startswith('>'):
                if header:
                    sequences[header] = ''.join(chunks).upper()
                header = line[1:].split()[0]
                chunks = []
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


parser = argparse.ArgumentParser()
parser.add_argument('--genome', required=True)
parser.add_argument('--gff', required=True)
parser.add_argument('--species-id', required=True)
parser.add_argument('--proteins', required=True)
parser.add_argument('--manifest', required=True)
parser.add_argument('--report', required=True)
args = parser.parse_args()

genome = fasta(args.genome)
transcripts = {}
cds_by_transcript = defaultdict(list)

with open(args.gff) as handle:
    for line in handle:
        if not line.strip() or line.startswith('#'):
            continue
        fields = line.rstrip().split('\t')
        if len(fields) != 9:
            continue
        contig, _, feature_type, start, end, _, strand, _, attributes = fields
        attributes = attr(attributes)
        feature_id = attributes.get('ID')
        if feature_type in ('mRNA', 'transcript') and feature_id and attributes.get('Parent'):
            transcripts[feature_id] = (
                attributes['Parent'], contig, strand, attributes.get('internal_id', '')
            )
        elif feature_type == 'CDS' and attributes.get('Parent'):
            try:
                start = int(start)
                end = int(end)
            except ValueError:
                start = end = None
            cds_by_transcript[attributes['Parent']].append((start, end, contig, strand))

accepted_rows = []
report_rows = []
with open(args.proteins, 'w') as protein_handle:
    for transcript_id, (gene_id, contig, strand, internal_id) in sorted(transcripts.items()):
        intervals = cds_by_transcript.get(transcript_id, [])
        reason = ''
        if not intervals or contig not in genome or strand not in ('+', '-'):
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
            report_rows.append([transcript_id, gene_id, 'excluded', reason])
            continue
        protein = protein.rstrip('*')
        protein_handle.write(f'>{transcript_id}\n{protein}\n')
        accepted_rows.append([
            transcript_id, gene_id, transcript_id, contig, args.species_id,
            len(protein), internal_id
        ])
        report_rows.append([transcript_id, gene_id, 'accepted', ''])

with open(args.manifest, 'w', newline='') as manifest_handle:
    writer = csv.writer(manifest_handle, delimiter='\t')
    writer.writerow([
        'sequence_id', 'gene_id', 'transcript_id', 'contig_id', 'species_id',
        'protein_length', 'irnotator_internal_id'
    ])
    writer.writerows(accepted_rows)

with open(args.report, 'w', newline='') as report_handle:
    writer = csv.writer(report_handle, delimiter='\t')
    writer.writerow(['transcript_id', 'gene_id', 'status', 'reason'])
    writer.writerows(report_rows)