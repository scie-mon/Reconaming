#!/usr/bin/env python3
import argparse, csv, os, sys

def records(path):
    header = None; sequence = []
    with open(path) as fh:
        for line in fh:
            if line.startswith('>'):
                if header is not None: yield header, ''.join(sequence)
                header = line[1:].rstrip('\n'); sequence = []
            else:
                sequence.append(line.strip())
    if header is not None: yield header, ''.join(sequence)

def main(a):
    with open(a.genome_map, newline='') as fh:
        rows = list(csv.DictReader(fh, delimiter='\t'))
    required = {'genome_file','species_id','role'}
    if not rows or not required.issubset(rows[0]): sys.exit('genome_to_species TSV requires genome_file, species_id, role')
    contigs = set(); manifest = []
    with open(a.concatenated_fasta, 'w') as out:
        for row in rows:
            path = row['genome_file']
            if not os.path.isabs(path): path = os.path.join(a.base_dir, path)
            path = os.path.abspath(path)
            if not os.path.isfile(path): sys.exit(f'Missing genome FASTA: {row["genome_file"]}')
            for header, sequence in records(path):
                contig = header.split()[0]
                if contig in contigs: sys.exit(f'Duplicate FASTA record identifier: {contig}')
                contigs.add(contig)
                out.write(f'>{header}\n{sequence}\n')
                manifest.append([row['genome_file'], row['species_id'], row['role'], contig, header, len(sequence)])
    with open(a.contig_map, 'w', newline='') as fh:
        writer = csv.writer(fh, delimiter='\t')
        writer.writerow(['genome_file','species_id','role','contig_id','original_defline','length'])
        writer.writerows(manifest)

p = argparse.ArgumentParser()
p.add_argument('--genome-map', required=True); p.add_argument('--base-dir', required=True)
p.add_argument('--concatenated-fasta', required=True); p.add_argument('--contig-map', required=True)
if __name__ == '__main__': main(p.parse_args())
