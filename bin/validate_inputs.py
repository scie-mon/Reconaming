#!/usr/bin/env python3
import argparse, csv, os, sys

def die(message):
    print(f'ERROR: {message}', file=sys.stderr); sys.exit(2)
def present(value): return bool(value and value.strip())
def main(a):
    annotation=present(a.genome_map) or present(a.gff)
    protein=present(a.protein_fasta)
    if annotation == protein: die('Provide exactly one mode: --genome-map with --gff, or --protein-fasta')
    mode='annotation' if annotation else 'protein'
    required=['genome_map','gff'] if mode=='annotation' else ['protein_fasta','sequence_mapping']
    for key in required:
        value=getattr(a,key)
        if not present(value) or not os.path.isfile(value): die(f'Missing required {key}: {value}')
    if mode=='annotation':
        with open(a.genome_map,newline='') as f: rows=list(csv.DictReader(f,delimiter='\t'))
        columns=set(rows[0]) if rows else set()
        if not {'genome_file','species_id','role'} <= columns: die('genome map requires genome_file, species_id, role')
        if not any(r['role']=='outgroup' for r in rows): die('genome map requires one outgroup row')
        if not any(r['role']=='focal' for r in rows): die('genome map requires one focal row')
        for r in rows:
            if r['role'] not in {'focal','outgroup'}: die(f'Invalid role: {r["role"]}')
            if not os.path.isfile(r['genome_file']): die(f'Genome file absent: {r["genome_file"]}')
    with open(a.settings,'w',newline='') as f:
        w=csv.writer(f,delimiter='\t'); w.writerow(['key','value'])
        for k,v in [('input_mode',mode),('translation_table','1'),('invalid_model_policy','exclude_and_report'),('protein_input_is_representative',str(mode=='protein').lower())]: w.writerow([k,v])
    with open(a.report,'w',newline='') as f:
        w=csv.writer(f,delimiter='\t'); w.writerow(['status','input_mode','message']); w.writerow(['accepted',mode,'input contract validated'])
p=argparse.ArgumentParser(); p.add_argument('--genome-map'); p.add_argument('--gff'); p.add_argument('--protein-fasta'); p.add_argument('--sequence-mapping'); p.add_argument('--settings',required=True); p.add_argument('--report',required=True)
if __name__=='__main__': main(p.parse_args())
