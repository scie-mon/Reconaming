#!/usr/bin/env python3
import argparse,csv,sys
AA=set('ACDEFGHIKLMNPQRSTVWYBXZJUO*')
def read_fasta(path):
 d={}; key=None; seq=[]
 for line in open(path):
  if line.startswith('>'):
   if key:
    if key in d: sys.exit(f'Duplicate FASTA identifier: {key}')
    d[key]=''.join(seq).upper()
   key=line[1:].strip().split()[0]; seq=[]
  else: seq.append(line.strip())
 if key:
  if key in d: sys.exit(f'Duplicate FASTA identifier: {key}')
  d[key]=''.join(seq).upper()
 return d
p=argparse.ArgumentParser(); p.add_argument('--proteins',required=True); p.add_argument('--sequence-species',required=True); p.add_argument('--output-proteins',required=True); p.add_argument('--output-manifest',required=True); p.add_argument('--report',required=True); a=p.parse_args()
proteins=read_fasta(a.proteins); mapping=list(csv.DictReader(open(a.sequence_species),delimiter='\t'))
if not mapping or not {'sequence_id','species_id'} <= set(mapping[0]): sys.exit('Mapping requires sequence_id and species_id')
seen=set(); rows=[]; report=[]
for r in mapping:
 sid=r['sequence_id']; species=r['species_id']
 if sid in seen: sys.exit(f'Duplicate mapping identifier: {sid}')
 seen.add(sid)
 if sid not in proteins: sys.exit(f'Mapping identifier absent from FASTA: {sid}')
 seq=proteins[sid].rstrip('*')
 reason='' if seq and set(seq)<=AA and '*' not in seq else 'invalid_protein_sequence'
 report.append([sid,species,'accepted' if not reason else 'excluded',reason])
 if not reason: rows.append([sid,sid,sid,species,len(seq)])
extra=set(proteins)-seen
if extra: sys.exit('FASTA records lack species mapping: '+', '.join(sorted(extra)))
with open(a.output_proteins,'w') as f:
 for sid,_,_,_,_ in rows: f.write(f'>{sid}\n{proteins[sid].rstrip("*")}\n')
with open(a.output_manifest,'w',newline='') as f:
 w=csv.writer(f,delimiter='\t'); w.writerow(['sequence_id','gene_id','transcript_id','species_id','protein_length']); w.writerows(rows)
with open(a.report,'w',newline='') as f:
 w=csv.writer(f,delimiter='\t'); w.writerow(['sequence_id','species_id','status','reason']); w.writerows(report)
