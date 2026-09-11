#!/usr/bin/env python3
import argparse, csv, statistics, sys

def fasta(path):
    d={}; h=None; s=[]
    for line in open(path):
        if line.startswith('>'):
            if h: d[h]=''.join(s)
            h=line[1:].strip().split()[0]; s=[]
        else: s.append(line.strip())
    if h: d[h]=''.join(s)
    return d

def main(a):
    rows=list(csv.DictReader(open(a.manifest),delimiter='\t'))
    required={'sequence_id','gene_id','transcript_id','protein_length'}
    if not rows or not required <= set(rows[0]): sys.exit('Manifest lacks required columns')
    proteins=fasta(a.proteins); by_gene={}
    for r in rows: by_gene.setdefault(r['gene_id'],[]).append(r)
    chosen={}; mode='median'
    if a.selection_table:
        selected=list(csv.DictReader(open(a.selection_table),delimiter='\t'))
        key='sequence_id' if selected and 'sequence_id' in selected[0] else 'transcript_id'
        if not selected or key not in selected[0]: sys.exit('Selection table requires sequence_id or transcript_id')
        lookup={r[key]:r for rs in by_gene.values() for r in rs}
        for s in selected:
            if s[key] not in lookup: sys.exit(f'Unknown selected ID: {s[key]}')
            r=lookup[s[key]]
            if r['gene_id'] in chosen: sys.exit(f'Multiple selected isoforms for gene: {r["gene_id"]}')
            chosen[r['gene_id']]=r
        missing=set(by_gene)-set(chosen)
        if missing: sys.exit('Selection table does not select every gene')
        mode='explicit'
    else:
        median=statistics.median([int(r['protein_length']) for r in rows])
        for gene,rs in by_gene.items():
            chosen[gene]=sorted(rs,key=lambda r:(abs(int(r['protein_length'])-median),-int(r['protein_length']),r['transcript_id']))[0]
    with open(a.representatives,'w') as f:
        for gene,r in sorted(chosen.items()):
            seq=proteins.get(r['sequence_id'])
            if not seq: sys.exit(f'Manifest sequence absent from FASTA: {r["sequence_id"]}')
            f.write(f'>{r["sequence_id"]}\n{seq}\n')
    with open(a.output_manifest,'w',newline='') as f:
        w=csv.DictWriter(f,fieldnames=list(rows[0])+['selection_mode'],delimiter='\t'); w.writeheader()
        for r in sorted(chosen.values(),key=lambda x:x['gene_id']): r['selection_mode']=mode; w.writerow(r)
    with open(a.report,'w',newline='') as f:
        w=csv.writer(f,delimiter='\t'); w.writerow(['gene_id','selected_sequence_id','selection_mode'])
        for gene,r in sorted(chosen.items()): w.writerow([gene,r['sequence_id'],mode])

p=argparse.ArgumentParser(); p.add_argument('--proteins',required=True); p.add_argument('--manifest',required=True); p.add_argument('--selection-table'); p.add_argument('--representatives',required=True); p.add_argument('--output-manifest',required=True); p.add_argument('--report',required=True)
if __name__=='__main__': main(p.parse_args())
