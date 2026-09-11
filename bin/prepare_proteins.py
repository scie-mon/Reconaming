#!/usr/bin/env python3
import argparse, csv
from collections import defaultdict
COMP=str.maketrans('ACGTNacgtn','TGCANtgcan')
CODONS={'TTT':'F','TTC':'F','TTA':'L','TTG':'L','TCT':'S','TCC':'S','TCA':'S','TCG':'S','TAT':'Y','TAC':'Y','TAA':'*','TAG':'*','TGT':'C','TGC':'C','TGA':'*','TGG':'W','CTT':'L','CTC':'L','CTA':'L','CTG':'L','CCT':'P','CCC':'P','CCA':'P','CCG':'P','CAT':'H','CAC':'H','CAA':'Q','CAG':'Q','CGT':'R','CGC':'R','CGA':'R','CGG':'R','ATT':'I','ATC':'I','ATA':'I','ATG':'M','ACT':'T','ACC':'T','ACA':'T','ACG':'T','AAT':'N','AAC':'N','AAA':'K','AAG':'K','AGT':'S','AGC':'S','AGA':'R','AGG':'R','GTT':'V','GTC':'V','GTA':'V','GTG':'V','GCT':'A','GCC':'A','GCA':'A','GCG':'A','GAT':'D','GAC':'D','GAA':'E','GAG':'E','GGT':'G','GGC':'G','GGA':'G','GGG':'G'}
def attr(s): return {x.split('=',1)[0]:x.split('=',1)[1] for x in s.split(';') if '=' in x}
def fasta(p):
 d={}; h=None; s=[]
 for x in open(p):
  if x.startswith('>'):
   if h:d[h]=''.join(s).upper()
   h=x[1:].split()[0];s=[]
  else:s.append(x.strip())
 if h:d[h]=''.join(s).upper()
 return d
def rc(s):return s.translate(COMP)[::-1]
def trans(s):return ''.join(CODONS.get(s[i:i+3],'X') for i in range(0,len(s),3))
p=argparse.ArgumentParser();p.add_argument('--genome',required=True);p.add_argument('--gff',required=True);p.add_argument('--species-id',required=True);p.add_argument('--proteins',required=True);p.add_argument('--manifest',required=True);p.add_argument('--report',required=True);a=p.parse_args()
gen=fasta(a.genome); genes={}; tx={}; cds=defaultdict(list)
for line in open(a.gff):
 if not line.strip() or line.startswith('#'):continue
 f=line.rstrip().split('\t')
 if len(f)!=9:continue
 c,_,t,st,en,_,strand,phase,ats=f; d=attr(ats); ident=d.get('ID')
 if t=='gene' and ident:genes[ident]=(c,strand)
 elif t in ('mRNA','transcript') and ident and d.get('Parent'):tx[ident]=(d['Parent'],c,strand,d.get('internal_id',''))
 elif t=='CDS' and d.get('Parent'):cds[d['Parent']].append((int(st),int(en),phase,c,strand))
rows=[];report=[]
with open(a.proteins,'w') as o:
 for tid,(gid,c,strand,iid) in sorted(tx.items()):
  ps=cds.get(tid,[]); reason=''
  if not ps or c not in gen:reason='missing_CDS_or_contig'
  else:
   ps.sort(key=lambda z:z[0],reverse=strand=='-'); seq=''
   for st,en,phase,cc,ss in ps:
    if cc!=c or ss!=strand or phase not in ('0','1','2') or en>len(gen[c]):reason='invalid_CDS';break
    expected=(3-len(seq)%3)%3
    if int(phase)!=expected:reason='inconsistent_CDS_phase';break
    frag=gen[c][st-1:en];seq+=rc(frag) if strand=='-' else frag
   if not reason and len(seq)%3:reason='incomplete_CDS'
   protein=trans(seq) if not reason else ''
   if not reason and '*' in protein[:-1]:reason='internal_stop_codon'
  if reason:report.append([tid,gid,'excluded',reason]);continue
  protein=protein.rstrip('*');o.write(f'>{tid}\n{protein}\n');rows.append([tid,gid,tid,c,a.species_id,len(protein),iid]);report.append([tid,gid,'accepted',''])
with open(a.manifest,'w',newline='') as f:
 w=csv.writer(f,delimiter='\t');w.writerow(['sequence_id','gene_id','transcript_id','contig_id','species_id','protein_length','irnotator_internal_id']);w.writerows(rows)
with open(a.report,'w',newline='') as f:
 w=csv.writer(f,delimiter='\t');w.writerow(['transcript_id','gene_id','status','reason']);w.writerows(report)
