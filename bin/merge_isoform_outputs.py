#!/usr/bin/env python3
import argparse,csv
p=argparse.ArgumentParser(); p.add_argument('--fastas',nargs='+',required=True); p.add_argument('--manifests',nargs='+',required=True); p.add_argument('--proteins',required=True); p.add_argument('--manifest',required=True); a=p.parse_args()
with open(a.proteins,'w') as out:
 for path in a.fastas:
  with open(path) as f: out.write(f.read().rstrip()+'\n')
rows=[]; header=None
for path in a.manifests:
 with open(path,newline='') as f:
  r=csv.DictReader(f,delimiter='\t'); header=header or r.fieldnames; rows.extend(r)
with open(a.manifest,'w',newline='') as f:
 w=csv.DictWriter(f,fieldnames=header,delimiter='\t'); w.writeheader(); w.writerows(rows)
