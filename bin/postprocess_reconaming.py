#!/usr/bin/env python3
import argparse, csv, hashlib
from pathlib import Path

def attrs(text):
    return dict(x.split("=", 1) for x in text.split(";") if "=" in x)
def attr_text(values):
    return ";".join(f"{k}={v}" for k,v in values.items())
def names(path):
    with open(path, newline="", encoding="utf-8") as h:
        return {r["Name"]: r["geneID"] for r in csv.DictReader(h) if r.get("geneID") and r["geneID"] != "XXXXXXXXXXX"}
def load_gff(paths):
    rows=[]; headers=[]
    for path in paths:
        for line in Path(path).read_text(encoding="utf-8").splitlines():
            if line.startswith("#"): headers.append(line)
            elif line:
                fields=line.split("\t")
                if len(fields)==9: rows.append(fields)
    return headers, rows
def descendants(rows, identifiers):
    kept=set(identifiers); changed=True
    while changed:
        changed=False
        for row in rows:
            a=attrs(row[8]); ident=a.get("ID"); parents=set(a.get("Parent", "").split(","))-{''}
            if ident and parents & kept and ident not in kept: kept.add(ident); changed=True
    return kept
def render(paths, mode, mapping, gene_key, tag):
    headers, rows=load_gff(paths); output=["##gff-version 3"]+[x for x in headers if x != "##gff-version 3"]
    selected={attrs(r[8]).get("ID") for r in rows if r[2]=="gene" and attrs(r[8]).get(gene_key) in mapping}-{None}
    keep=descendants(rows,selected) if mode=="minimal" else None; count=0
    for row in rows:
        a=attrs(row[8]); locus=a.get(gene_key)
        if mode=="minimal" and a.get("ID") not in keep: continue
        if row[2]=="gene" and locus in mapping:
            old=a.get("Name")
            if mode=="minimal": a={"ID":a["ID"]}
            if old: a.setdefault("source_name",old)
            a["Name"]=mapping[locus]; a["reconaming_revision"]=tag; count+=1
        elif mode=="minimal": a={k:a[k] for k in ("ID","Parent") if a.get(k)}
        row[8]=attr_text(a); output.append("\t".join(row))
    return output,count
def digest(path):
    h=hashlib.sha256()
    with open(path,"rb") as f:
        for b in iter(lambda:f.read(1048576),b""): h.update(b)
    return h.hexdigest()
def main():
    p=argparse.ArgumentParser()
    for k in ("named_tree","name_table","revision_registry","temporary_partition","core_report","revision_tag"): p.add_argument("--"+k.replace("_","-"),required=True)
    p.add_argument("--gff-files",nargs="*"); p.add_argument("--gff-output",choices=("update","minimal","none"),required=True); p.add_argument("--gene-id-attribute",default="ID"); p.add_argument("--input-mode",choices=("annotation","protein"),required=True)
    a=p.parse_args()
    if a.input_mode=="protein" and a.gff_output!="none": raise SystemExit("ERROR: protein input mode requires gff_output=none")
    if a.gff_output!="none" and not a.gff_files: raise SystemExit("ERROR: GFF output requires source GFF input")
    mapping=names(a.name_table); outputs=[a.named_tree,a.name_table,a.revision_registry,a.temporary_partition,a.core_report]; gff_count=0
    if a.gff_output!="none":
        lines,gff_count=render(a.gff_files,a.gff_output,mapping,a.gene_id_attribute,a.revision_tag); gff=f"reconaming_{a.revision_tag}.gff3"; Path(gff).write_text("\n".join(lines)+"\n",encoding="utf-8"); outputs.append(gff)
    with open("reconaming_results_manifest.tsv","w",newline="",encoding="utf-8") as h:
        w=csv.writer(h,delimiter="\t"); w.writerow(["file","sha256"]); [w.writerow([Path(x).name,digest(x)]) for x in outputs]
    with open("reconaming_validation_report.tsv","w",newline="",encoding="utf-8") as h:
        w=csv.writer(h,delimiter="\t"); w.writerow(["metric","value"]); w.writerows([["revision_tag",a.revision_tag],["named_loci",len(mapping)],["gff_output",a.gff_output],["gff_named_genes",gff_count]])
if __name__=="__main__": main()
