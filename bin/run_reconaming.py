#!/usr/bin/env python3
import argparse, csv, re, shutil, subprocess, sys
from pathlib import Path

PIPELINE_OPTIONS = {"input_mode","genome_fasta","ir_gff","protein_fasta","sequence_mapping","isoform_selection","minimum_protein_length","invalid_protein_policy","allow_terminal_stop","automatic_isoform_tiebreak","outdir","gff_output"}
CORE_OPTIONS = {"prefix","threshold","aliases","ignore_sco","minimal_id","include_outgroup","sort_tree","root_by","pgignore_nodes","max_pg_size","allow_revive_pg","force_pg_root","show_single_para","force_boundary_true","force_boundary_false"}

def fail(message): raise SystemExit(f"ERROR: {message}")
def read_ctl(path):
    values = {}
    for number, raw in enumerate(Path(path).read_text(encoding="utf-8").splitlines(), 1):
        line = raw.strip()
        if not line or line.startswith("#"): continue
        if "=" not in line: fail(f"{path}:{number}: expected key=value")
        key, value = (x.strip() for x in line.split("=", 1))
        if not re.fullmatch(r"[a-z][a-z0-9_]*", key): fail(f"{path}:{number}: invalid option '{key}'")
        if key in values: fail(f"{path}:{number}: duplicate option '{key}'")
        values[key] = value
    unknown = set(values) - PIPELINE_OPTIONS - CORE_OPTIONS
    if unknown: fail("unknown reconaming_opt.ctl option(s): " + ", ".join(sorted(unknown)))
    return values
def boolean(options, key):
    value = options.get(key, "").lower()
    if not value: return None
    if value in {"true","yes","1"}: return True
    if value in {"false","no","0"}: return False
    fail(f"{key} must be true or false")
def outgroup_tips(path):
    return [x.strip() for x in Path(path).read_text(encoding="utf-8").splitlines() if x.strip() and not x.lstrip().startswith("#")]
def core_args(options, outgroups):
    args=[]
    for key, flag in {"prefix":"--prefix","threshold":"--threshold","aliases":"--aliases","minimal_id":"--minimal-id","pgignore_nodes":"--PGignore-nodes","max_pg_size":"--max-pg-size","force_pg_root":"--force-pg-root","force_boundary_true":"--force-boundary-true","force_boundary_false":"--force-boundary-false"}.items():
        if options.get(key, ""): args += [flag, options[key]]
    roots = options.get("root_by", "").split() or outgroups
    if roots: args += ["--root-by", *roots]
    if boolean(options,"ignore_sco") is False: args += ["--no-ignore-sco"]
    if boolean(options,"include_outgroup") is True: args += ["--include-outgroup"]
    if boolean(options,"sort_tree") is False: args += ["--no-sort-tree"]
    if boolean(options,"allow_revive_pg") is False: args += ["--no-allow-revive-pg"]
    if boolean(options,"show_single_para") is True: args += ["--show-single-para"]
    return args, roots
def main():
    p=argparse.ArgumentParser()
    for key in ("core_script","tree","outgroup_genes","revision_registry","revision_tag","options"): p.add_argument("--"+key.replace("_","-"), required=True)
    a=p.parse_args()
    if not re.fullmatch(r"[A-Za-z0-9][A-Za-z0-9_.-]*", a.revision_tag): fail("revision tag contains invalid characters")
    options=read_ctl(a.options); args, roots=core_args(options,outgroup_tips(a.outgroup_genes))
    command=[sys.executable,a.core_script,"--infile",a.tree,"--outfile","named_reconaming.nhx","--revlog",a.revision_registry,"--revtag",a.revision_tag,"--temp-partition-out","temporary_partition.csv",*args]
    result=subprocess.run(command,text=True,capture_output=True)
    Path("reconaming_core.stdout.log").write_text(result.stdout,encoding="utf-8"); Path("reconaming_core.stderr.log").write_text(result.stderr,encoding="utf-8")
    if result.returncode: sys.stderr.write(result.stderr); raise SystemExit(result.returncode)
    produced=Path(f"revlog_{a.revision_tag}.csv")
    if not produced.exists() or not Path("named_reconaming.csv").exists(): fail("Reconaming did not produce all required outputs")
    shutil.move(str(produced),f"revision_registry_{a.revision_tag}.csv")
    with open("resolved_reconaming_options.tsv","w",newline="",encoding="utf-8") as h:
        w=csv.writer(h,delimiter="\t"); w.writerow(["option","value"])
        for k in sorted(options): w.writerow([k,options[k]])
        w.writerow(["effective_root_by"," ".join(roots)]); w.writerow(["revision_tag",a.revision_tag])
    with open("reconaming_run_report.tsv","w",newline="",encoding="utf-8") as h:
        w=csv.writer(h,delimiter="\t"); w.writerow(["metric","value"]); w.writerows([["revision_tag",a.revision_tag],["root_by_count",len(roots)],["named_tree","named_reconaming.nhx"],["name_table","named_reconaming.csv"],["revision_registry",f"revision_registry_{a.revision_tag}.csv"]])
if __name__ == "__main__": main()
