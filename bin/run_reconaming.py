#!/usr/bin/env python3
import argparse
import csv
import re
import shutil
import subprocess
import sys
from pathlib import Path

PIPELINE_OPTIONS = {
    "input_mode", "genome_fasta", "ir_gff", "protein_fasta", "sequence_mapping",
    "isoform_selection", "minimum_protein_length", "invalid_protein_policy",
    "allow_terminal_stop", "automatic_isoform_tiebreak", "outdir"
}
CORE_OPTIONS = {
    "prefix", "threshold", "aliases", "ignore_sco", "minimal_id",
    "include_outgroup", "sort_tree", "root_by", "pgignore_nodes",
    "max_pg_size", "allow_revive_pg", "force_pg_root", "show_single_para",
    "force_boundary_true", "force_boundary_false"
}


def fail(message):
    raise SystemExit(f"ERROR: {message}")


def read_ctl(path):
    values = {}
    for number, raw in enumerate(Path(path).read_text(encoding="utf-8").splitlines(), 1):
        line = raw.strip()
        if not line or line.startswith("#"):
            continue
        if "=" not in line:
            fail(f"{path}:{number}: expected key=value")
        key, value = (part.strip() for part in line.split("=", 1))
        if not re.fullmatch(r"[a-z][a-z0-9_]*", key):
            fail(f"{path}:{number}: invalid option name '{key}'")
        if key in values:
            fail(f"{path}:{number}: duplicate option '{key}'")
        values[key] = value
    unknown = set(values) - PIPELINE_OPTIONS - CORE_OPTIONS
    if unknown:
        fail("unknown reconaming_opt.ctl option(s): " + ", ".join(sorted(unknown)))
    return values


def bool_value(options, key):
    value = options.get(key, "").lower()
    if not value:
        return None
    if value in {"true", "yes", "1"}:
        return True
    if value in {"false", "no", "0"}:
        return False
    fail(f"{key} must be true or false")


def outgroup_tips(path):
    return [line.strip() for line in Path(path).read_text(encoding="utf-8").splitlines()
            if line.strip() and not line.lstrip().startswith("#")]


def option_args(options, outgroup):
    args = []
    scalar = {
        "prefix": "--prefix", "threshold": "--threshold", "aliases": "--aliases",
        "minimal_id": "--minimal-id", "pgignore_nodes": "--PGignore-nodes",
        "max_pg_size": "--max-pg-size", "force_pg_root": "--force-pg-root",
        "force_boundary_true": "--force-boundary-true",
        "force_boundary_false": "--force-boundary-false",
    }
    for key, flag in scalar.items():
        if options.get(key, ""):
            args.extend([flag, options[key]])
    roots = options.get("root_by", "").split() or outgroup
    if roots:
        args.extend(["--root-by", *roots])
    if bool_value(options, "ignore_sco") is False:
        args.append("--no-ignore-sco")
    if bool_value(options, "include_outgroup") is True:
        args.append("--include-outgroup")
    if bool_value(options, "sort_tree") is False:
        args.append("--no-sort-tree")
    if bool_value(options, "allow_revive_pg") is False:
        args.append("--no-allow-revive-pg")
    if bool_value(options, "show_single_para") is True:
        args.append("--show-single-para")
    return args, roots


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--core-script", required=True)
    parser.add_argument("--tree", required=True)
    parser.add_argument("--outgroup-genes", required=True)
    parser.add_argument("--revision-registry", required=True)
    parser.add_argument("--revision-tag", required=True)
    parser.add_argument("--options", required=True)
    args = parser.parse_args()
    if not re.fullmatch(r"[A-Za-z0-9][A-Za-z0-9_.-]*", args.revision_tag):
        fail("revision tag may contain only letters, numbers, dots, underscores, and hyphens")

    options = read_ctl(args.options)
    core_args, roots = option_args(options, outgroup_tips(args.outgroup_genes))
    command = [sys.executable, args.core_script, "--infile", args.tree,
               "--outfile", "named_reconaming.nhx", "--revlog", args.revision_registry,
               "--revtag", args.revision_tag, "--temp-partition-out", "temporary_partition.csv",
               *core_args]
    completed = subprocess.run(command, text=True, capture_output=True)
    Path("reconaming_core.stdout.log").write_text(completed.stdout, encoding="utf-8")
    Path("reconaming_core.stderr.log").write_text(completed.stderr, encoding="utf-8")
    if completed.returncode:
        sys.stderr.write(completed.stderr)
        raise SystemExit(completed.returncode)

    produced_registry = Path(f"revlog_{args.revision_tag}.csv")
    if not produced_registry.exists() or not Path("named_reconaming.csv").exists():
        fail("Reconaming completed without producing all required outputs")
    shutil.move(str(produced_registry), f"revision_registry_{args.revision_tag}.csv")

    with open("resolved_reconaming_options.tsv", "w", newline="", encoding="utf-8") as handle:
        writer = csv.writer(handle, delimiter="\t")
        writer.writerow(["option", "value"])
        for key in sorted(options):
            writer.writerow([key, options[key]])
        writer.writerow(["effective_root_by", " ".join(roots)])
        writer.writerow(["revision_tag", args.revision_tag])
    with open("reconaming_run_report.tsv", "w", newline="", encoding="utf-8") as handle:
        writer = csv.writer(handle, delimiter="\t")
        writer.writerow(["metric", "value"])
        writer.writerows([["revision_tag", args.revision_tag], ["root_by_count", len(roots)],
                          ["named_tree", "named_reconaming.nhx"],
                          ["name_table", "named_reconaming.csv"],
                          ["revision_registry", f"revision_registry_{args.revision_tag}.csv"]])


if __name__ == "__main__":
    main()
