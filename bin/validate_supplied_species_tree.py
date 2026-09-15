#!/usr/bin/env python3
import argparse
import csv
import shutil
from pathlib import Path

from ete3 import Tree


def read_focal_species(manifest_path):
    with open(manifest_path, newline="") as handle:
        sample = handle.read(4096)
        handle.seek(0)
        try:
            dialect = csv.Sniffer().sniff(sample, delimiters="\t,")
        except csv.Error:
            dialect = csv.excel_tab
        reader = csv.DictReader(handle, dialect=dialect)
        if not reader.fieldnames or "species_id" not in reader.fieldnames:
            raise SystemExit("Representative manifest must contain a species_id column.")
        focal = set()
        for row in reader:
            species_id = (row.get("species_id") or "").strip()
            role = (row.get("role") or "focal").strip().lower()
            if species_id and role != "outgroup":
                focal.add(species_id)
    if len(focal) < 2:
        raise SystemExit("Representative manifest must contain at least two focal species.")
    return focal


def require_binary_rooted(tree):
    if tree.is_leaf() or len(tree.children) != 2:
        raise SystemExit("Species tree must be rooted and strictly bifurcating.")
    for node in tree.traverse():
        if not node.is_leaf() and node is not tree and len(node.children) != 2:
            raise SystemExit("Species tree contains a polytomy or unary internal node.")


def validate_leaf_names(tree):
    names = tree.get_leaf_names()
    duplicates = sorted({name for name in names if names.count(name) > 1})
    if duplicates:
        raise SystemExit("Species tree contains duplicate leaf labels: " + ", ".join(duplicates))
    return set(names)


def require_equal(observed, expected, label):
    missing = sorted(expected - observed)
    extra = sorted(observed - expected)
    if missing or extra:
        parts = []
        if missing:
            parts.append("missing " + ", ".join(missing))
        if extra:
            parts.append("unexpected " + ", ".join(extra))
        raise SystemExit(f"{label} does not match focal species IDs: " + "; ".join(parts))


parser = argparse.ArgumentParser(description="Validate and optionally root/prune a supplied species tree.")
parser.add_argument("--tree", required=True)
parser.add_argument("--manifest", required=True)
parser.add_argument("--outgroup", default="")
parser.add_argument("--output", required=True)
parser.add_argument("--report", required=True)
args = parser.parse_args()

focal = read_focal_species(args.manifest)
try:
    tree = Tree(args.tree, format=1)
except Exception as exc:
    raise SystemExit(f"Unable to parse supplied species tree as Newick: {exc}")

leaves = validate_leaf_names(tree)
outgroup = args.outgroup.strip()

if outgroup:
    if outgroup in focal:
        raise SystemExit("species_tree_outgroup must not be a focal species_id.")
    require_equal(leaves, focal | {outgroup}, "Supplied tree leaves before outgroup pruning")
    tree.set_outgroup(tree & outgroup)
    tree.prune(sorted(focal), preserve_branch_length=True)
    validate_leaf_names(tree)
    require_equal(set(tree.get_leaf_names()), focal, "Supplied tree leaves after outgroup pruning")
    require_binary_rooted(tree)
    tree.write(outfile=args.output, format=1, format_root_node=True)
    action = "rooted_on_outgroup_and_pruned"
else:
    require_equal(leaves, focal, "Supplied tree leaves")
    require_binary_rooted(tree)
    shutil.copyfile(args.tree, args.output)
    action = "validated_without_modification"

with open(args.report, "w", newline="") as handle:
    writer = csv.writer(handle, delimiter="\t")
    writer.writerow(["status", "action", "focal_species_count", "outgroup"])
    writer.writerow(["valid", action, len(focal), outgroup])
