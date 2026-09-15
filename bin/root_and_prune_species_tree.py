#!/usr/bin/env python3
import argparse
import csv
from ete3 import Tree

parser = argparse.ArgumentParser()
parser.add_argument('--tree', required=True)
parser.add_argument('--species-tree-inputs', required=True)
parser.add_argument('--output', required=True)
args = parser.parse_args()

with open(args.species_tree_inputs, newline='') as handle:
    rows = list(csv.DictReader(handle, delimiter='\t'))
focal = {row['species_id'].strip() for row in rows if row['role'].strip() == 'focal'}
outgroups = [row['species_id'].strip() for row in rows if row['role'].strip() == 'outgroup']
if len(outgroups) != 1:
    raise SystemExit('Expected exactly one outgroup in selected inputs.')
outgroup = outgroups[0]

tree = Tree(args.tree, format=1)
leaves = set(tree.get_leaf_names())
expected = focal | {outgroup}
missing = expected - leaves
extra = leaves - expected
if missing:
    raise SystemExit('Species absent from inferred tree: ' + ', '.join(sorted(missing)))
if extra:
    raise SystemExit('Unexpected species in inferred tree: ' + ', '.join(sorted(extra)))

tree.set_outgroup(tree & outgroup)
tree.prune(sorted(focal), preserve_branch_length=True)
if set(tree.get_leaf_names()) != focal:
    raise SystemExit('Pruned tree does not contain exactly the required focal species.')
if len(focal) < 2:
    raise SystemExit('At least two focal species are required after outgroup pruning.')
tree.write(outfile=args.output, format=1, format_root_node=True)
