#!/bin/bash -euo pipefail
root_and_prune_species_tree.py --tree species_tree_unrooted.treefile --species-tree-inputs selected_species_inputs.tsv --output species_tree.rooted_pruned.nwk
