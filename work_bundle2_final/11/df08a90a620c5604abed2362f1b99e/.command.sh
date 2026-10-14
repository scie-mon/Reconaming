#!/bin/bash -euo pipefail
mkdir BUSCO_results
for result in Species_2 Species_3 Root_species Species_1; do cp -aL "$result" "BUSCO_results/$(basename "$result")"; done
BUSCO_phylogenomics.py -i BUSCO_results -o busco_phylogenomics --gene_tree_program iqtree --min_species_gene_tree 3 -t 64
