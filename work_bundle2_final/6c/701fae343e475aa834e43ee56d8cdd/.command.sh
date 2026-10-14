#!/bin/bash -euo pipefail
iqtree -s SUPERMATRIX.phylip -p SUPERMATRIX.partitions.nex -B 1000 -T 64 --prefix species_tree_unrooted
