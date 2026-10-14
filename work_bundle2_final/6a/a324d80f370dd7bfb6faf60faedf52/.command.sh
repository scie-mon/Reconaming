#!/bin/bash -euo pipefail
python3 /home/primeline/Data/Reconaming/bin/prepare_focal_genome_map.py \
  --sequence-inputs test_species_tree_inputs.tsv \
  --output focal_genome_map.tsv
