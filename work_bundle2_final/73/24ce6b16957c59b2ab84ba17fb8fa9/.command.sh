#!/bin/bash -euo pipefail
python3 /home/primeline/Data/Reconaming/bin/build_input_manifest.py \
  --genome-map focal_genome_map.tsv \
  --base-dir /home/primeline/Data/Reconaming \
  --concatenated-fasta combined_genomes.fna \
  --contig-map contig_to_species.tsv
