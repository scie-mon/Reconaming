#!/bin/bash -euo pipefail
mkdir -p busco
busco -i CM061115.fa -o Root_species --out_path busco -m genome --auto-lineage -c 16
