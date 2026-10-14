#!/bin/bash -euo pipefail
mkdir -p busco
busco -i CM023194.1.fa -o Species_2 --out_path busco -m genome --auto-lineage -c 16
