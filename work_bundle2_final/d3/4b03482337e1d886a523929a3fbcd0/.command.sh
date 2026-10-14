#!/bin/bash -euo pipefail
mkdir -p busco
busco -i CM058173.1.fa -o Species_1 --out_path busco -m genome --auto-lineage -c 16
