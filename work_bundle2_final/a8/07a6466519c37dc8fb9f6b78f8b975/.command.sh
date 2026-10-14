#!/bin/bash -euo pipefail
mkdir -p busco
busco -i NC_081337.1.fa -o Species_3 --out_path busco -m genome --auto-lineage -c 16
