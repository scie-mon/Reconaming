#!/bin/bash -euo pipefail
python run_reconaming.py \
    --core-script reconaming_core.py \
    --tree reconciled_gene_tree.nhx \
    --outgroup-genes 'outgroup_genes.txt' \
     \
    --revision-tag 'bundle2_final' \
    --options reconaming_opt.ctl
