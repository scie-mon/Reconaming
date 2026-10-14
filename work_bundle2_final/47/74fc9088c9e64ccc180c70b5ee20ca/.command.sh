#!/bin/bash -euo pipefail
python postprocess_reconaming.py --named-tree named_reconaming.nhx --name-table named_reconaming.csv \
  --revision-registry revision_registry_bundle2_final.csv --temporary-partition temporary_partition.csv \
  --core-report reconaming_run_report.tsv --gff-files merged.rev0.gff3 --gff-output update \
  --gene-id-attribute 'ID' --revision-tag 'bundle2_final' --input-mode annotation
