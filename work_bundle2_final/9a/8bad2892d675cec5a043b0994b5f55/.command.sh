#!/bin/bash -euo pipefail
python3 select_representative_isoforms.py \
  --proteins isoform_proteins.faa \
  --manifest isoform_manifest.tsv \
   \
  --selection-table-supplied false \
  --representatives representative_proteins.faa \
  --output-manifest representative_manifest.tsv \
  --report isoform_selection_report.tsv \
  --validation-report selection_validation_report.tsv \
  --template selection_template.tsv \
  --status selection_status.tsv
