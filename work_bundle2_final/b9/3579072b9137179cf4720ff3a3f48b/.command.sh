#!/bin/bash -euo pipefail
status=$(awk -F '\t' 'NR == 2 { print $1 }' selection_status.tsv)
    if [[ "$status" != 'complete' ]]; then
        cat <<'MESSAGE'
Manual isoform selection is incomplete.

A template and validation report were written to:
results_bundle2_final/manual_selection/selection_template.tsv
results_bundle2_final/manual_selection/selection_validation_report.tsv

For each multi-isoform gene, enter one gene_id<TAB>transcript_id selection.
Rows with an empty gene_id are ignored and may remain as candidate notes.
Then rerun with the revised table using --selection_table and -resume.
MESSAGE
    fi
    touch manual_selection_ready.txt
