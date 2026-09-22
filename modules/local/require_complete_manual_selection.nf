process REQUIRE_COMPLETE_MANUAL_SELECTION {
    tag 'manual isoform selection gate'
    debug true

    input:
    path selection_status
    val template_destination

    output:
    path 'manual_selection_ready.txt', emit: ready

    script:
    """
    status=\$(awk -F '\\t' 'NR == 2 { print \$1 }' ${selection_status})
    if [[ "\$status" != 'complete' ]]; then
        cat <<'MESSAGE'
Manual isoform selection is incomplete.

A template and validation report were written to:
${template_destination}/selection_template.tsv
${template_destination}/selection_validation_report.tsv

For each multi-isoform gene, enter one gene_id<TAB>transcript_id selection.
Rows with an empty gene_id are ignored and may remain as candidate notes.
Then rerun with the revised table using --selection_table and -resume.
MESSAGE
    fi
    touch manual_selection_ready.txt
    """
}
