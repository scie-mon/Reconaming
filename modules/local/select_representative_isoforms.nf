process SELECT_REPRESENTATIVE_ISOFORMS {
    tag 'annotation isoform selection'

    publishDir "${params.outdir}/manual_selection", mode: 'copy', saveAs: { filename ->
        filename in ['selection_template.tsv', 'selection_validation_report.tsv'] ? filename : null
    }

    input:
    path proteins
    path manifest
    path selection_table
    path selector_script
    val selection_table_supplied

    output:
    path 'representative_proteins.faa', emit: proteins, optional: true
    path 'representative_manifest.tsv', emit: manifest, optional: true
    path 'isoform_selection_report.tsv', emit: report
    path 'selection_validation_report.tsv', emit: validation_report
    path 'selection_status.tsv', emit: status
    path 'selection_template.tsv', emit: template, optional: true

    script:
    """
    python3 select_representative_isoforms.py \\
      --proteins ${proteins} \\
      --manifest ${manifest} \\
      --selection-table ${selection_table} \\
      --selection-table-supplied ${selection_table_supplied} \\
      --representatives representative_proteins.faa \\
      --output-manifest representative_manifest.tsv \\
      --report isoform_selection_report.tsv \\
      --validation-report selection_validation_report.tsv \\
      --template selection_template.tsv \\
      --status selection_status.tsv
    """
}
