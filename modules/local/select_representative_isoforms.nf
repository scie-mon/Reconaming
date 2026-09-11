process SELECT_REPRESENTATIVE_ISOFORMS {
    tag 'annotation isoform selection'
    input:
    path proteins
    path manifest
    val selection_table
    output:
    path 'representative_proteins.faa', emit: proteins
    path 'representative_manifest.tsv', emit: manifest
    path 'isoform_selection_report.tsv', emit: report
    script:
    def selection = selection_table ? "--selection-table ${selection_table}" : ''
    """
    python3 ${projectDir}/bin/select_representative_isoforms.py \\
      --proteins ${proteins} --manifest ${manifest} ${selection} \\
      --representatives representative_proteins.faa \\
      --output-manifest representative_manifest.tsv \\
      --report isoform_selection_report.tsv
    """
}
