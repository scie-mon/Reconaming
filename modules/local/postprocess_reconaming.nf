process POSTPROCESS_RECONAMING_ANNOTATION {
    tag "postprocess:${revision_tag}"
    container { core_container }
    publishDir { "${outdir}/reconaming/${revision_tag}" }, mode: 'copy', overwrite: true
    input:
    tuple path(named_tree), path(name_table), path(revision_registry), path(temporary_partition), path(core_report), path(gff_files)
    path postprocess_script
    val gff_output
    val gene_id_attribute
    val revision_tag
    val outdir
    val core_container
    output:
    path 'reconaming_results_manifest.tsv', emit: manifest
    path 'reconaming_validation_report.tsv', emit: validation
    path 'reconaming_*.gff3', optional: true, emit: gff
    path 'named_reconaming.nhx', includeInputs: true, emit: named_tree
    path 'named_reconaming.csv', includeInputs: true, emit: name_table
    path 'revision_registry_*.csv', includeInputs: true, emit: revision_registry
    path 'temporary_partition.csv', includeInputs: true, emit: temporary_partition
    path 'reconaming_run_report.tsv', includeInputs: true, emit: core_report
    script:
    """
    python ${postprocess_script} --named-tree ${named_tree} --name-table ${name_table} \\
      --revision-registry ${revision_registry} --temporary-partition ${temporary_partition} \\
      --core-report ${core_report} --gff-files ${gff_files} --gff-output ${gff_output} \\
      --gene-id-attribute '${gene_id_attribute}' --revision-tag '${revision_tag}' --input-mode annotation
    """
}

process POSTPROCESS_RECONAMING_PROTEIN {
    tag "postprocess:${revision_tag}"
    container { core_container }
    publishDir { "${outdir}/reconaming/${revision_tag}" }, mode: 'copy', overwrite: true
    input:
    tuple path(named_tree), path(name_table), path(revision_registry), path(temporary_partition), path(core_report)
    path postprocess_script
    val gff_output
    val revision_tag
    val outdir
    val core_container
    output:
    path 'reconaming_results_manifest.tsv', emit: manifest
    path 'reconaming_validation_report.tsv', emit: validation
    path 'named_reconaming.nhx', includeInputs: true, emit: named_tree
    path 'named_reconaming.csv', includeInputs: true, emit: name_table
    path 'revision_registry_*.csv', includeInputs: true, emit: revision_registry
    path 'temporary_partition.csv', includeInputs: true, emit: temporary_partition
    path 'reconaming_run_report.tsv', includeInputs: true, emit: core_report
    script:
    """
    python ${postprocess_script} --named-tree ${named_tree} --name-table ${name_table} \\
      --revision-registry ${revision_registry} --temporary-partition ${temporary_partition} \\
      --core-report ${core_report} --gff-output ${gff_output} \\
      --revision-tag '${revision_tag}' --input-mode protein
    """
}
