process RUN_RECONAMING {
    tag "reconaming:${revision_tag}"
    container { core_container }

    input:
    tuple path(reconciled_tree), path(outgroup_genes), path(revision_registry)
    path core_script
    path runner_script
    path options_file
    val revision_tag
    val core_container

    output:
    path 'named_reconaming.nhx', emit: named_tree
    path 'named_reconaming.csv', emit: name_table
    path 'revision_registry_*.csv', emit: revision_registry
    path 'temporary_partition.csv', emit: temporary_partition
    path 'reconaming_run_report.tsv', emit: report
    path 'resolved_reconaming_options.tsv', emit: resolved_options

    script:
    def outgroup_args = outgroup_genes ? "--outgroup-genes '${outgroup_genes}'" : ''
    def registry_args = revision_registry ? "--revision-registry '${revision_registry}'" : ''
    """
    python ${runner_script} \\
        --core-script ${core_script} \\
        --tree ${reconciled_tree} \\
        ${outgroup_args} \\
        ${registry_args} \\
        --revision-tag '${revision_tag}' \\
        --options ${options_file}
    """
}
