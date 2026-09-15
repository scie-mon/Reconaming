process VALIDATE_SUPPLIED_SPECIES_TREE {
    tag { tree.baseName }
    container { params.species_tree_container ?: null }

    input:
    path tree
    path manifest
    val outgroup

    output:
    path 'validated_species_tree.nwk', emit: species_tree
    path 'species_tree_validation.tsv', emit: validation_report

    script:
    def outgroup_arg = outgroup ? "--outgroup '${outgroup}'" : ''
    """
    validate_supplied_species_tree.py \\
        --tree '${tree}' \\
        --manifest '${manifest}' \\
        ${outgroup_arg} \\
        --output validated_species_tree.nwk \\
        --report species_tree_validation.tsv
    """
}
