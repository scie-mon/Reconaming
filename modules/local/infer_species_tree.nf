process SELECT_SPECIES_TREE_INPUTS {
    tag 'select BUSCO species'
    input:
    path representative_manifest
    path species_inputs
    output:
    path 'selected_species_inputs.tsv', emit: selected_inputs
    script:
    """
    select_species_tree_inputs.py --representative-manifest ${representative_manifest} --species-tree-inputs ${species_inputs} --busco-mode ${params.busco_mode} --output selected_species_inputs.tsv
    """
}

process RUN_BUSCO {
    tag { species_id }
    cpus { params.busco_threads as int }
    container { params.species_tree_container ?: null }
    input:
    tuple val(species_id), path(input_file), val(input_type)
    output:
    tuple val(species_id), path("busco/${species_id}"), emit: busco_dir
    script:
    if (!params.busco_lineage && !params.busco_auto_lineage) error 'Set --busco_lineage or enable --busco_auto_lineage.'
    def lineage_arg = params.busco_lineage ? "-l ${params.busco_lineage}" : '--auto-lineage'
    def mode = input_type == 'genome' ? 'genome' : 'proteins'
    """
    mkdir -p busco
    busco -i ${input_file} -o ${species_id} --out_path busco -m ${mode} ${lineage_arg} -c ${task.cpus}
    """
}

process ASSESS_BUSCO_COMPLETENESS {
    tag 'BUSCO completeness'
    container { params.species_tree_container ?: null }
    input:
    path busco_dirs
    output:
    path 'busco_completeness.tsv', emit: completeness
    script:
    """
    assess_busco_completeness.py --busco-dirs ${busco_dirs} --warning-threshold ${params.busco_warn_complete} --output busco_completeness.tsv
    """
}

process RUN_BUSCO_PHYLOGENOMICS {
    tag 'BUSCO phylogenomics'
    cpus { params.species_tree_threads as int }
    container { params.species_tree_container ?: null }
    input:
    path busco_dirs
    output:
    path 'busco_phylogenomics/supermatrix/SUPERMATRIX.phylip', emit: supermatrix
    path 'busco_phylogenomics/supermatrix/SUPERMATRIX.partitions.nex', emit: partitions
    script:
    """
    mkdir BUSCO_results
    for result in ${busco_dirs}; do
        cp -aL "\$result" "BUSCO_results/\$(basename "\$result")"
    done
    BUSCO_phylogenomics.py -i BUSCO_results -o busco_phylogenomics --gene_tree_program iqtree --min_species_gene_tree ${params.species_tree_min_taxa} -t ${task.cpus}
    """
}

process INFER_SUPERMATRIX_TREE {
    tag 'IQ-TREE species tree'
    cpus { params.species_tree_threads as int }
    container { params.species_tree_container ?: null }
    input:
    path supermatrix
    path partitions
    output:
    path 'species_tree_unrooted.treefile', emit: tree
    script:
    def bootstrap_arg = (params.species_tree_bootstrap as int) > 0 ? "-B ${params.species_tree_bootstrap}" : ''
    """
    iqtree -s ${supermatrix} -p ${partitions} ${bootstrap_arg} -T ${task.cpus} --prefix species_tree_unrooted
    """
}

process ROOT_AND_PRUNE_SPECIES_TREE {
    tag 'root and prune auxiliary outgroup'
    container { params.species_tree_container ?: null }
    input:
    path unrooted_tree
    path selected_inputs
    output:
    path 'species_tree.rooted_pruned.nwk', emit: species_tree
    script:
    """
    root_and_prune_species_tree.py --tree ${unrooted_tree} --species-tree-inputs ${selected_inputs} --output species_tree.rooted_pruned.nwk
    """
}

workflow INFER_SPECIES_TREE {
    take:
    representative_manifest
    species_inputs
    main:
    SELECT_SPECIES_TREE_INPUTS(representative_manifest, species_inputs)
    busco_rows = SELECT_SPECIES_TREE_INPUTS.out.selected_inputs
        .splitCsv(header: true, sep: '\t')
        .map { row -> tuple(row.species_id, file(row.input_file), row.input_type) }
    RUN_BUSCO(busco_rows)
    busco_dirs = RUN_BUSCO.out.busco_dir.map { species_id, directory -> directory }.collect()
    ASSESS_BUSCO_COMPLETENESS(busco_dirs)
    RUN_BUSCO_PHYLOGENOMICS(busco_dirs)
    INFER_SUPERMATRIX_TREE(RUN_BUSCO_PHYLOGENOMICS.out.supermatrix, RUN_BUSCO_PHYLOGENOMICS.out.partitions)
    ROOT_AND_PRUNE_SPECIES_TREE(INFER_SUPERMATRIX_TREE.out.tree, SELECT_SPECIES_TREE_INPUTS.out.selected_inputs)
    emit:
    species_tree = ROOT_AND_PRUNE_SPECIES_TREE.out.species_tree
    busco_completeness = ASSESS_BUSCO_COMPLETENESS.out.completeness
}