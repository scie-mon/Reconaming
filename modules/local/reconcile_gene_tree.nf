process RECONCILE_GENE_TREE {
    label 'generax'
    tag 'GeneRax reconciliation'

    cpus { generax_cpus }
    container { generax_container }

    input:
    tuple path(gene_tree), path(iqtree_report), path(trimmed_alignment), path(species_tree), path(gene_to_species), path(representative_manifest), path(outgroup_genes)
    path prepare_generax_inputs
    path normalise_generax_output
    val generax_cpus
    val generax_rec_model
    val generax_strategy
    val generax_seed
    val generax_args
    val generax_container

    output:
    path 'reconciled_gene_tree.nhx', emit: reconciled_tree
    path 'generax_raw_reconciled_tree.nhx', emit: raw_reconciled_tree
    path 'generax_gene_tree.nwk', emit: generax_gene_tree
    path 'generax_alignment.faa', emit: generax_alignment
    path 'families.txt', emit: families
    path 'generax_mapping.tsv', emit: mapping
    path 'reconciliation_summary.tsv', emit: summary
    path 'generax.log', emit: log

    script:
    def args = generax_args ?: ''
    """
    set -euo pipefail

    python3 ${prepare_generax_inputs} \\
        --gene-tree ${gene_tree} \\
        --iqtree-report ${iqtree_report} \\
        --alignment ${trimmed_alignment} \\
        --species-tree ${species_tree} \\
        --gene-to-species ${gene_to_species} \\
        --representative-manifest ${representative_manifest} \\
        --outgroup-genes ${outgroup_genes} \\
        --families-out families.txt \\
        --mapping-out generax_mapping.tsv \\
        --generax-tree-out generax_gene_tree.nwk \\
        --generax-alignment-out generax_alignment.faa \\
        --summary-out reconciliation_input_summary.tsv

    mpiexec --oversubscribe -np ${task.cpus} generax \\
        --families families.txt \\
        --species-tree ${species_tree} \\
        --rec-model '${generax_rec_model}' \\
        --strategy '${generax_strategy}' \\
        --seed '${generax_seed}' \\
        ${args} \\
        > generax.log 2>&1

    python3 ${normalise_generax_output} \\
        --generax-results . \\
        --expected-gene-tree generax_gene_tree.nwk \\
        --input-summary reconciliation_input_summary.tsv \\
        --raw-out generax_raw_reconciled_tree.nhx \\
        --normalised-out reconciled_gene_tree.nhx \\
        --summary-out reconciliation_summary.tsv
    """
}
