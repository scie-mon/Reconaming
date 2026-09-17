process ALIGN_PROTEINS {
    label 'gene_tree'
    tag 'FAMSA protein alignment'

    input:
    path representative_proteins
    val famsa_args

    output:
    path 'representative_proteins.famsa.faa', emit: alignment

    script:
    def args = famsa_args ?: ''
    """
famsa ${args} '${representative_proteins}' representative_proteins.famsa.faa
"""
}
