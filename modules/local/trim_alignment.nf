process TRIM_ALIGNMENT {
    label 'gene_tree'
    tag 'ClipKIT smart-gap trimming'

    input:
    path alignment

    output:
    path 'representative_proteins.clipkit.faa', emit: trimmed_alignment

    script:
    """
clipkit '${alignment}' -m smart-gap -o representative_proteins.clipkit.faa
"""
}
