process PREPARE_PROTEINS {
    tag { species_id }

    input:
    tuple val(species_id), path(genome), path(gff)
    path prepare_proteins_script
    val gene_id_attribute
    val transcript_id_attribute
    val protein_id_attribute

    output:
    path 'isoform_proteins.faa', emit: proteins
    path 'isoform_manifest.tsv', emit: manifest
    path 'protein_preparation_report.tsv', emit: report

    script:
    """
    python3 prepare_proteins.py \\
        --genome ${genome} \\
        --gff ${gff} \\
        --species-id '${species_id}' \\
        --gene-id-attribute '${gene_id_attribute}' \\
        --transcript-id-attribute '${transcript_id_attribute}' \\
        --protein-id-attribute '${protein_id_attribute}' \\
        --proteins isoform_proteins.faa \\
        --manifest isoform_manifest.tsv \\
        --report protein_preparation_report.tsv
    """
}
