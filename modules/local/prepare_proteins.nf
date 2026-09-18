process PREPARE_PROTEINS {
    tag { species_id }

    input:
    tuple val(species_id), path(genome), path(gff)
    path prepare_proteins_script

    output:
    path 'isoform_proteins.faa', emit: proteins
    path 'isoform_manifest.tsv', emit: manifest
    path 'protein_preparation_report.tsv', emit: report

    script:
    """
    python3 prepare_proteins.py \
        --genome ${genome} \
        --gff ${gff} \
        --species-id '${species_id}' \
        --proteins isoform_proteins.faa \
        --manifest isoform_manifest.tsv \
        --report protein_preparation_report.tsv
    """
}