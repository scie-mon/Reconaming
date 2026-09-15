process PREPARE_FOCAL_GENOME_MAP {
    tag 'focal genome inputs'

    input:
    path species_inputs

    output:
    path 'focal_genome_map.tsv', emit: focal_genome_map

    script:
    """
    python3 ${projectDir}/bin/prepare_focal_genome_map.py \\
      --species-inputs ${species_inputs} \\
      --output focal_genome_map.tsv
    """
}
