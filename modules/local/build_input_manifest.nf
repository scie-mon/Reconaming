process BUILD_INPUT_MANIFEST {
    tag 'multi-genome annotation inputs'
    input:
    path genome_map
    output:
    path 'combined_genomes.fna', emit: genomes
    path 'contig_to_species.tsv', emit: contigs
    script:
    """
    python3 ${projectDir}/bin/build_input_manifest.py \\
      --genome-map ${genome_map} \\
      --base-dir ${projectDir} \\
      --concatenated-fasta combined_genomes.fna \\
      --contig-map contig_to_species.tsv
    """
}
