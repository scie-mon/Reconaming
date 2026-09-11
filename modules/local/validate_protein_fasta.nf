process VALIDATE_PROTEIN_FASTA {
    tag { proteins.baseName }
    input:
    tuple path(proteins), path(sequence_species)
    output:
    path 'representative_proteins.faa', emit: proteins
    path 'representative_manifest.tsv', emit: manifest
    path 'protein_validation_report.tsv', emit: report
    script:
    """
    python3 ${projectDir}/bin/validate_protein_fasta.py \\
      --proteins ${proteins} --sequence-species ${sequence_species} \\
      --output-proteins representative_proteins.faa \\
      --output-manifest representative_manifest.tsv \\
      --report protein_validation_report.tsv
    """
}
