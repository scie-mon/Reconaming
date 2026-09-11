process VALIDATE_INPUTS {
    tag 'input contract'

    input:
    path genome_map
    path gff
    path protein_fasta
    path sequence_species

    output:
    path 'resolved_settings.tsv', emit: settings
    path 'input_validation_report.tsv', emit: report

    script:
    def genomeArg = genome_map ? "--genome-map ${genome_map}" : ''
    def gffArg = gff ? "--gff ${gff}" : ''
    def proteinArg = protein_fasta ? "--protein-fasta ${protein_fasta}" : ''
    def mappingArg = sequence_species ? "--sequence-mapping ${sequence_species}" : ''
    """
    python3 ${projectDir}/bin/validate_inputs.py \\
        ${genomeArg} ${gffArg} ${proteinArg} ${mappingArg} \\
        --settings resolved_settings.tsv \\
        --report input_validation_report.tsv
    """
}
