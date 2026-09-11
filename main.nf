#!/usr/bin/env nextflow
nextflow.enable.dsl=2

params.input_mode        = params.input_mode ?: 'annotation'
params.genome_to_species = params.genome_to_species ?: null
params.ir_gff            = params.ir_gff ?: null
params.protein_fasta     = params.protein_fasta ?: null
params.sequence_species  = params.sequence_species ?: null
params.isoform_selection = params.isoform_selection ?: null

include { BUILD_INPUT_MANIFEST }          from './modules/local/build_input_manifest'
include { PREPARE_PROTEINS }               from './modules/local/prepare_proteins'
include { MERGE_ISOFORM_OUTPUTS }          from './modules/local/merge_isoform_outputs'
include { VALIDATE_PROTEIN_FASTA }         from './modules/local/validate_protein_fasta'
include { SELECT_REPRESENTATIVE_ISOFORMS } from './modules/local/select_representative_isoforms'

workflow {
    if (params.input_mode == 'protein') {
        if (!params.protein_fasta || !params.sequence_species) error 'Protein mode requires --protein_fasta and --sequence_species.'
        VALIDATE_PROTEIN_FASTA(Channel.of(tuple(file(params.protein_fasta), file(params.sequence_species))))
        representative_proteins = VALIDATE_PROTEIN_FASTA.out.proteins
        representative_manifest = VALIDATE_PROTEIN_FASTA.out.manifest
    }
    else if (params.input_mode == 'annotation') {
        if (!params.genome_to_species || !params.ir_gff) error 'Annotation mode requires --genome_to_species and --ir_gff.'
        BUILD_INPUT_MANIFEST(Channel.fromPath(params.genome_to_species, checkIfExists: true))
        gffs = Channel.fromPath(params.ir_gff, checkIfExists: true).collect()
        focal_genomes = Channel.fromPath(params.genome_to_species)
            .splitCsv(header: true, sep: '\t')
            .filter { row -> row.role == 'focal' }
            .map { row -> tuple(row.species_id, file(row.genome_file)) }
        PREPARE_PROTEINS(focal_genomes.combine(gffs).map { s, genome, files -> tuple(s, genome, files) })
        MERGE_ISOFORM_OUTPUTS(PREPARE_PROTEINS.out.proteins.collect(), PREPARE_PROTEINS.out.manifest.collect())
        selection = params.isoform_selection ? file(params.isoform_selection).toAbsolutePath().toString() : ''
        SELECT_REPRESENTATIVE_ISOFORMS(MERGE_ISOFORM_OUTPUTS.out.proteins, MERGE_ISOFORM_OUTPUTS.out.manifest, selection)
        representative_proteins = SELECT_REPRESENTATIVE_ISOFORMS.out.proteins
        representative_manifest = SELECT_REPRESENTATIVE_ISOFORMS.out.manifest
    }
    else error "Unknown --input_mode '${params.input_mode}'; use annotation or protein."

    // Jobs 05–18 consume representative_proteins and representative_manifest.
}
