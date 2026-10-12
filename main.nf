#!/usr/bin/env nextflow
nextflow.enable.dsl=2

params.sequence_inputs    = params.sequence_inputs ?: null
params.annotation        = params.annotation ?: null
params.protein_fasta     = params.protein_fasta ?: null
params.id_species        = params.id_species ?: null
params.selection_table   = params.selection_table ?: null
params.outgroup_genes    = params.outgroup_genes ?: null
params.outdir            = params.outdir ?: 'results'
params.famsa_args        = params.famsa_args ?: null
params.iqtree_args       = params.iqtree_args ?: null
params.species_tree        = params.species_tree ?: null
params.species_tree_outgroup = params.species_tree_outgroup ?: null
params.busco_lineage       = params.busco_lineage ?: null
params.busco_auto_lineage  = params.busco_auto_lineage ?: true
params.busco_threads       = params.busco_threads ?: 8
params.busco_mode          = params.busco_mode ?: 'auto'
params.busco_warn_complete = params.busco_warn_complete ?: 90
params.species_tree_bootstrap = params.species_tree_bootstrap ?: 1000
params.species_tree_threads   = params.species_tree_threads ?: params.busco_threads
params.species_tree_min_taxa  = params.species_tree_min_taxa ?: 3
params.species_tree_container = params.species_tree_container ?: null
params.generax_cpus          = params.generax_cpus ?: 64
params.generax_rec_model     = params.generax_rec_model ?: 'UndatedDL'
params.generax_strategy      = params.generax_strategy ?: 'SPR'
params.generax_seed          = params.generax_seed ?: 12345
params.generax_args          = params.generax_args ?: ''
params.generax_container     = params.generax_container ?: 'reconaming-generax:2.0.4'
params.gene_id_attribute       = params.gene_id_attribute ?: 'ID'
params.transcript_id_attribute = params.transcript_id_attribute ?: 'ID'
params.protein_id_attribute    = params.protein_id_attribute ?: ''
params.revision_registry       = params.revision_registry ?: null
params.revision_tag            = params.revision_tag ?: null
params.reconaming_opt          = params.reconaming_opt ?: "${projectDir}/conf/reconaming_opt.ctl"
params.reconaming_container    = params.reconaming_container ?: 'reconaming-core:6.2'
params.gff_output              = params.gff_output ?: 'update'

include { PREPARE_FOCAL_GENOME_MAP }      from './modules/local/prepare_focal_genome_map'
include { BUILD_INPUT_MANIFEST }          from './modules/local/build_input_manifest'
include { PREPARE_PROTEINS }               from './modules/local/prepare_proteins'
include { MERGE_ISOFORM_OUTPUTS }          from './modules/local/merge_isoform_outputs'
include { VALIDATE_PROTEIN_FASTA }         from './modules/local/validate_protein_fasta'
include { SELECT_REPRESENTATIVE_ISOFORMS } from './modules/local/select_representative_isoforms'
include { REQUIRE_COMPLETE_MANUAL_SELECTION } from './modules/local/require_complete_manual_selection'
include { INFER_SPECIES_TREE }             from './modules/local/infer_species_tree'
include { VALIDATE_SUPPLIED_SPECIES_TREE } from './modules/local/validate_supplied_species_tree'
include { DERIVE_GENE_TO_SPECIES }         from './modules/local/derive_gene_to_species'
include { ALIGN_PROTEINS }                 from './modules/local/align_proteins'
include { TRIM_ALIGNMENT }                 from './modules/local/trim_alignment'
include { INFER_GENE_TREE }                from './modules/local/infer_gene_tree'
include { RECONCILE_GENE_TREE }            from './modules/local/reconcile_gene_tree'
include { RUN_RECONAMING }                 from './modules/local/run_reconaming'
include { POSTPROCESS_RECONAMING_ANNOTATION } from './modules/local/postprocess_reconaming'
include { POSTPROCESS_RECONAMING_PROTEIN }    from './modules/local/postprocess_reconaming'

workflow {
    // BEGIN INPUT SELECTION
    // Resolve analysis mode and species-tree source once; reuse below.
    inputProvided = { Object value -> value != null && value.toString().trim() != '' }
    boolean suppliedTree = inputProvided.call(params.get('species_tree'))
    boolean annotationRequested = inputProvided.call(params.get('annotation'))
    sequence_rows = []
    if (inputProvided.call(params.get('sequence_inputs')) && (annotationRequested || !suppliedTree)) {
        sequence_table = file(params.sequence_inputs, checkIfExists: true)
        String tableText = sequence_table.text
        List<String> lines = tableText.readLines().findAll { it.trim() }
        if (!lines) error '--sequence_inputs is empty.'
        List<String> header = lines[0].replaceFirst('^\uFEFF', '').split('\t', -1).collect { it.trim() }
        List<String> required = ['species_id', 'input_file', 'input_type', 'role']
        if (header.size() != header.toSet().size() || !header.containsAll(required)) {
            error '--sequence_inputs requires unique columns: species_id, input_file, input_type, role.'
        }
        Set<String> seen = new HashSet<String>()
        lines.drop(1).eachWithIndex { String line, int index ->
            List<String> values = line.split('\t', -1).collect { it.trim() }
            if (values.size() != header.size()) {
                error "--sequence_inputs: invalid column count at row ${index + 2}."
            }
            Map row = [:]
            header.eachWithIndex { String key, int column -> row[key] = values[column] }
            row.input_type = row.input_type.toLowerCase()
            // A supplied tree does not need auxiliary inference inputs.
            if (!suppliedTree || row.role == 'focal') {
                if (!row.species_id || !row.input_file) {
                    error '--sequence_inputs contains an empty species_id or input_file.'
                }
                if (!seen.add(row.species_id)) {
                    error "Duplicate species_id in --sequence_inputs: ${row.species_id}"
                }
                if (!(row.input_type in ['genome', 'proteome'])) {
                    error "Unsupported input_type for ${row.species_id}: ${row.input_type}"
                }
                if (!(row.role in ['focal', 'outgroup'])) {
                    error "Unsupported sequence-input role: ${row.role}"
                }
                sequence_rows.add(row)
            }
        }
    }
    focal_rows = sequence_rows.findAll { it.role == 'focal' }
    boolean annotationComplete = annotationRequested && !focal_rows.isEmpty() && focal_rows.every { it.input_type == 'genome' }
    boolean proteinComplete = inputProvided.call(params.get('protein_fasta')) && inputProvided.call(params.get('id_species'))
    if (annotationComplete == proteinComplete) {
        error(annotationComplete ?
            'Both analysis input sets are complete; supply either focal genomes + annotations or protein FASTA + --id_species.' :
            'No complete analysis input set: supply focal genomes (input_type=genome, role=focal) + --annotation, or --protein_fasta + --id_species.')
    }
    if (!suppliedTree) {
        outgroup_rows = sequence_rows.findAll { it.role == 'outgroup' }
        if (focal_rows.isEmpty() || outgroup_rows.size() != 1) {
            error 'Without --species_tree, inference requires focal genome/proteome inputs and exactly one auxiliary outgroup in --sequence_inputs.'
        }
    }
    input_mode = annotationComplete ? 'annotation' : 'protein'
    tree_source = suppliedTree ? 'supplied' : 'inferred'
    // END INPUT SELECTION

    sequence_rows.findAll { tree_source == 'inferred' || (input_mode == 'annotation' && it.role == 'focal') }
        .each { row -> file(row.input_file, checkIfExists: true) }
    if (input_mode == 'protein') {
        file(params.protein_fasta, checkIfExists: true)
        file(params.id_species, checkIfExists: true)
    }
    if (tree_source == 'supplied') file(params.species_tree, checkIfExists: true)
    log.info "Resolved routes: ${input_mode} input; ${tree_source} species tree"
    if (input_mode == 'protein') {
        VALIDATE_PROTEIN_FASTA(Channel.of(tuple(file(params.protein_fasta), file(params.id_species))))
        representative_proteins = VALIDATE_PROTEIN_FASTA.out.proteins
        representative_manifest = VALIDATE_PROTEIN_FASTA.out.manifest
    }
    else if (input_mode == 'annotation') {
        PREPARE_FOCAL_GENOME_MAP(Channel.fromPath(params.sequence_inputs, checkIfExists: true))
        BUILD_INPUT_MANIFEST(PREPARE_FOCAL_GENOME_MAP.out.focal_genome_map)
        gffs = Channel.fromPath(params.annotation, checkIfExists: true).collect()
        source_gffs = gffs
        focal_genomes = PREPARE_FOCAL_GENOME_MAP.out.focal_genome_map
            .splitCsv(header: true, sep: '\t')
            .map { row -> tuple(row.species_id, file(row.genome_file)) }
        prepare_proteins_script = file("${projectDir}/bin/prepare_proteins.py", checkIfExists: true)
        PREPARE_PROTEINS(
            focal_genomes.combine(gffs).map { s, genome, files -> tuple(s, genome, files) },
            prepare_proteins_script,
            params.gene_id_attribute,
            params.transcript_id_attribute,
            params.protein_id_attribute
        )
        MERGE_ISOFORM_OUTPUTS(PREPARE_PROTEINS.out.proteins.collect(), PREPARE_PROTEINS.out.manifest.collect())
        selection_table = params.selection_table ? file(params.selection_table, checkIfExists: true) : []
        selection_table_supplied = params.selection_table ? true : false
        SELECT_REPRESENTATIVE_ISOFORMS(
            MERGE_ISOFORM_OUTPUTS.out.proteins,
            MERGE_ISOFORM_OUTPUTS.out.manifest,
            selection_table,
            file("${projectDir}/bin/select_representative_isoforms.py", checkIfExists: true),
            selection_table_supplied
        )
        REQUIRE_COMPLETE_MANUAL_SELECTION(
            SELECT_REPRESENTATIVE_ISOFORMS.out.status, "${params.outdir}/manual_selection")
        representative_proteins = SELECT_REPRESENTATIVE_ISOFORMS.out.proteins
        representative_manifest = SELECT_REPRESENTATIVE_ISOFORMS.out.manifest
    }

    if (tree_source == 'supplied') {
        VALIDATE_SUPPLIED_SPECIES_TREE(Channel.of(file(params.species_tree)), representative_manifest, params.species_tree_outgroup ?: '')
        species_tree = VALIDATE_SUPPLIED_SPECIES_TREE.out.species_tree
    }
    else {
        INFER_SPECIES_TREE(representative_manifest, Channel.fromPath(params.sequence_inputs, checkIfExists: true))
        species_tree = INFER_SPECIES_TREE.out.species_tree
    }

    DERIVE_GENE_TO_SPECIES(representative_manifest)
    gene_to_species = DERIVE_GENE_TO_SPECIES.out.gene_to_species
    ALIGN_PROTEINS(representative_proteins, params.famsa_args ?: '')
    TRIM_ALIGNMENT(ALIGN_PROTEINS.out.alignment)
    INFER_GENE_TREE(TRIM_ALIGNMENT.out.trimmed_alignment, params.iqtree_args ?: '')
    gene_tree = INFER_GENE_TREE.out.gene_tree

    outgroup_genes = params.outgroup_genes ? file(params.outgroup_genes, checkIfExists: true) : []
    reconciliation_input = gene_tree
        .combine(INFER_GENE_TREE.out.report)
        .combine(TRIM_ALIGNMENT.out.trimmed_alignment)
        .combine(species_tree)
        .combine(gene_to_species)
        .combine(representative_manifest)
        .map { tree, report, alignment, tree_species, mapping, manifest -> tuple(tree, report, alignment, tree_species, mapping, manifest, outgroup_genes) }

    RECONCILE_GENE_TREE(
        reconciliation_input,
        file("${projectDir}/bin/prepare_generax_inputs.py", checkIfExists: true),
        file("${projectDir}/bin/normalise_generax_output.py", checkIfExists: true),
        params.generax_cpus, params.generax_rec_model, params.generax_strategy,
        params.generax_seed, params.generax_args, params.generax_container
    )

    if (!params.revision_tag) error 'Job 15 requires --revision_tag (for example rev_20260922).'
    revision_registry = params.revision_registry ? file(params.revision_registry, checkIfExists: true) : []
    reconaming_input = RECONCILE_GENE_TREE.out.reconciled_tree
        .map { reconciled_tree -> tuple(reconciled_tree, outgroup_genes, revision_registry) }

    RUN_RECONAMING(
        reconaming_input,
        file("${projectDir}/bin/reconaming_core.py", checkIfExists: true),
        file("${projectDir}/bin/run_reconaming.py", checkIfExists: true),
        file(params.reconaming_opt, checkIfExists: true),
        params.revision_tag,
        params.reconaming_container
    )
    named_gene_tree = RUN_RECONAMING.out.named_tree
    updated_revision_registry = RUN_RECONAMING.out.revision_registry

    postprocess_common = RUN_RECONAMING.out.named_tree
        .combine(RUN_RECONAMING.out.name_table)
        .combine(RUN_RECONAMING.out.revision_registry)
        .combine(RUN_RECONAMING.out.temporary_partition)
        .combine(RUN_RECONAMING.out.report)
        .map { tree, table, registry, partition, report -> tuple(tree, table, registry, partition, report) }

    if (input_mode == 'annotation') {
        annotation_postprocess_input = postprocess_common.combine(source_gffs)
            .map { tree, table, registry, partition, report, source_gff_list -> tuple(tree, table, registry, partition, report, source_gff_list) }
        POSTPROCESS_RECONAMING_ANNOTATION(
            annotation_postprocess_input,
            file("${projectDir}/bin/postprocess_reconaming.py", checkIfExists: true),
            params.gff_output, params.gene_id_attribute, params.revision_tag,
            params.outdir, params.reconaming_container
        )
    }
    else {
        POSTPROCESS_RECONAMING_PROTEIN(
            postprocess_common,
            file("${projectDir}/bin/postprocess_reconaming.py", checkIfExists: true),
            'none', params.revision_tag, params.outdir, params.reconaming_container
        )
    }
}
