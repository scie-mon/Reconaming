#!/usr/bin/env nextflow
nextflow.enable.dsl=2

params.input_mode        = params.input_mode ?: 'annotation'
params.species_inputs    = params.species_inputs ?: null
params.ir_gff            = params.ir_gff ?: null
params.protein_fasta     = params.protein_fasta ?: null
params.id_species        = null
params.isoform_selection = params.isoform_selection ?: null
params.outgroup_genes    = null
params.famsa_args        = null
params.iqtree_args       = null
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

include { PREPARE_FOCAL_GENOME_MAP }      from './modules/local/prepare_focal_genome_map'
include { BUILD_INPUT_MANIFEST }          from './modules/local/build_input_manifest'
include { PREPARE_PROTEINS }               from './modules/local/prepare_proteins'
include { MERGE_ISOFORM_OUTPUTS }          from './modules/local/merge_isoform_outputs'
include { VALIDATE_PROTEIN_FASTA }         from './modules/local/validate_protein_fasta'
include { SELECT_REPRESENTATIVE_ISOFORMS } from './modules/local/select_representative_isoforms'
include { INFER_SPECIES_TREE }             from './modules/local/infer_species_tree'
include { VALIDATE_SUPPLIED_SPECIES_TREE } from './modules/local/validate_supplied_species_tree'
include { DERIVE_GENE_TO_SPECIES }         from './modules/local/derive_gene_to_species'
include { ALIGN_PROTEINS }                 from './modules/local/align_proteins'
include { TRIM_ALIGNMENT }                 from './modules/local/trim_alignment'
include { INFER_GENE_TREE }                from './modules/local/infer_gene_tree'
include { RECONCILE_GENE_TREE }            from './modules/local/reconcile_gene_tree'

workflow {
    if (params.input_mode == 'protein') {
        if (!params.protein_fasta || !params.id_species) error 'Protein mode requires --protein_fasta and --id_species (gene_id to species_id mapping).'
        VALIDATE_PROTEIN_FASTA(Channel.of(tuple(file(params.protein_fasta), file(params.id_species))))
        representative_proteins = VALIDATE_PROTEIN_FASTA.out.proteins
        representative_manifest = VALIDATE_PROTEIN_FASTA.out.manifest
    }
    else if (params.input_mode == 'annotation') {
        if (!params.species_inputs || !params.ir_gff) error 'Annotation mode requires --species_inputs and --ir_gff.'
        PREPARE_FOCAL_GENOME_MAP(Channel.fromPath(params.species_inputs, checkIfExists: true))
        BUILD_INPUT_MANIFEST(PREPARE_FOCAL_GENOME_MAP.out.focal_genome_map)
        gffs = Channel.fromPath(params.ir_gff, checkIfExists: true).collect()
        focal_genomes = PREPARE_FOCAL_GENOME_MAP.out.focal_genome_map
            .splitCsv(header: true, sep: '\t')
            .map { row -> tuple(row.species_id, file(row.genome_file)) }
        prepare_proteins_script = file("${projectDir}/bin/prepare_proteins.py", checkIfExists: true)
        PREPARE_PROTEINS(
            focal_genomes.combine(gffs).map { s, genome, files -> tuple(s, genome, files) },
            prepare_proteins_script
        )
        MERGE_ISOFORM_OUTPUTS(PREPARE_PROTEINS.out.proteins.collect(), PREPARE_PROTEINS.out.manifest.collect())
        selection = params.isoform_selection ? file(params.isoform_selection).toAbsolutePath().toString() : ''
        SELECT_REPRESENTATIVE_ISOFORMS(MERGE_ISOFORM_OUTPUTS.out.proteins, MERGE_ISOFORM_OUTPUTS.out.manifest, selection)
        representative_proteins = SELECT_REPRESENTATIVE_ISOFORMS.out.proteins
        representative_manifest = SELECT_REPRESENTATIVE_ISOFORMS.out.manifest
    }
    else error "Unknown --input_mode '${params.input_mode}'; use annotation or protein."

    if (params.species_tree) {
        VALIDATE_SUPPLIED_SPECIES_TREE(Channel.of(file(params.species_tree)), representative_manifest, params.species_tree_outgroup ?: '')
        species_tree = VALIDATE_SUPPLIED_SPECIES_TREE.out.species_tree
    }
    else {
        if (!params.species_inputs) error 'Internal species-tree inference requires --species_inputs; alternatively provide --species_tree.'
        INFER_SPECIES_TREE(representative_manifest, Channel.fromPath(params.species_inputs, checkIfExists: true))
        species_tree = INFER_SPECIES_TREE.out.species_tree
    }

    DERIVE_GENE_TO_SPECIES(representative_manifest)
    gene_to_species = DERIVE_GENE_TO_SPECIES.out.gene_to_species

    ALIGN_PROTEINS(representative_proteins, params.famsa_args ?: '')
    TRIM_ALIGNMENT(ALIGN_PROTEINS.out.alignment)
    INFER_GENE_TREE(TRIM_ALIGNMENT.out.trimmed_alignment, params.iqtree_args ?: '')
    gene_tree = INFER_GENE_TREE.out.gene_tree

    outgroup_genes = params.outgroup_genes ? file(params.outgroup_genes, checkIfExists: true) : file("${projectDir}/assets/empty_outgroup_genes.txt", checkIfExists: true)
    reconciliation_input = gene_tree
        .combine(INFER_GENE_TREE.out.report)
        .combine(TRIM_ALIGNMENT.out.trimmed_alignment)
        .combine(species_tree)
        .combine(gene_to_species)
        .map { tree, report, alignment, tree_species, mapping -> tuple(tree, report, alignment, tree_species, mapping, outgroup_genes) }

    RECONCILE_GENE_TREE(
        reconciliation_input,
        file("${projectDir}/bin/prepare_generax_inputs.py", checkIfExists: true),
        file("${projectDir}/bin/normalise_generax_output.py", checkIfExists: true),
        params.generax_cpus,
        params.generax_rec_model,
        params.generax_strategy,
        params.generax_seed,
        params.generax_args,
        params.generax_container
    )
    reconciled_gene_tree = RECONCILE_GENE_TREE.out.reconciled_tree

    // Jobs 15–18 consume reconciled_gene_tree and the existing upstream manifests.
}
