# Reconaming

Reconaming (**Recon**ciliated **naming**) is an ortholog-aware multigene-family reconciliation and naming pipeline. Starting from either gene annotations and genomic sequences or a prepared protein set, it infers or validates a species tree, infers and reconciles a gene tree, and assigns consistent ortholog-aware gene names. A revision registry carries forward unambiguous historical orthology group assignments between annotation or dataset revisions.

![Reconaming workflow](assets/Reconaming_FlowChart.png)

## Contents
...

## Installation and prerequisites

### Software requirements

- [Nextflow](https://www.nextflow.io/)
- Java, as required by the installed Nextflow release
- Either Docker or Apptainer/Singularity to run the configured containers

The workflow currently uses Nextflow's `local` executor by default. Docker and Singularity/Apptainer profiles are provided. Scheduler-based HPC execution is described as experimental below.

### Obtain Reconaming

```bash
git clone https://github.com/scie-mon/Reconaming.git
cd Reconaming
```

### Container runtime

Use one of the supplied profiles:

```bash
# Docker
nextflow run main.nf -profile docker [options]

# Apptainer/Singularity
nextflow run main.nf -profile singularity [options]
```

The four container build recipes are supplied in the repository:

- Species-tree inference: `containers/Dockerfile.species-tree`
- Gene-tree inference: `containers/Dockerfile.gene-tree`
- GeneRax reconciliation: `containers/Dockerfile.generax`
- Reconaming core: `containers/Dockerfile.reconaming-core`

Prebuilt Reconaming images are not currently published. The profiles do not build images automatically: before running the pipeline, build or otherwise provide the required images under the names configured for the workflow. The Dockerfiles are recipes, not images that Apptainer/Singularity can run directly; building them on a system without Docker requires a suitable container-build facility.

### Run the bundled test dataset

The repository includes small integration-test datasets under `test/`, including a species-tree-inference case in `test/spec_tree/` and supplied-species-tree cases in `test/provide_spec_tree/`. These datasets include test genome FASTA files, a GFF3, a species-input table, and—where relevant—outgroup genes or a supplied species tree.

The following is the intended smoke-test pattern. Run it in a new output directory and set a deliberately modest GeneRax CPU allocation appropriate to the machine:

```bash
nextflow run main.nf \
  -profile docker \
  --sequence_inputs test/spec_tree/test_species_tree_inputs.tsv \
  --annotation test/spec_tree/merged.rev0.gff3 \
  --outgroup_genes test/spec_tree/outgroup_genes.txt \
  --revision_tag test_rev_001 \
  --generax_cpus 4 \
  --outdir results_test
```

## Quick start

For a standard analysis that starts from genome FASTA files and annotations, run:

```bash
nextflow run main.nf \
  -profile docker \
  --sequence_inputs path/to/sequence_inputs.tsv \
  --annotation 'path/to/annotations/*.gff3' \
  --revision_tag rev_001 \
  --generax_cpus 8 \
  --outdir results \
  -resume
```

Replace the placeholder paths. `--revision_tag` is required. Use `-resume` when repeating or restarting a run with unchanged inputs and configuration; it lets Nextflow reuse valid completed tasks retained in the work directory.

A protein-route run, using pre-extracted representative proteins, has the following form:

```bash
nextflow run main.nf \
  -profile docker \
  --protein_fasta path/to/representative_proteins.faa \
  --id_species path/to/gene_to_species.tsv \
  --species_tree path/to/species_tree.nwk \
  --revision_tag rev_001 \
  --generax_cpus 8 \
  --gff_output none \
  --outdir results \
  -resume
```

## Input routes

Reconaming has two input routes. Both converge on protein alignment, gene-tree inference, reconciliation, and naming.

### Annotation route

Annotation mode is selected automatically when `--annotation` is supplied and `--sequence_inputs` contains focal genome inputs. All focal inputs used for annotation preparation must have `input_type=genome`. There is no explicit mode selector.

Required inputs:

| Parameter | Description |
|---|---|
| `--sequence_inputs` | Tab-separated sequence-input table identifying focal species and their genome FASTA files |
| `--annotation` | annotation GFF/GFF3 file or path/glob resolving to the annotation file(s) |
| `--revision_tag` | Identifier for this naming revision round |
| `--generax_cpus` | CPU allocation to use for GeneRax reconciliation |

The pipeline prepares proteins from the annotations, builds an input manifest, and selects one representative isoform per gene. Supply `--selection_table` when manual representative-isoform choices are required. Without it, the pipeline applies its automatic representative-selection procedure; if a required manual selection remains unresolved, the run stops and writes material for completion under the output directory.

### Protein route

Protein mode is selected automatically when both `--protein_fasta` and `--id_species` are supplied.

Required inputs:

| Parameter | Description |
|---|---|
| `--protein_fasta` | Protein FASTA containing the representative sequences to analyse |
| `--id_species` | Gene-ID-to-species-ID mapping file corresponding to the protein FASTA |
| `--revision_tag` | Identifier for this naming revision round |
| `--generax_cpus` | CPU allocation to use for GeneRax reconciliation |

`--id_species` is a TSV with `sequence_id` and `species_id` columns; sequence IDs must match protein FASTA identifiers. Protein mode validates the supplied proteins and mapping, but does not have source GFF/GFF3 records to rewrite. Use `--gff_output none` in this mode.

### Shared optional inputs

| Parameter | Purpose |
|---|---|
| `--species_tree` | Use a supplied species tree instead of inferring one |
| `--species_tree_outgroup` | Outgroup information used when validating or rooting a supplied species tree |
| `--outgroup_genes` | File of gene IDs used for gene-tree rooting; an empty project asset is used when omitted |
| `--revision_registry` | Revision registry from a previous naming round; an empty project registry is used for a first round |
| `--reconaming_opt` | Path to the advanced Reconaming-core control file; defaults to `conf/reconaming_opt.ctl` |
| `--outdir` | Destination for final results; defaults to `results` |
| `--sequence_inputs` (protein route) | In protein route this can be also used to configure the species tree inference input |

Analysis-input selection and species-tree selection are independent.

| Analysis route | Complete required input set |
|---|---|
| Annotation | `--annotation` + focal genomes in `--sequence_inputs` |
| Protein | `--protein_fasta` + `--id_species` |

Exactly one analysis input set must be complete. If neither or both are complete, startup validation stops the run. Incomplete alternative input sets do not select another route.

Species-tree selection uses the following precedence:

A species tree supplied through `--species_tree` takes precedence over internal inference. The pipeline validates the supplied tree before using it; validation failure stops the run rather than triggering inference.

When no species tree is supplied, the pipeline infers one from suitable focal-species and auxiliary-outgroup inputs listed in --sequence_inputs. If the required inference inputs are missing or incomplete, the run stops with an error identifying the missing prerequisites.

`--sequence_inputs` is a TSV with these columns:

| Column | Meaning |
|---|---|
| `species_id` | Species identifier |
| `input_file` | Genome or proteome FASTA path, accessible from the run directory |
| `input_type` | `genome` or `proteome` |
| `role` | `focal` or `outgroup` |

For inference, supply one input for every analysed focal species and exactly one auxiliary outgroup. Species identifiers must be unique, the focal set must match the representative-protein manifest exactly, and the auxiliary outgroup must not be a gene-family species. With `--busco_mode auto`, each row's `input_type` is used; genome, proteome and mixed inputs are supported. Here, `auto` follows the table rather than detecting sequence type. The inference workflow roots on the auxiliary outgroup and then prunes it.

The per-species proteomes used for species-tree inference are separate from the representative gene-family `--protein_fasta`. Protein mode can therefore infer a species tree when the required per-species sequence inputs are supplied.

Annotation preparation still requires focal genomes. If a tree is supplied, auxiliary inference rows are not required or used; a table containing them is not a route conflict. Protein mode with a supplied tree does not need `--sequence_inputs`, and an unused sequence table is not read unless annotation inputs also need evaluation.

## Parameters

Provide Nextflow parameters as `--name value`. The parameters below are grouped by purpose. Values in `nextflow.config` and `main.nf` are the authoritative defaults for the checked-out pipeline version.

### Core run parameters

| Parameter | Default | Meaning |
|---|---:|---|
| `--outdir` | `results` | Final-results directory |
| `--revision_tag` | none; required | Name for the current revision round, for example `rev_001` |
| `--revision_registry` | empty project registry | Previous revision registry to continue from |
| `--gff_output` | `update` | Annotation-mode output policy: `update`, `minimal`, or `none` |

### Annotation and protein inputs

| Parameter | Default | Meaning |
|---|---:|---|
| `--sequence_inputs` | none | Species-labelled genome/proteome table; focal genomes prepare annotations, while focal + outgroup sequences enable inference |
| `--annotation` | none | annotation GFF/GFF3 for annotation mode |
| `--selection_table` | empty project table | Optional manual representative-isoform selection table |
| `--protein_fasta` | none | Representative protein FASTA for protein mode |
| `--id_species` | none | Gene-to-species mapping for protein mode |
| `--gene_id_attribute` | `ID` | GFF attribute used as the gene identifier |
| `--transcript_id_attribute` | `ID` | GFF attribute used as the transcript identifier |
| `--protein_id_attribute` | empty | Optional GFF attribute used as the protein identifier |

### Species-tree parameters

| Parameter | Default | Meaning |
|---|---:|---|
| `--species_tree` | none | Supplied species tree; skips internal species-tree inference |
| `--species_tree_outgroup` | none | Outgroup for supplied-tree validation/rooting |
| `--busco_lineage` | none | BUSCO lineage dataset to use, if specified |
| `--busco_auto_lineage` | `true` | Let BUSCO choose a lineage automatically |
| `--busco_mode` | `auto` | Follow each row's `input_type`; `genome` or `proteome` forces one type for all inference inputs |
| `--busco_threads` | `4` | Threads for BUSCO |
| `--busco_warn_complete` | `90` | Completeness threshold that triggers a warning |
| `--species_tree_threads` | `4` | Threads for species-tree inference |
| `--species_tree_bootstrap` | `1000` | Number of species-tree bootstrap replicates |
| `--species_tree_min_taxa` | `3` | Minimum number of taxa required for species-tree inference |

### Gene-tree and reconciliation parameters

| Parameter | Default | Meaning |
|---|---:|---|
| `--famsa_args` | empty | Additional arguments passed to FAMSA |
| `--iqtree_args` | empty | Additional arguments passed to IQ-TREE |
| `--gene_tree_container` | `reconaming-gene-tree:0.1.0` | Container for the labelled gene-tree process |
| `--generax_cpus` | **set explicitly** | CPUs allocated to GeneRax reconciliation |
| `--generax_rec_model` | `UndatedDL` | GeneRax reconciliation model |
| `--generax_strategy` | `SPR` | GeneRax search strategy |
| `--generax_seed` | `12345` | Random seed for GeneRax |
| `--generax_args` | empty | Additional arguments passed to GeneRax |
| `--generax_container` | `reconaming-generax:2.0.4` | Container for GeneRax reconciliation |

`--generax_cpus` should always be selected for the resources available to the run. The README intentionally treats it as an explicit user decision rather than recommending a high fixed default.

### Advanced naming configuration

| Parameter | Default | Meaning |
|---|---:|---|
| `--reconaming_opt` | `conf/reconaming_opt.ctl` | Advanced Reconaming-core control file |
| `--reconaming_container` | `reconaming-core:6.2` | Container for the Reconaming core |
| `--outgroup_genes` | empty project file | Gene IDs used to root the reconciled gene tree |

## Revision rounds and revision registry

### Why revisioning exists

Gene-family naming might be required to be repeated as annotations improve, taxa are added, gene models change, or phylogenetic evidence is updated. Reconaming treats each such analysis as a revision round. The goal is to preserve unambiguous historical ParentGroup assignments where possible while allowing names to change when the underlying reconciliation genuinely changes. Deprecated gene numbers are never re-used and stay associated with their gene as an alias. A revision matrix keeps track of this, declaring legacy names of each gene and orthology of all contained genes with the current name set.

### Starting a revision round

Every run must have a `--revision_tag`. Choose a stable, human-readable tag that identifies the round, for example `rev_001`, `rev_002`, or a date-based tag such as `rev_20260924`.

For a first round, omit `--revision_registry`; the workflow starts from the project's empty registry. For a later round, supply the registry produced by the previous accepted round:

```bash
nextflow run main.nf -profile docker \
  --revision_tag rev_002 \
  --revision_registry path/to/rev_001_revision_registry.csv \
  [other input and run parameters]
```

Do not reuse a revision tag for a biologically distinct round. Retain the inputs, configuration, registry, named tree, name table, and report associated with every accepted revision.

### The revision registry

The revision registry records historical ParentGroup assignments. It is an input to a new round and an output from the Reconaming stage. Carrying it forward gives the naming core the information needed to restore unambiguous historical assignments rather than treating each revision as an unrelated naming exercise.

The actual policy is controlled by the Reconaming control file. In particular:

- `minimal_id` defines the lowest numeric ParentGroup identifier for a new empty registry.
- `allow_revive_pg` permits restoration of unambiguous historical ParentGroup assignments.
- `prefix` supplies the optional gene-family label embedded in generated names.
- `aliases` provides species-to-alias mappings.

### Advanced control file

`conf/reconaming_opt.ctl` is an advanced configuration file for the Reconaming core. Empty values are ignored. It contains options for support thresholds, aliases, outgroup treatment, ParentGroup limits and overrides, rooting, tree ordering, and boundary behaviour.

Do not use its input-preparation-looking entries as a substitute for the Nextflow command-line interface. Invoke the pipeline with `--sequence_inputs`, `--annotation`, `--protein_fasta`, `--id_species`, and the other documented Nextflow options; use `--reconaming_opt` only to select the naming-core configuration file.

## Output

Reconaming writes final, user-facing outputs under `--outdir` (default: `results`). Exact filenames and subdirectory names are part of the pipeline version and should be recorded with the run.

### Primary outputs

The Reconaming stage produces the material needed to interpret and continue a naming round:

- A named reconciled gene tree
- A gene-name table
- An updated revision registry for a subsequent revision round
- A Reconaming report
- A temporary partition artifact used by post-processing

### Annotation-mode outputs

When annotation mode is selected and `--gff_output update` is used, post-processing applies Reconaming results to the supplied annotation data. `--gff_output minimal` requests a reduced annotation-oriented output, while `--gff_output none` suppresses GFF/GFF3 output.

### Protein-mode outputs

Protein mode produces the naming outputs but has no source GFF/GFF3 records to update. Set `--gff_output none`.

### Run state and intermediate files

The Nextflow `work/` directory and `.nextflow/` contain execution state, staged inputs, logs, task scripts, caches, and intermediates. Remove them only after deciding that the associated run no longer needs to be resumed.

If reports, traces, or timelines are enabled for a run, treat them as execution diagnostics. They are useful for troubleshooting and provenance, but they are distinct from the final named tree, name table, revision registry, and annotation outputs.

## HPC and container configuration

HPC schedulers and non-local executors are not yet fully supported or fully validated by Reconaming. The checked-in configuration defaults to Nextflow's local executor and supplies Docker and Singularity/Apptainer profiles as starting points.

Before production use on an HPC system:

- Run the bundled test dataset with the intended container runtime.
- Configure the executor and process-resource settings for the local scheduler policy.
- Verify that input, output, temporary, and container-cache filesystems are visible on compute nodes.
- Set `--generax_cpus` to the CPU allocation granted to the GeneRax task; do not assume a machine-wide or fixed pipeline default.
- Validate memory limits, walltime, container mounts, and writable output ownership.

A scheduler-specific profile—for example, Slurm, PBS, LSF, or SGE—should be considered site configuration and tested before it is presented as supported.

## License

Reconaming is distributed under the GNU General Public License v3.0. See [`LICENSE`](LICENSE).

## Citation

Citation metadata is provided in [`CITATION.cff`](CITATION.cff). Please cite the software using that file and cite the underlying methods and software used by your analysis where appropriate.
