# Reconaming pipeline configuration
# Empty values are ignored. Command-line --params take precedence.
# This is a shared Maker-style control file: input-preparation entries are
# consumed upstream; naming entries are consumed by Job 15.

# -----------------------------------------------------------------------------
# Input preparation
# -----------------------------------------------------------------------------

# Input mode: auto, annotation, protein
input_mode=auto

# Annotation mode: genome FASTA plus IR annotation GFF/GFF3
genome_fasta=
ir_gff=

# Protein mode: pre-extracted isoform protein FASTA plus TSV mapping
protein_fasta=
sequence_mapping=

# Optional TSV selecting one sequence/transcript per gene
isoform_selection=

# Protein validation
minimum_protein_length=1
invalid_protein_policy=fail
allow_terminal_stop=true

# Automatic representative isoform choice
# On equal distance from the global median: longest protein, then lexical transcript ID.
automatic_isoform_tiebreak=longest_then_lexical

# Output
outdir=results

# Annotation output after Reconaming: update, minimal, or none.
# Protein input mode requires none.
gff_output=update

# -----------------------------------------------------------------------------
# Reconaming (Job 15)
# -----------------------------------------------------------------------------
# All values below are passed to the Reconaming core only when non-empty.

# Gene-family label inserted between species alias and ParentGroup number.
prefix=

# Minimum internal-node support required to avoid a naming boundary.
threshold=

# Path to a two-column, tab-separated species-to-alias file, or inline mapping.
aliases=

# true: ignore low-support single-copy ortholog nodes; false: treat them as boundaries.
ignore_sco=

# Lowest numeric ParentGroup identifier for a new empty revision registry.
minimal_id=

# true: name the rooted outgroup too; false/default: exclude it from naming.
include_outgroup=

# true/default: ladderize/core-order the tree; false: preserve existing ordering.
sort_tree=

# Optional whitespace-separated reconciled-tree leaf IDs used to root the tree.
# When blank, non-comment entries in --outgroup_genes are used automatically.
root_by=

# Node IDs excluded as ParentGroup roots. Comma/whitespace list or file path.
pgignore_nodes=

# Largest permitted ParentGroup subtree; 0/default means no limit.
max_pg_size=

# true/default: restore unambiguous historical ParentGroup assignments.
allow_revive_pg=

# Comma-separated Node IDs forced to be ParentGroup roots.
force_pg_root=

# true: display '-1' for singleton paralog groups; false/default: omit it.
show_single_para=

# Comma-separated Node IDs forced to be boundaries or non-boundaries.
force_boundary_true=
force_boundary_false=
