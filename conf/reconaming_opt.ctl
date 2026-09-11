# Reconaming input-preparation configuration
# Empty values are ignored. Command-line --params take precedence.

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
