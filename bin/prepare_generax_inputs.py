#!/usr/bin/env python3
import argparse
import csv
import re
from collections import defaultdict
from pathlib import Path


def leaves(newick_path):
    text = Path(newick_path).read_text().strip()
    labels = re.findall(r"(?:\(|,)\s*(?:'([^']+)'|([^\s():,;\[\]]+))", text)
    result = [a or b for a, b in labels]
    if not result or len(result) != len(set(result)):
        raise ValueError(f"Could not obtain a unique leaf set from {newick_path}")
    return set(result)


def fasta_ids(fasta_path):
    result = set()
    for line in Path(fasta_path).read_text().splitlines():
        if line.startswith('>'):
            identifier = line[1:].split()[0]
            if not identifier or identifier in result:
                raise ValueError(f"Invalid or duplicate FASTA identifier: {identifier!r}")
            result.add(identifier)
    if not result:
        raise ValueError(f"No FASTA records found in {fasta_path}")
    return result


def two_column_mapping(mapping_path):
    with open(mapping_path, newline='') as handle:
        rows = list(csv.DictReader(handle, delimiter='\t'))
    if not rows or set(rows[0]) != {'sequence_id', 'species_id'}:
        raise ValueError(f"{mapping_path} must have exactly sequence_id and species_id columns")
    result = {row['sequence_id'].strip(): row['species_id'].strip() for row in rows}
    if not all(result) or len(result) != len(rows):
        raise ValueError(f"Invalid or duplicate sequence_id in {mapping_path}")
    return result


def manifest_mapping(manifest_path):
    with open(manifest_path, newline='') as handle:
        rows = list(csv.DictReader(handle, delimiter='\t'))
    required = {'sequence_id', 'gene_id', 'species_id'}
    if not rows or required - set(rows[0]):
        raise ValueError(f"{manifest_path} must contain {', '.join(sorted(required))}")
    result = {}
    for row in rows:
        sequence_id = row['sequence_id'].strip()
        gene_id = row['gene_id'].strip()
        species_id = row['species_id'].strip()
        if not sequence_id or not gene_id or not species_id:
            raise ValueError('Representative manifest contains an empty sequence_id, gene_id, or species_id')
        if sequence_id in result:
            raise ValueError(f'Duplicate sequence_id in representative manifest: {sequence_id}')
        result[sequence_id] = (gene_id, species_id)
    gene_ids = [gene_id for gene_id, _ in result.values()]
    duplicates = sorted({x for x in gene_ids if gene_ids.count(x) > 1})
    if duplicates:
        raise ValueError(f"User-facing gene_id is not globally unique and cannot label a GeneRax tree: {duplicates[:10]}")
    invalid = [x for x in gene_ids if re.search(r"[\s():,;\[\]']", x)]
    if invalid:
        raise ValueError(f"GeneRax leaf labels contain unsupported Newick characters: {invalid[:10]}")
    return result


def outgroup_ids(outgroup_path):
    lines = [x.strip() for x in Path(outgroup_path).read_text().splitlines() if x.strip() and not x.lstrip().startswith('#')]
    if lines and lines[0].lower().split('\t')[0] in {'gene', 'gene_id', 'sequence', 'sequence_id', 'id'}:
        lines = lines[1:]
    return {x.split('\t')[0].split()[0] for x in lines}


def selected_model(report_path):
    text = Path(report_path).read_text()
    patterns = [r'Best-fit model according to (?:BIC|AIC|AICc)\s*:\s*([^\s]+)', r'Best-fit model\s*:\s*([^\s]+)', r'Model of substitution\s*:\s*([^\s]+)']
    found = list(dict.fromkeys(x for p in patterns for x in re.findall(p, text, flags=re.IGNORECASE)))
    if len(found) != 1:
        raise ValueError(f"Expected one selected IQ-TREE substitution model in {report_path}; found {found}")
    return found[0]


def translate_fasta(source, target, labels):
    with open(source) as input_handle, open(target, 'w') as output_handle:
        for line in input_handle:
            if line.startswith('>'):
                old = line[1:].split()[0]
                suffix = line[1 + len(old):]
                output_handle.write(f'>{labels[old]}{suffix}')
            else:
                output_handle.write(line)


def translate_newick(source, target, labels):
    text = Path(source).read_text()
    pattern = re.compile(r"(?P<prefix>[\(,]\s*)(?:'(?P<quoted>[^']+)'|(?P<bare>[^\s():,;\[\]]+))")
    def replace(match):
        old = match.group('quoted') or match.group('bare')
        if old not in labels:
            raise ValueError(f'Gene-tree leaf absent from canonical manifest: {old}')
        return f"{match.group('prefix')}{labels[old]}"
    target.write_text(pattern.sub(replace, text))


parser = argparse.ArgumentParser()
parser.add_argument('--gene-tree', required=True)
parser.add_argument('--iqtree-report', required=True)
parser.add_argument('--alignment', required=True)
parser.add_argument('--species-tree', required=True)
parser.add_argument('--gene-to-species', required=True)
parser.add_argument('--representative-manifest', required=True)
parser.add_argument('--outgroup-genes', required=True)
parser.add_argument('--families-out', required=True)
parser.add_argument('--mapping-out', required=True)
parser.add_argument('--generax-tree-out', required=True)
parser.add_argument('--generax-alignment-out', required=True)
parser.add_argument('--summary-out', required=True)
args = parser.parse_args()

gene_leaves = leaves(args.gene_tree)
alignment_leaves = fasta_ids(args.alignment)
species_leaves = leaves(args.species_tree)
canonical_mapping = two_column_mapping(args.gene_to_species)
manifest = manifest_mapping(args.representative_manifest)
model = selected_model(args.iqtree_report)
if gene_leaves != alignment_leaves or gene_leaves != set(canonical_mapping) or gene_leaves != set(manifest):
    raise ValueError('Canonical IDs differ among alignment, gene tree, gene-to-species mapping, and representative manifest')
for sequence_id, (_, species_id) in manifest.items():
    if canonical_mapping[sequence_id] != species_id:
        raise ValueError(f'Species mismatch for canonical sequence_id: {sequence_id}')
if set(canonical_mapping.values()) - species_leaves:
    raise ValueError(f"Mapped species absent from species tree: {sorted(set(canonical_mapping.values()) - species_leaves)}")

labels = {sequence_id: gene_id for sequence_id, (gene_id, _) in manifest.items()}
generax_gene_ids = set(labels.values())
outgroups = outgroup_ids(args.outgroup_genes)
if outgroups - generax_gene_ids:
    raise ValueError(f"Outgroup gene IDs absent from user-facing GeneRax tree: {sorted(outgroups - generax_gene_ids)[:10]}")

translate_fasta(args.alignment, args.generax_alignment_out, labels)
translate_newick(args.gene_tree, Path(args.generax_tree_out), labels)
if leaves(args.generax_tree_out) != fasta_ids(args.generax_alignment_out):
    raise ValueError('User-facing GeneRax tree/alignment IDs differ after translation')

by_species = defaultdict(list)
for sequence_id, gene_id in labels.items():
    by_species[canonical_mapping[sequence_id]].append(gene_id)
with open(args.mapping_out, 'w') as handle:
    for species in sorted(by_species):
        handle.write(f"{species}:" + ';'.join(sorted(by_species[species])) + '\n')
with open(args.families_out, 'w') as handle:
    handle.write('[FAMILIES]\n- gene_family\n')
    handle.write(f"starting_gene_tree = {Path(args.generax_tree_out).name}\n")
    handle.write(f"alignment = {Path(args.generax_alignment_out).name}\n")
    handle.write(f"mapping = {Path(args.mapping_out).name}\n")
    handle.write(f"subst_model = {model}\n")
with open(args.summary_out, 'w', newline='') as handle:
    writer = csv.writer(handle, delimiter='\t')
    writer.writerow(['metric', 'value'])
    writer.writerows([('canonical_sequence_ids', len(gene_leaves)), ('user_facing_gene_ids', len(generax_gene_ids)), ('mapped_species', len(by_species)), ('outgroup_genes', len(outgroups)), ('subst_model', model), ('validation', 'passed')])
