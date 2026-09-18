#!/usr/bin/env python3
import argparse
import csv
import re
from collections import defaultdict
from pathlib import Path


def leaves(newick_path):
    text = Path(newick_path).read_text().strip()
    labels = re.findall(r"(?:\(|,)\s*(?:'([^']+)'|([^\s():,;]+))", text)
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


def mapping_rows(mapping_path):
    lines = [x for x in Path(mapping_path).read_text().splitlines() if x.strip() and not x.lstrip().startswith('#')]
    if not lines:
        raise ValueError(f"No mappings found in {mapping_path}")
    dialect = csv.excel_tab if any('\t' in x for x in lines[:2]) else csv.excel
    rows = list(csv.reader(lines, dialect=dialect))
    header = [x.strip().lower() for x in rows[0]]
    gene_names = {'gene', 'gene_id', 'sequence', 'sequence_id', 'id'}
    species_names = {'species', 'species_id', 'taxon', 'taxon_id'}
    if len(header) >= 2 and set(header) & gene_names and set(header) & species_names:
        gene_i = next(i for i, x in enumerate(header) if x in gene_names)
        species_i = next(i for i, x in enumerate(header) if x in species_names)
        rows = rows[1:]
    else:
        gene_i, species_i = 0, 1
    result = {}
    for row in rows:
        if len(row) <= max(gene_i, species_i):
            raise ValueError(f"Malformed mapping row: {row}")
        gene, species = row[gene_i].strip(), row[species_i].strip()
        if not gene or not species or ':' in species or ';' in species:
            raise ValueError(f"Invalid GeneRax mapping entry: {row}")
        if gene in result and result[gene] != species:
            raise ValueError(f"Conflicting species assignments for {gene}")
        result[gene] = species
    return result


def outgroup_ids(outgroup_path):
    lines = [x.strip() for x in Path(outgroup_path).read_text().splitlines() if x.strip() and not x.lstrip().startswith('#')]
    if lines and lines[0].lower().split('\t')[0] in {'gene', 'gene_id', 'sequence', 'sequence_id', 'id'}:
        lines = lines[1:]
    return {x.split('\t')[0].split()[0] for x in lines}


def selected_model(report_path):
    text = Path(report_path).read_text()
    patterns = [
        r'Best-fit model according to (?:BIC|AIC|AICc)\s*:\s*([^\s]+)',
        r'Best-fit model\s*:\s*([^\s]+)',
        r'Model of substitution\s*:\s*([^\s]+)',
    ]
    found = []
    for pattern in patterns:
        found.extend(re.findall(pattern, text, flags=re.IGNORECASE))
    found = list(dict.fromkeys(found))
    if len(found) != 1:
        raise ValueError(f"Expected one selected IQ-TREE substitution model in {report_path}; found {found}")
    return found[0]


parser = argparse.ArgumentParser()
parser.add_argument('--gene-tree', required=True)
parser.add_argument('--iqtree-report', required=True)
parser.add_argument('--alignment', required=True)
parser.add_argument('--species-tree', required=True)
parser.add_argument('--gene-to-species', required=True)
parser.add_argument('--outgroup-genes', required=True)
parser.add_argument('--families-out', required=True)
parser.add_argument('--mapping-out', required=True)
parser.add_argument('--summary-out', required=True)
args = parser.parse_args()

gene_leaves = leaves(args.gene_tree)
alignment_leaves = fasta_ids(args.alignment)
species_leaves = leaves(args.species_tree)
gene_species = mapping_rows(args.gene_to_species)
outgroups = outgroup_ids(args.outgroup_genes)
model = selected_model(args.iqtree_report)

if gene_leaves != alignment_leaves:
    raise ValueError(f"Gene tree/alignment leaf mismatch: tree_only={sorted(gene_leaves - alignment_leaves)[:10]}, alignment_only={sorted(alignment_leaves - gene_leaves)[:10]}")
missing_mapping = gene_leaves - set(gene_species)
extra_mapping = set(gene_species) - gene_leaves
if missing_mapping or extra_mapping:
    raise ValueError(f"Gene/species mapping mismatch: unmapped={sorted(missing_mapping)[:10]}, unused={sorted(extra_mapping)[:10]}")
missing_species = set(gene_species.values()) - species_leaves
if missing_species:
    raise ValueError(f"Mapped species absent from species tree: {sorted(missing_species)}")
missing_outgroups = outgroups - gene_leaves
if missing_outgroups:
    raise ValueError(f"Outgroup IDs absent from gene tree: {sorted(missing_outgroups)[:10]}")

by_species = defaultdict(list)
for gene in sorted(gene_leaves):
    by_species[gene_species[gene]].append(gene)
with open(args.mapping_out, 'w') as handle:
    for species in sorted(by_species):
        handle.write(f"{species}:" + ';'.join(by_species[species]) + '\n')

with open(args.families_out, 'w') as handle:
    handle.write('[FAMILIES]\n- gene_family\n')
    handle.write(f"starting_gene_tree = {Path(args.gene_tree).name}\n")
    handle.write(f"alignment = {Path(args.alignment).name}\n")
    handle.write(f"mapping = {Path(args.mapping_out).name}\n")
    handle.write(f"subst_model = {model}\n")

with open(args.summary_out, 'w', newline='') as handle:
    writer = csv.writer(handle, delimiter='\t')
    writer.writerow(['metric', 'value'])
    writer.writerows([
        ('gene_tree_leaves', len(gene_leaves)),
        ('alignment_sequences', len(alignment_leaves)),
        ('mapped_species', len(by_species)),
        ('outgroup_genes', len(outgroups)),
        ('subst_model', model),
        ('validation', 'passed'),
    ])
