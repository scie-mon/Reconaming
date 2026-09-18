#!/usr/bin/env python3
import argparse
import csv
import re
import shutil
from pathlib import Path


def leaves(tree_path):
    text = Path(tree_path).read_text().strip()
    labels = re.findall(r"(?:\(|,)\s*(?:'([^']+)'|([^\s():,;\[\]]+))", text)
    result = {a or b for a, b in labels}
    if not result:
        raise ValueError(f"No leaf labels found in {tree_path}")
    return result


def find_raw_nhx(results_dir, excluded):
    candidates = [p for p in Path(results_dir).rglob('*.nhx') if p.resolve() not in excluded]
    preferred = [p for p in candidates if re.search(r'reconcil', p.name, re.IGNORECASE)]
    candidates = preferred or candidates
    if len(candidates) != 1:
        names = ', '.join(str(p) for p in candidates) or 'none'
        raise ValueError(f"Expected exactly one GeneRax reconciled NHX result; found: {names}")
    return candidates[0]


def read_summary(path):
    rows = []
    with open(path, newline='') as handle:
        reader = csv.reader(handle, delimiter='\t')
        header = next(reader, None)
        if header != ['metric', 'value']:
            raise ValueError(f"Invalid input summary header in {path}")
        rows.extend(reader)
    return rows


parser = argparse.ArgumentParser()
parser.add_argument('--generax-results', required=True)
parser.add_argument('--expected-gene-tree', required=True)
parser.add_argument('--input-summary', required=True)
parser.add_argument('--raw-out', required=True)
parser.add_argument('--normalised-out', required=True)
parser.add_argument('--summary-out', required=True)
args = parser.parse_args()

raw_out = Path(args.raw_out).resolve()
normalised_out = Path(args.normalised_out).resolve()
source = find_raw_nhx(args.generax_results, {raw_out, normalised_out})
expected = leaves(args.expected_gene_tree)
raw_content = source.read_text()
raw_leaves = leaves(source)
if raw_leaves != expected:
    raise ValueError(
        'GeneRax/output leaf mismatch: '
        f"missing={sorted(expected - raw_leaves)[:10]}, "
        f"unexpected={sorted(raw_leaves - expected)[:10]}"
    )

shutil.copyfile(source, raw_out)
normalised_content, substitutions = re.subn(r'\)n[0-9]+', ')100', raw_content)
if not normalised_content.endswith('\n'):
    normalised_content += '\n'
normalised_out.write_text(normalised_content)
normalised_leaves = leaves(normalised_out)
if normalised_leaves != expected:
    raise ValueError('NHX normalisation changed the gene-tree leaf set')

rows = read_summary(args.input_summary)
rows.extend([
    ('raw_generax_nhx', str(source)),
    ('raw_generax_nhx_leaves', len(raw_leaves)),
    ('nhx_internal_labels_replaced', substitutions),
    ('normalised_nhx_leaves', len(normalised_leaves)),
    ('validation', 'passed'),
])
with open(args.summary_out, 'w', newline='') as handle:
    writer = csv.writer(handle, delimiter='\t')
    writer.writerow(['metric', 'value'])
    writer.writerows(rows)
