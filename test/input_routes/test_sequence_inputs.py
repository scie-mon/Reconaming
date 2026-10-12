#!/usr/bin/env python3
import csv
import re
import subprocess
import sys
import tempfile
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
HEADER = 'species_id\tinput_file\tinput_type\trole\n'
FOCAL = 'A\ta.fa\tgenome\tfocal\nB\tb.fa\tgenome\tfocal\n'
OUTGROUP = 'O\to.fa\tgenome\toutgroup\n'

class SequenceInputTests(unittest.TestCase):
    def run_script(self, script, table, species=('A', 'B'), mode='auto'):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            sequences = root / 'inputs.tsv'
            sequences.write_text(table)
            output = root / 'out.tsv'
            command = [sys.executable, str(ROOT / 'bin' / script), '--sequence-inputs', str(sequences), '--output', str(output)]
            if script == 'select_species_tree_inputs.py':
                manifest = root / 'manifest.tsv'
                manifest.write_text('species_id\n' + ''.join(s + '\n' for s in species))
                command += ['--representative-manifest', str(manifest), '--busco-mode', mode]
            result = subprocess.run(command, capture_output=True, text=True)
            rows = []
            if output.exists():
                with output.open(newline='') as handle:
                    rows = list(csv.DictReader(handle, delimiter='\t'))
            return result, rows

    def test_annotation_focal_only(self):
        result, rows = self.run_script('prepare_focal_genome_map.py', HEADER + FOCAL)
        self.assertEqual(result.returncode, 0, result.stderr)
        self.assertEqual([r['species_id'] for r in rows], ['A', 'B'])

    def test_annotation_ignores_outgroup(self):
        result, rows = self.run_script('prepare_focal_genome_map.py', HEADER + FOCAL + OUTGROUP)
        self.assertEqual(result.returncode, 0, result.stderr)
        self.assertEqual([r['species_id'] for r in rows], ['A', 'B'])

    def test_annotation_rejects_focal_proteomes(self):
        result, _ = self.run_script('prepare_focal_genome_map.py', (HEADER + FOCAL).replace('genome', 'proteome'))
        self.assertNotEqual(result.returncode, 0)
        self.assertIn('requires input_type=genome', result.stderr)

    def test_inference_genomes(self):
        result, rows = self.run_script('select_species_tree_inputs.py', HEADER + FOCAL + OUTGROUP)
        self.assertEqual(result.returncode, 0, result.stderr)
        self.assertEqual([r['input_type'] for r in rows], ['genome'] * 3)

    def test_inference_proteomes(self):
        result, rows = self.run_script('select_species_tree_inputs.py', (HEADER + FOCAL + OUTGROUP).replace('genome', 'proteome'))
        self.assertEqual(result.returncode, 0, result.stderr)
        self.assertEqual([r['input_type'] for r in rows], ['proteome'] * 3)

    def test_inference_mixed_types(self):
        table = HEADER + 'A\ta.fa\tgenome\tfocal\nB\tb.faa\tproteome\tfocal\nO\to.faa\tproteome\toutgroup\n'
        result, rows = self.run_script('select_species_tree_inputs.py', table)
        self.assertEqual(result.returncode, 0, result.stderr)
        self.assertEqual([r['input_type'] for r in rows], ['genome', 'proteome', 'proteome'])

    def test_annotation_with_proteome_outgroup(self):
        result, rows = self.run_script('prepare_focal_genome_map.py', HEADER + FOCAL + OUTGROUP.replace('genome', 'proteome'))
        self.assertEqual(result.returncode, 0, result.stderr)
        self.assertEqual([r['species_id'] for r in rows], ['A', 'B'])

    def test_preserved_proteome_override(self):
        result, rows = self.run_script('select_species_tree_inputs.py', HEADER + FOCAL + OUTGROUP, mode='proteome')
        self.assertEqual(result.returncode, 0, result.stderr)
        self.assertEqual([r['input_type'] for r in rows], ['proteome'] * 3)

    def test_preserved_genome_override(self):
        result, rows = self.run_script('select_species_tree_inputs.py', (HEADER + FOCAL + OUTGROUP).replace('genome', 'proteome'), mode='genome')
        self.assertEqual(result.returncode, 0, result.stderr)
        self.assertEqual([r['input_type'] for r in rows], ['genome'] * 3)

    def test_missing_outgroup(self):
        result, _ = self.run_script('select_species_tree_inputs.py', HEADER + FOCAL)
        self.assertNotEqual(result.returncode, 0)
        self.assertIn('exactly one outgroup', result.stderr)

    def test_extra_outgroup(self):
        result, _ = self.run_script('select_species_tree_inputs.py', HEADER + FOCAL + OUTGROUP + 'O2\to2.fa\tgenome\toutgroup\n')
        self.assertNotEqual(result.returncode, 0)
        self.assertIn('exactly one outgroup', result.stderr)

    def test_missing_focal_species(self):
        result, _ = self.run_script('select_species_tree_inputs.py', HEADER + FOCAL + OUTGROUP, species=('A', 'B', 'C'))
        self.assertNotEqual(result.returncode, 0)
        self.assertIn('Missing focal', result.stderr)

    def test_unused_focal_species(self):
        result, _ = self.run_script('select_species_tree_inputs.py', HEADER + FOCAL + OUTGROUP, species=('A',))
        self.assertNotEqual(result.returncode, 0)
        self.assertIn('Unused focal', result.stderr)

    def test_duplicate_species(self):
        result, _ = self.run_script('select_species_tree_inputs.py', HEADER + FOCAL + OUTGROUP + 'A\tx.fa\tgenome\tfocal\n')
        self.assertNotEqual(result.returncode, 0)
        self.assertIn('Duplicate species_id', result.stderr)

    def test_unsupported_input_type(self):
        result, _ = self.run_script('select_species_tree_inputs.py', (HEADER + FOCAL + OUTGROUP).replace('genome', 'unknown'))
        self.assertNotEqual(result.returncode, 0)
        self.assertIn('Unsupported input_type', result.stderr)

class InlineSelectionStructureTests(unittest.TestCase):
    def setUp(self):
        self.source = (ROOT / 'main.nf').read_text()
        self.selection = self.source.split('// BEGIN INPUT SELECTION', 1)[1].split('// END INPUT SELECTION', 1)[0]
        self.branches = self.source.split('// END INPUT SELECTION', 1)[1]

    def test_one_labelled_section_without_helper(self):
        self.assertEqual(self.source.count('// BEGIN INPUT SELECTION'), 1)
        self.assertEqual(self.source.count('// END INPUT SELECTION'), 1)
        self.assertNotIn('InputRoutes', self.source)
        self.assertFalse((ROOT / 'lib' / 'InputRoutes.groovy').exists())

    def test_each_decision_assigned_once(self):
        for name in ('input_mode', 'tree_source'):
            self.assertEqual(len(re.findall(r'^\s*' + name + r'\s*=', self.source, re.M)), 1)
            self.assertRegex(self.selection, r'\b' + name + r'\s*=')

    def test_branches_reuse_resolved_decisions(self):
        self.assertNotIn('params.input_mode', self.source)
        self.assertNotIn('if (params.species_tree)', self.branches)
        self.assertIn("if (tree_source == 'supplied')", self.branches)
        self.assertEqual(len(re.findall(r"^    if \(input_mode == 'annotation'\)", self.branches, re.M)), 1)
        self.assertIn("else if (input_mode == 'annotation')", self.branches)

    def test_validation_and_inference_are_exclusive_no_fallback(self):
        match = re.search(r"if \(tree_source == 'supplied'\) \{(.*?)\n    \}\n    else \{(.*?)\n    \}", self.branches, re.S)
        self.assertIsNotNone(match)
        supplied, inferred = match.groups()
        self.assertIn('VALIDATE_SUPPLIED_SPECIES_TREE(', supplied)
        self.assertNotIn('INFER_SPECIES_TREE(', supplied)
        self.assertIn('INFER_SPECIES_TREE(', inferred)
        self.assertNotIn('catch', supplied + inferred)

    def test_parameter_initialisation_fixes_retained(self):
        for name in ('id_species', 'outgroup_genes', 'famsa_args', 'iqtree_args'):
            self.assertRegex(self.source, r'params\.' + name + r'\s*=\s*params\.' + name + r'\s*\?:\s*null')
        self.assertIn('params.sequence_inputs', self.source)
        self.assertNotIn('params.species_inputs', self.source)
        self.assertNotIn('params.genome_inputs', self.source)

if __name__ == '__main__':
    unittest.main()
