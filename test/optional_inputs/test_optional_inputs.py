#!/usr/bin/env python3
import argparse
import ast
import csv
import json
import os
from pathlib import Path
import subprocess
import sys
import tempfile
import unittest

ROOT = Path(__file__).resolve().parents[2]

class OptionalInputs(unittest.TestCase):
    def setUp(self):
        self.temp = tempfile.TemporaryDirectory()
        self.addCleanup(self.temp.cleanup)
        self.work = Path(self.temp.name)

    def write(self, name, text):
        path = self.work / name
        path.write_text(text)
        return path

    def run_script(self, name, args):
        return subprocess.run([sys.executable, str(ROOT / 'bin' / name), *map(str, args)],
                              cwd=self.work, text=True, capture_output=True)

    def selector(self, supplied=False, table=None):
        proteins = self.write('proteins.faa', '>c1\nAAA\n>c2\nAAAAA\n>c3\nAAAA\n')
        manifest = self.write('manifest.tsv', 'sequence_id\tgene_key\tgene_id\ttranscript_id\tprotein_length\tspecies_id\n'
                              'c1\tS1:g1\tg1\tt1\t3\tS1\n'
                              'c2\tS1:g1\tg1\tt2\t5\tS1\n'
                              'c3\tS2:g2\tg2\tt3\t4\tS2\n')
        args = ['--proteins', proteins, '--manifest', manifest,
                '--selection-table-supplied', str(supplied).lower(),
                '--representatives', 'reps.faa', '--output-manifest', 'reps.tsv',
                '--report', 'report.tsv', '--validation-report', 'validation.tsv',
                '--template', 'template.tsv', '--status', 'status.tsv']
        if table is not None:
            args += ['--selection-table', table]
        return self.run_script('select_representative_isoforms.py', args)

    def test_automatic_selection_without_table(self):
        result = self.selector()
        self.assertEqual(result.returncode, 0, result.stderr)
        self.assertIn('complete', (self.work / 'status.tsv').read_text())
        self.assertEqual((self.work / 'reps.faa').read_text(), '>c2\nAAAAA\n>c3\nAAAA\n')
        self.assertFalse((self.work / 'template.tsv').exists())

    def test_manual_selection_with_table(self):
        table = self.write('manual.tsv', 'gene_id\ttranscript_id\ng1\tt1\n')
        result = self.selector(True, table)
        self.assertEqual(result.returncode, 0, result.stderr)
        self.assertIn('>c1\nAAA', (self.work / 'reps.faa').read_text())
        self.assertIn('manual', (self.work / 'reps.tsv').read_text())

    def test_incomplete_manual_table_still_generates_template(self):
        table = self.write('manual.tsv', 'gene_id\ttranscript_id\n')
        result = self.selector(True, table)
        self.assertEqual(result.returncode, 0, result.stderr)
        self.assertIn('manual_selection_required', (self.work / 'status.tsv').read_text())
        self.assertTrue((self.work / 'template.tsv').is_file())
        self.assertFalse((self.work / 'reps.faa').exists())

    def test_manual_mode_without_table_is_error(self):
        result = self.selector(True)
        self.assertNotEqual(result.returncode, 0)
        self.assertIn('--selection-table is required', result.stderr)

    def test_missing_supplied_selection_table_is_error(self):
        result = self.selector(True, self.work / 'missing.tsv')
        self.assertNotEqual(result.returncode, 0)
        self.assertFalse((self.work / 'reps.faa').exists())

    def generax(self, outgroup=None):
        fixtures = {
            'gene.nwk': '(c1:0.1,c2:0.1);\n',
            'iqtree.txt': 'Best-fit model according to BIC: LG+G4\n',
            'alignment.faa': '>c1\nAAA\n>c2\nAAA\n',
            'species.nwk': '(S1:0.1,S2:0.1);\n',
            'mapping.tsv': 'sequence_id\tspecies_id\nc1\tS1\nc2\tS2\n',
            'representatives.tsv': 'sequence_id\tgene_id\tspecies_id\nc1\tg1\tS1\nc2\tg2\tS2\n',
        }
        for name, text in fixtures.items():
            self.write(name, text)
        args = ['--gene-tree', 'gene.nwk', '--iqtree-report', 'iqtree.txt',
                '--alignment', 'alignment.faa', '--species-tree', 'species.nwk',
                '--gene-to-species', 'mapping.tsv', '--representative-manifest', 'representatives.tsv',
                '--families-out', 'families.txt', '--mapping-out', 'generax_mapping.tsv',
                '--generax-tree-out', 'translated.nwk', '--generax-alignment-out', 'translated.faa',
                '--summary-out', 'summary.tsv']
        if outgroup is not None:
            args += ['--outgroup-genes', outgroup]
        return self.run_script('prepare_generax_inputs.py', args)

    def summary(self):
        with (self.work / 'summary.tsv').open() as handle:
            return dict(list(csv.reader(handle, delimiter='\t'))[1:])

    def test_generax_preparation_without_outgroups(self):
        result = self.generax()
        self.assertEqual(result.returncode, 0, result.stderr)
        self.assertEqual(self.summary()['outgroup_genes'], '0')
        self.assertIn('g1', (self.work / 'translated.nwk').read_text())

    def test_generax_preparation_with_outgroups(self):
        outgroup = self.write('outgroups.txt', 'g1\n')
        result = self.generax(outgroup)
        self.assertEqual(result.returncode, 0, result.stderr)
        self.assertEqual(self.summary()['outgroup_genes'], '1')

    def test_generax_rejects_unknown_outgroup(self):
        result = self.generax(self.write('outgroups.txt', 'unknown\n'))
        self.assertNotEqual(result.returncode, 0)
        self.assertIn('Outgroup gene IDs absent', result.stderr)

    def test_generax_rejects_missing_supplied_outgroup_file(self):
        result = self.generax(self.work / 'missing.txt')
        self.assertNotEqual(result.returncode, 0)

    def runner(self, outgroup=None, registry=None, options=''):
        fake_core = self.write('fake_core.py',
            "import json, sys\nfrom pathlib import Path\n"
            "args = sys.argv[1:]\nPath('core_args.json').write_text(json.dumps(args))\n"
            "tag = args[args.index('--revtag') + 1]\n"
            "Path(f'revlog_{tag}.csv').write_text('uid,' + tag + '\\ng1,1-1\\n')\n"
            "Path('named_reconaming.nhx').write_text('(g1,g2);\\n')\n"
            "Path('named_reconaming.csv').write_text('gene_id,name\\ng1,1-1\\n')\n"
            "Path('temporary_partition.csv').write_text('gene_id,group\\ng1,1\\n')\n")
        tree = self.write('tree.nhx', '(g1,g2);\n')
        ctl = self.write('options.ctl', options)
        args = ['--core-script', fake_core, '--tree', tree, '--revision-tag', 'smoke', '--options', ctl]
        if outgroup is not None:
            args += ['--outgroup-genes', outgroup]
        if registry is not None:
            args += ['--revision-registry', registry]
        return self.run_script('run_reconaming.py', args)

    def core_arguments(self):
        return json.loads((self.work / 'core_args.json').read_text())

    def test_runner_omits_absent_inputs(self):
        result = self.runner()
        self.assertEqual(result.returncode, 0, result.stderr)
        args = self.core_arguments()
        self.assertNotIn('--revlog', args)
        self.assertNotIn('--root-by', args)
        self.assertTrue((self.work / 'revision_registry_smoke.csv').exists())

    def test_runner_passes_supplied_inputs(self):
        outgroup = self.write('outgroups.txt', 'g1\n')
        registry = self.write('registry.csv', 'uid,old\ng1,1-1\n')
        result = self.runner(outgroup, registry)
        self.assertEqual(result.returncode, 0, result.stderr)
        args = self.core_arguments()
        self.assertEqual(args[args.index('--revlog') + 1], str(registry))
        self.assertEqual(args[args.index('--root-by') + 1], 'g1')
        self.assertEqual(registry.read_text(), 'uid,old\ng1,1-1\n')

    def test_explicit_ctl_root_without_outgroup_file(self):
        result = self.runner(options='root_by=g2\n')
        self.assertEqual(result.returncode, 0, result.stderr)
        args = self.core_arguments()
        self.assertEqual(args[args.index('--root-by') + 1], 'g2')

    def test_runner_rejects_missing_registry(self):
        result = self.runner(registry=self.work / 'missing.csv')
        self.assertNotEqual(result.returncode, 0)
        self.assertIn('--revision-registry file does not exist', result.stderr)
        self.assertFalse((self.work / 'core_args.json').exists())

    def test_runner_rejects_missing_outgroup(self):
        result = self.runner(outgroup=self.work / 'missing.txt')
        self.assertNotEqual(result.returncode, 0)
        self.assertIn('--outgroup-genes file does not exist', result.stderr)
        self.assertFalse((self.work / 'core_args.json').exists())

    def core_loader(self):
        try:
            import pandas as pd
        except ImportError:
            self.skipTest('pandas unavailable: isolated core registry-loader test')
        tree = ast.parse((ROOT / 'bin/reconaming_core.py').read_text())
        node = next(n for n in tree.body if isinstance(n, ast.FunctionDef) and n.name == '_load_revlog')
        namespace = {'pd': pd, 'os': os}
        exec(compile(ast.Module(body=[node], type_ignores=[]), 'core_loader', 'exec'), namespace)
        return namespace['_load_revlog']

    def test_core_loader_initialises_empty_registry(self):
        registry = self.core_loader()(None)
        self.assertTrue(registry.empty)
        self.assertEqual(len(registry.columns), 0)

    def test_core_loader_reads_supplied_registry(self):
        path = self.write('registry.csv', 'uid,old\ng1,1-1\n')
        registry = self.core_loader()(str(path))
        self.assertEqual(registry.loc['g1', 'old'], '1-1')

    def test_core_revlog_argument_is_optional(self):
        tree = ast.parse((ROOT / 'bin/reconaming_core.py').read_text())
        call = next(n for n in ast.walk(tree) if isinstance(n, ast.Call)
                    and isinstance(n.func, ast.Attribute) and n.func.attr == 'add_argument'
                    and any(isinstance(a, ast.Constant) and a.value == '--revlog' for a in n.args))
        self.assertFalse(any(k.arg == 'required' and isinstance(k.value, ast.Constant)
                             and k.value.value is True for k in call.keywords))

    def test_main_has_no_dummy_file_dependencies(self):
        main = (ROOT / 'main.nf').read_text()
        for param, filename in [('selection_table', 'empty_selection_table.tsv'),
                                ('outgroup_genes', 'empty_outgroup_genes.txt'),
                                ('revision_registry', 'empty_revision_registry.csv')]:
            self.assertNotIn(filename, main)
            self.assertIn(f'{param} = params.{param} ? file(params.{param}, checkIfExists: true) : []', main)

    def test_config_declares_previously_undefined_defaults(self):
        config = (ROOT / 'nextflow.config').read_text()
        for name in ['revision_registry', 'revision_tag', 'reconaming_opt', 'reconaming_container', 'gff_output']:
            self.assertRegex(config, rf'(?m)^\s*{name}\s*=')

    def test_modules_use_conditional_optional_arguments(self):
        selector = (ROOT / 'modules/local/select_representative_isoforms.nf').read_text()
        self.assertIn('selection_table_supplied ? "--selection-table', selector)
        self.assertIn('${selection_args}', selector)
        for filename in ['reconcile_gene_tree.nf', 'run_reconaming.nf']:
            text = (ROOT / 'modules/local' / filename).read_text()
            self.assertIn('outgroup_genes ? "--outgroup-genes', text)
            self.assertIn('${outgroup_args}', text)
        runner = (ROOT / 'modules/local/run_reconaming.nf').read_text()
        self.assertIn('revision_registry ? "--revision-registry', runner)
        self.assertIn('${registry_args}', runner)

if __name__ == '__main__':
    unittest.main()
