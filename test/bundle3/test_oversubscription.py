#!/usr/bin/env python3
"""Static wiring checks; Nextflow previews and smoke runs remain necessary."""
from pathlib import Path
import unittest

ROOT = Path(__file__).resolve().parents[2]
CONFIG = (ROOT / 'nextflow.config').read_text()
MAIN = (ROOT / 'main.nf').read_text()
MODULE = (ROOT / 'modules/local/reconcile_gene_tree.nf').read_text()
README = (ROOT / 'README.md').read_text()

class OversubscriptionWiringTests(unittest.TestCase):
    def test_default_enabled(self):
        self.assertRegex(CONFIG, r'(?m)^\s*generax_oversubscribe\s*=\s*true\s*$')

    def test_validation_and_normalisation_before_dispatch(self):
        self.assertIn('generax_oversubscribe_raw instanceof Boolean', MAIN)
        self.assertIn('generax_oversubscribe_raw instanceof CharSequence', MAIN)
        self.assertIn("generax_oversubscribe_text in ['true', 'false']", MAIN)
        self.assertIn('generax_oversubscribe_raw?.toString()?.trim()?.toLowerCase()', MAIN)
        normalise = "def generax_oversubscribe = generax_oversubscribe_text == 'true'"
        self.assertIn(normalise, MAIN)
        self.assertLess(MAIN.index(normalise), MAIN.index('    RECONCILE_GENE_TREE('))
        self.assertIn('must be true or false (Boolean or CLI string)', MAIN)

    def test_explicit_normalised_module_input(self):
        self.assertRegex(MAIN, r'params\.generax_container,\s*generax_oversubscribe\s*\)')
        self.assertIn('val generax_oversubscribe', MODULE)

    def test_conditional_flag_and_unchanged_rank_allocation(self):
        self.assertIn("def mpi_oversubscribe = generax_oversubscribe ? '--oversubscribe' : ''", MODULE)
        command = next(line.strip() for line in MODULE.splitlines() if line.strip().startswith('mpiexec '))
        self.assertTrue(command.startswith('mpiexec ${mpi_oversubscribe} -np ${task.cpus} generax '))
        self.assertIn('cpus { params.process_cpus.RECONCILE_GENE_TREE }', MODULE)
        self.assertEqual(MODULE.count('--oversubscribe'), 1)

    def test_documentation_and_removed_obsolete_parameter(self):
        self.assertNotIn('--generax_cpus', README)
        self.assertIn('--generax_oversubscribe false', README)
        self.assertIn('not currently developed or validated', README)

if __name__ == '__main__':
    unittest.main()
