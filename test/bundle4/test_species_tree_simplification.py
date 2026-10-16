#!/usr/bin/env python3
"""Static wiring checks only; these do not execute Nextflow, BUSCO or IQ-TREE."""
from pathlib import Path
import unittest

ROOT = Path(__file__).resolve().parents[2]
CONFIG = (ROOT / 'nextflow.config').read_text()
MODULE = (ROOT / 'modules/local/infer_species_tree.nf').read_text()
README = (ROOT / 'README.md').read_text()

class SpeciesTreeSimplificationTests(unittest.TestCase):
    def test_obsolete_public_parameters_removed(self):
        for text in (CONFIG, MODULE, README):
            for key in ('busco_auto_lineage', 'species_tree_min_taxa'):
                self.assertNotIn(key, text)

    def test_lineage_normalised_before_selection(self):
        self.assertIn('def lineage = params.busco_lineage?.toString()?.trim()', MODULE)
        self.assertIn('def lineage_arg = lineage ? "-l ${quoted_lineage}" : \'--auto-lineage\'', MODULE)

    def test_explicit_lineage_shell_quoted(self):
        expected = r'''def quoted_lineage = "'" + (lineage ?: '').replace("'", "'\"'\"'") + "'"'''
        self.assertIn(expected, MODULE)

    def test_busco_failure_not_masked_by_fallback(self):
        command = next(line.strip() for line in MODULE.splitlines() if line.strip().startswith('busco -i '))
        self.assertIn('${lineage_arg}', command)
        self.assertIn('-c ${task.cpus}', command)
        self.assertNotIn('||', command)
        self.assertEqual(MODULE.count('--auto-lineage'), 1)
        self.assertIn("shell = ['/bin/bash', '-euo', 'pipefail']", CONFIG)

    def test_only_supermatrix_requested(self):
        command = next(line.strip() for line in MODULE.splitlines() if line.strip().startswith('BUSCO_phylogenomics.py '))
        self.assertIn('--supermatrix_only', command)
        self.assertNotIn('--min_species_gene_tree', command)
        self.assertNotIn('--gene_tree_program', command)
        self.assertIn('-t ${task.cpus}', command)
        self.assertIn('INFER_SUPERMATRIX_TREE(RUN_BUSCO_PHYLOGENOMICS.out.supermatrix, RUN_BUSCO_PHYLOGENOMICS.out.partitions)', MODULE)

    def test_documentation_describes_lineage_and_supermatrix(self):
        self.assertIn('omitted or blank selects automatic lineage selection', README)
        self.assertIn('invalid explicit lineage stops the run', README)
        self.assertIn('individual BUSCO marker trees are not generated', README)

if __name__ == '__main__':
    unittest.main()
