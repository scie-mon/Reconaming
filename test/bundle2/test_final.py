import importlib.util, tempfile
from pathlib import Path
import unittest
ROOT=Path(__file__).resolve().parents[2]
spec=importlib.util.spec_from_file_location('core',ROOT/'bin/run_reconaming.py');core=importlib.util.module_from_spec(spec);spec.loader.exec_module(core)
class TestFinal(unittest.TestCase):
 def test_cpu_map(self):
  t=(ROOT/'nextflow.config').read_text();b=t.split('process_cpus = [',1)[1].split(']',1)[0]
  for k,v in {'PREPARE_PROTEINS':1,'RUN_BUSCO':16,'RUN_BUSCO_PHYLOGENOMICS':64,'INFER_SUPERMATRIX_TREE':64,'ALIGN_PROTEINS':64,'INFER_GENE_TREE':64,'RECONCILE_GENE_TREE':64}.items():self.assertRegex(b,rf'{k}\s*:\s*{v}')
 def test_iqthread_setting(self):self.assertIn("iqtree_threads = 'AUTO'",(ROOT/'nextflow.config').read_text())
 def test_core_only(self):
  with tempfile.TemporaryDirectory() as d:
   p=Path(d)/'x.ctl';p.write_text('annotation=\n')
   with self.assertRaises(SystemExit):core.read_ctl(p)
 def test_shipped_core_keys(self):self.assertEqual(set(core.read_ctl(ROOT/'conf/reconaming_opt.ctl')),core.CORE_OPTIONS)
 def test_no_auto_launcher(self):
  for p in ['main.nf','nextflow.config']:
   t=(ROOT/p).read_text();self.assertNotIn('resolved_local',t);self.assertNotIn('cpu_resources_resolved',t)
if __name__=='__main__':unittest.main()
