import subprocess
import tempfile
import unittest
from pathlib import Path
from unittest.mock import patch
from indexed_gam_pipeline.gam_reader import IndexedGam
from indexed_gam_pipeline.vg_reader import VgGam
from test_pipeline import fixture

class VgReaderTests(unittest.TestCase):
    def test_exact_nodes_order_duplicates_and_empty(self):
        with tempfile.TemporaryDirectory() as d:
            gam,_=fixture(d)
            def fake(command,stdout,stderr):
                stdout.write(gam.read_bytes())
                self.assertEqual(command[-1],'10:30')
                return subprocess.CompletedProcess(command,0)
            reader=VgGam(gam,vg='/unused/vg')
            with patch('indexed_gam_pipeline.vg_reader.subprocess.run',side_effect=fake) as run:
                actual=[a.SerializeToString() for a in reader.fetch({10,30})]
                expected=[a.SerializeToString() for a in IndexedGam(gam).fetch({10,30})]
                self.assertEqual(actual,expected)
                self.assertEqual(list(reader.fetch([])),[])
                self.assertEqual(run.call_count,1)
            self.assertEqual(reader.cache_stats['queries'],1)
    def test_command_failure_is_not_empty_result(self):
        with tempfile.TemporaryDirectory() as d:
            gam,_=fixture(d)
            with patch('indexed_gam_pipeline.vg_reader.subprocess.run',return_value=subprocess.CompletedProcess([],2)):
                with self.assertRaisesRegex(RuntimeError,'vg find failed'):list(VgGam(gam).fetch({10}))
    def test_reject_alternate_index_and_changed_input(self):
        with tempfile.TemporaryDirectory() as d:
            gam,_=fixture(d)
            with self.assertRaisesRegex(ValueError,'adjacent'):VgGam(gam,Path(d)/'other.gai')
            reader=VgGam(gam)
            with gam.open('ab') as f:f.write(b'changed')
            with self.assertRaisesRegex(ValueError,'changed'):list(reader.fetch({10}))
