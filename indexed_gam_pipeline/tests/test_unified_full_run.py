import json
import os
from pathlib import Path
import tempfile
import unittest
from unittest.mock import patch
from indexed_gam_pipeline.unified_full_run import partition_nodes, MemoryRecorder
from indexed_gam_pipeline.full_reader_comparison import run_builders


class UnifiedFullTests(unittest.TestCase):
    def test_twenty_disjoint_contiguous_partitions(self):
        with tempfile.TemporaryDirectory() as d:
            root=Path(d);source=root/'nodes';source.write_text(''.join(f'{i*3}\n' for i in range(1,104)))
            parts=partition_nodes(source,root,20,103,2)
            self.assertEqual(len(parts),20)
            self.assertEqual(''.join(Path(p['nodes_file']).read_text() for p in parts),source.read_text())
            self.assertEqual(sum(p['nodes'] for p in parts),103)
            self.assertTrue(all(parts[i]['last_node']<parts[i+1]['first_node'] for i in range(19)))
            self.assertLessEqual(max(p['nodes'] for p in parts)-min(p['nodes'] for p in parts),1)

    def test_reject_bad_discovery(self):
        for contents,total in [('1\n1\n',2),('2\n1\n',2),('1\n',2),('1\n2\n3\n',2)]:
            with tempfile.TemporaryDirectory() as d:
                root=Path(d);source=root/'nodes';source.write_text(contents)
                with self.assertRaises(ValueError): partition_nodes(source,root,1,total)

    def test_twenty_builders_share_readonly_index(self):
        from test_pipeline import fixture
        from test_graph_index import fixture_index
        with tempfile.TemporaryDirectory() as d:
            root=Path(d);gam,_=fixture(d)
            graph=root/'graph.sqlite';fixture_index(graph,[(n,'AAAAAA',86) for n in [10,20,30,1000]])
            files=[]
            for i in range(20):
                p=root/f'nodes_{i}';p.write_text('10\n20\n30\n1000\n');files.append(p)
            cfg=dict(config=dict(gam=str(gam),index=str(gam)+'.gai',unified_graph_index=str(graph),vg='/unused'),
                     prior_run=str(root),params=dict(gam_cache_mb_per_process=1,max_tensors=2))
            # Use production commands/controller, with fixture's lower support threshold.
            from indexed_gam_pipeline.full_reader_comparison import build_command
            def command(*args):
                cmd=build_command(*args);cmd[cmd.index('--min-variants')+1]='1';return cmd
            with patch.dict(os.environ,dict(OPENBLAS_NUM_THREADS='1',OMP_NUM_THREADS='1',MKL_NUM_THREADS='1')), \
                 patch('indexed_gam_pipeline.full_reader_comparison.build_command',side_effect=command), \
                 MemoryRecorder(root,lambda:'test',interval=.1) as memory:
                parts=run_builders(cfg,root/'run',files,'python')
            self.assertEqual(len(parts),20)
            self.assertTrue(all(json.loads((p/'manifest.json').read_text())['tensors']>0 for p in parts))
            self.assertFalse(list((root/'run').glob('gbwt_*.sqlite')))
            self.assertGreater(memory.peak_rss_kib,0)
