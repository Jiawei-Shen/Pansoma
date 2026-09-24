import json
from pathlib import Path
import sys
import tempfile
import unittest

from indexed_gam_pipeline.dynamic_tensor_run import execute_queue
from indexed_gam_pipeline.unified_full_run import partition_nodes


class DynamicTests(unittest.TestCase):
    def test_dynamic_refill_and_fresh_processes(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            script = root/'worker.py'
            script.write_text('import os,sys,time,json\nfrom pathlib import Path\np=Path(sys.argv[1])\nstart=time.time()\ntime.sleep(float(sys.argv[2]))\np.write_text(json.dumps(dict(pid=os.getpid(),start=start,end=time.time())))\n')
            commands = [[sys.executable, str(script), str(root/f'{i}.json'), str(delay)]
                        for i, delay in enumerate([.8,.05,.05,.05])]
            result = execute_queue(commands, root/'output', 2, .01)
            rows = [json.loads((root/f'{i}.json').read_text()) for i in range(4)]
            self.assertEqual(result['completed_tasks'], 4)
            self.assertEqual(len({r['pid'] for r in rows}), 4)
            self.assertLess(rows[2]['start'], rows[0]['end'])
            events = sorted([(r['start'],1) for r in rows]+[(r['end'],-1) for r in rows])
            active = 0
            for _, delta in events:
                active += delta
                self.assertLessEqual(active, 2)

    def test_failure_cancels_active_and_does_not_start_pending(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            commands = [[sys.executable,'-c','import time;time.sleep(.1);raise SystemExit(7)'],
                        [sys.executable,'-c','import time;time.sleep(30)'],
                        [sys.executable,'-c','raise SystemExit(0)']]
            with self.assertRaisesRegex(RuntimeError, 'code 7'):
                execute_queue(commands, root/'output', 2, .01)
            status = json.loads((root/'output/status.json').read_text())
            self.assertEqual(status['status'], 'failed')
            self.assertEqual([t['status'] for t in status['tasks']], ['failed','cancelled','pending'])

    def test_real_build_tasks_match_unsplit_tensor(self):
        import os
        import numpy as np
        from unittest.mock import patch
        from indexed_gam_pipeline.tests.test_pipeline import fixture
        from indexed_gam_pipeline.tests.test_graph_index import fixture_index
        from indexed_gam_pipeline.dynamic_tensor_run import task
        with tempfile.TemporaryDirectory() as tmp:
            root=Path(tmp); gam,_=fixture(tmp)
            graph=root/'graph.sqlite'
            fixture_index(graph, [(n,'AAAAAA',86) for n in [10,20,30,1000]])
            parts=[]
            for i,nodes in enumerate([[10,20,30,1000],[10,20],[30,1000]]):
                path=root/f'nodes_{i}.txt';path.write_text(''.join(f'{n}\n' for n in nodes))
                parts.append(dict(nodes_file=str(path)))
            config=dict(config=dict(gam=str(gam),index=str(gam)+'.gai',
                        unified_graph_index=str(graph),vg='/unused'), prior_run=str(root),
                        params=dict(gam_cache_mb_per_process=1),parts=parts)
            (root/'config.json').write_text(json.dumps(config))
            with patch.dict(os.environ,dict(OPENBLAS_NUM_THREADS='1',OMP_NUM_THREADS='1')):
                commands=[[sys.executable,'-m','indexed_gam_pipeline.dynamic_tensor_run',
                           'task','--root',str(root),'--index',str(i)] for i in range(3)]
                result=execute_queue(commands,root/'python',2,.01)
            self.assertEqual(result['completed_tasks'],3)
            outputs=[]
            for i in range(3):
                dest=root/'python'/f'task_{i:04d}'
                report=json.loads((dest/'validation_report.json').read_text())
                self.assertTrue(report['passed'])
                outputs.append([np.load(f) for f in sorted(dest.glob('shard_*_data.npy'))])
            self.assertTrue(outputs[0])
            np.testing.assert_array_equal(np.concatenate(outputs[0]),np.concatenate(outputs[1]+outputs[2]))

    def test_split_outputs_apply_independent_af_thresholds(self):
        import numpy as np
        import pysam
        from test_pipeline import alignment, vi
        from test_graph_index import fixture_index
        from indexed_gam_pipeline.gam_reader import build_index
        from indexed_gam_pipeline.dynamic_tensor_run import task
        from indexed_gam_pipeline.tensor_storage import STORAGE_VERSION
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp); gam = root/'tiny.gam'
            with pysam.BGZFile(str(gam), 'wb') as stream:
                for node, alt_count in [(10, 3), (20, 3), (30, 4), (40, 3), (50, 4)]:
                    for i in range(50):
                        a = alignment((node,), mutation=(node == 10 and i < alt_count))
                        a.name = f'{node}-{i}'
                        if node != 10 and i < alt_count:
                            m = a.path.mapping[0]; del m.edit[:]
                            m.edit.add(from_length=2, to_length=2)
                            if node < 40:
                                m.edit.add(from_length=1, to_length=0)
                                m.edit.add(from_length=3, to_length=3)
                                a.sequence = 'AAAAA'
                            else:
                                m.edit.add(from_length=0, to_length=1, sequence='T')
                                m.edit.add(from_length=4, to_length=4)
                                a.sequence = 'AATAAAA'
                            a.quality = bytes([30]*len(a.sequence))
                        raw = a.SerializeToString()
                        stream.write(vi(2)+vi(3)+b'GAM'+vi(len(raw))+raw)
                        stream.flush()
            build_index(gam, str(gam)+'.gai')
            graph = root/'graph.sqlite'; fixture_index(graph, [(n,'AAAAAA',331) for n in (10,20,30,40,50)])
            nodes = root/'nodes.txt'; nodes.write_text('10\n20\n30\n40\n50\n')
            config = dict(config=dict(gam=str(gam),index=str(gam)+'.gai',unified_graph_index=str(graph),vg='/unused'),
                prior_run=str(root), params=dict(gam_cache_mb_per_process=1,decode_policy='full'),
                parts=[dict(nodes_file=str(nodes))], tensor_storage_version=STORAGE_VERSION,
                variant_outputs=[dict(folder='SNV',variant_type='snp',min_af=.06),
                                 dict(folder='INDEL',variant_type='indel',min_af=.08)])
            (root/'config.json').write_text(json.dumps(config))
            task(root, 0)
            for folder, expected_node, kind, threshold in [('SNV',10,'SNP',.06),('INDEL',30,'DEL',.08)]:
                dest = root/folder/'task_0000'
                records = [json.loads(line) for line in (dest/'variant_summary.ndjson').read_text().splitlines()]
                self.assertEqual([(m['node_id'],m['event_type']) for m in records],
                                 [(10,'SNP')] if folder == 'SNV' else [(30,'DEL'),(50,'INS')])
                self.assertEqual(records[0]['af'], threshold)
                x = np.load(dest/'shard_00000_data.npy')
                self.assertEqual(x.dtype,np.int8)
                self.assertEqual(int(x[:,6].max()),82)
                manifest = json.loads((dest/'manifest.json').read_text())
                self.assertEqual(manifest['parameters']['min_af'],threshold)

            # Compare one shared pass with two independent builds, with tiny shards
            # to exercise full-shard flushes, final shards and per-type numbering.
            from unittest.mock import patch
            from contextlib import redirect_stdout
            import io
            from indexed_gam_pipeline import run, build_v2
            from indexed_gam_pipeline.full_reader_comparison import build_command
            from indexed_gam_pipeline.gam_reader import IndexedGam
            def invoke(dest, extra):
                cmd = build_command(config, str(nodes), dest, None, 'python', 1)
                cmd[cmd.index('--shard-size')+1] = '1'
                with patch.object(sys,'argv',cmd[1:]+['--decode-mode','full']+extra), redirect_stdout(io.StringIO()):
                    run.main()
            fetch = IndexedGam.fetch
            with patch.object(build_v2,'decode_alignment',wraps=build_v2.decode_alignment) as decode, \
                 patch.object(IndexedGam,'fetch',autospec=True,side_effect=fetch) as read:
                invoke(root/'shared', ['--snv-output',str(root/'shared_snv'),
                    '--indel-output',str(root/'shared_indel'),'--snv-min-af','.06','--indel-min-af','.08'])
                self.assertEqual(decode.call_count,250)
                self.assertEqual(read.call_count,1)
            with patch.object(build_v2,'decode_alignment',wraps=build_v2.decode_alignment) as decode, \
                 patch.object(IndexedGam,'fetch',autospec=True,side_effect=fetch) as read:
                for kind,af in [('snp','.06'),('indel','.08')]:
                    invoke(root/('baseline_'+kind),['--variant-type',kind,'--min-af',af])
                self.assertEqual(decode.call_count,500)
                self.assertEqual(read.call_count,2)
            for folder,kind in [('snv','snp'),('indel','indel')]:
                actual=root/('shared_'+folder);expected=root/('baseline_'+kind)
                self.assertEqual((actual/'variant_summary.ndjson').read_bytes(),
                                 (expected/'variant_summary.ndjson').read_bytes())
                actual_files=sorted(actual.glob('shard_*_data.npy'))
                self.assertEqual(len(actual_files),len(list(expected.glob('shard_*_data.npy'))))
                for file in actual_files:
                    self.assertEqual(file.read_bytes(),(expected/file.name).read_bytes())
            self.assertFalse(list((root/'shared').glob('shard_*_data.npy')))
            shared = json.loads((root/'shared/manifest.json').read_text())
            self.assertEqual(shared['tensors_by_type'],dict(SNV=1,INDEL=2))
            for name,extra,expected in [('empty', ['--snv-min-af','1','--indel-min-af','1'],0),
                                        ('limited',['--max-tensors','1'],1),
                                        ('workers',['--workers','2'],3)]:
                invoke(root/name, ['--snv-output',str(root/(name+'_SNV')),
                    '--indel-output',str(root/(name+'_INDEL')),'--snv-min-af','.06','--indel-min-af','.08']+extra)
                coordinator = json.loads((root/name/'manifest.json').read_text())
                self.assertEqual(coordinator['status'],'complete')
                self.assertEqual(coordinator['tensors'],expected)
                from indexed_gam_pipeline.full_run import validate_shards
                reports=[validate_shards(root/(name+'_'+kind),1) for kind in ('SNV','INDEL')]
                self.assertEqual(sum(r['tensors'] for r in reports),expected)
                if name == 'workers':
                    for kind in ('SNV','INDEL'):
                        actual=root/(name+'_'+kind);reference=root/('shared_'+kind.lower())
                        self.assertEqual((actual/'variant_summary.ndjson').read_bytes(),
                                         (reference/'variant_summary.ndjson').read_bytes())
                        for file in actual.glob('shard_*_data.npy'):
                            self.assertEqual(file.read_bytes(),(reference/file.name).read_bytes())


    def test_512_tasks_exact_coverage(self):
        with tempfile.TemporaryDirectory() as tmp:
            root=Path(tmp); source=root/'input'
            source.write_text(''.join(f'{i*3}\n' for i in range(1,10001)))
            parts=partition_nodes(source,root,512,10000)
            self.assertEqual(len(parts),512)
            self.assertEqual(''.join(Path(p['nodes_file']).read_text() for p in parts),source.read_text())


if __name__ == '__main__':
    unittest.main()
