import json
import os
from pathlib import Path
import tempfile
import unittest
from unittest.mock import patch
import sys
import numpy as np
from indexed_gam_pipeline.benchmark_partition_reader import merge_parts
from indexed_gam_pipeline.full_reader_comparison import logical_hash,run_builders

class FullComparisonTests(unittest.TestCase):
    def part(self,root,name,values):
        p=root/name;p.mkdir();n=len(values);summary=[]
        for i,v in enumerate(values):
            summary.append(dict(candidate_id=str(v),coverage=1,alt_count=1,ref_count=0,other_count=0,af=1.,selected_alignments=1,
                event_type='SNP',tensor_format_version='indexed-gam-candidate-v4',shard_index=i//2,index_within_shard=i%2))
        for i in range(0,n,2):
            np.save(p/f'shard_{i//2:05d}_data.npy',np.stack([np.full((7,200,101),v,dtype=np.int32) for v in values[i:i+2]]))
        (p/'variant_summary.ndjson').write_text(''.join(json.dumps(x)+'\n' for x in summary))
        (p/'manifest.json').write_text(json.dumps(dict(status='complete',tensor_format_version='indexed-gam-candidate-v4',tensors=n,shards=(n+1)//2,nodes=n)))
        return p
    def test_merge_partial_boundary_and_zero_part(self):
        with tempfile.TemporaryDirectory() as d:
            root=Path(d)
            for i,values in enumerate([([1,2,3],[4,5,6]),([],[1,2,3]),([1,2],[])]):
                parts=[self.part(root,f'{i}_{j}',v) for j,v in enumerate(values)]
                result=root/f'merged_{i}';merge_parts(parts,result,shard_size=2)
                self.assertEqual(logical_hash(parts),logical_hash([result]))
                if values[0] and len(values[0])>=2:
                    self.assertEqual((parts[0]/'shard_00000_data.npy').stat().st_ino,(result/'shard_00000_data.npy').stat().st_ino)
    def test_failed_worker_stops_peer_and_records_failure(self):
        with tempfile.TemporaryDirectory() as d:
            root=Path(d);seed=root/'seed';seed.write_bytes(b'cache')
            commands=[[sys.executable,'-c','import time; time.sleep(0.2); raise SystemExit(3)'],[sys.executable,'-c','import time; time.sleep(60)']]
            with patch('indexed_gam_pipeline.full_reader_comparison.build_command',side_effect=commands):
                with self.assertRaisesRegex(RuntimeError,'failed'):
                    run_builders({'seed_cache':str(seed)},root/'run',[root/'a',root/'b'],'python')
            status=json.loads((root/'run/status.json').read_text())
            self.assertEqual(status['status'],'failed')
            for proc in status['processes']:
                with self.assertRaises(ProcessLookupError):os.kill(proc['pid'],0)
