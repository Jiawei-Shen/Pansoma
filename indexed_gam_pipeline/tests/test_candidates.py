"""Candidate-v2 semantics, synthetic indexed integration, and output contracts."""
import argparse
from contextlib import redirect_stdout
import io
import json
from pathlib import Path
import tempfile
import unittest

import numpy as np

from indexed_gam_pipeline import vg_pb2
from indexed_gam_pipeline.candidates import (BASES, OPS, Candidate, decode_alignment,
                                            overlap, make_tensor, rc)
from indexed_gam_pipeline.run import build
from test_pipeline import fixture


def alignment(specs, sequences, name="same-name"):
    a = vg_pb2.Alignment(name=name, mapping_quality=60)
    for nid, offset, reverse, edits in specs:
        m = a.path.mapping.add()
        m.position.node_id, m.position.offset, m.position.is_reverse = nid, offset, reverse
        ref = rc(sequences[nid]) if reverse else sequences[nid]
        cursor = offset
        for f, t, seq in edits:
            m.edit.add(from_length=f, to_length=t, sequence=seq)
            a.sequence += seq if seq else ref[cursor:cursor+f] if t else ""
            cursor += f
    a.quality = bytes([30]*len(a.sequence))
    return a


def read(specs, sequences):
    return decode_alignment(alignment(specs, sequences), sequences)[0]


def eligible(candidate, reads):
    return [(r, *overlap(r, candidate, 10)) for r in reads if overlap(r, candidate, 10)]


class CandidatesTest(unittest.TestCase):
    def test_window_edit_bp_excludes_distant_edits_and_counts_indel_bases(self):
        seq = {1: 'A'*60}
        candidate = Candidate(1, 30, 'A', 'T', 'SNP')
        specs = [
            [(10,10,'T'*10),(20,20,''),(1,1,'T'),(29,29,'')],
            [(28,28,''),(0,2,'GG'),(2,2,''),(1,1,'T'),(29,29,'')],
            [(28,28,''),(2,0,''),(1,1,'T'),(29,29,'')],
            [(29,29,''),(1,1,'C'),(1,1,'T'),(29,29,'')],
        ]
        reads = [decode_alignment(alignment([(1,0,False,edits)],seq,name=str(i)),seq)[0]
                 for i, edits in enumerate(specs)]
        _, meta = make_tensor(candidate, eligible(candidate, reads), width=11, debug=True)
        scores = {row['read_name']:row['window_mismatch_bp'] for row in meta['rows']}
        self.assertEqual(scores, {'0':1,'1':3,'2':3,'3':2})
        self.assertEqual(meta['window_mismatch_bp'], [3,3,2,1])

    def test_group_before_uniform_200_sampling_and_keep_seven_channels(self):
        seq = {1:'AC',2:'GG',3:'TT'}
        reads = []
        for i in range(401):
            branch = 3 if i % 2 else 2
            # Branch 3 has the highest score, so it must precede path 2.
            first = [(1,1,'A'),(1,1,'')] if branch == 3 else [(2,2,'')]
            a = alignment([(branch,0,False,first),(1,0,False,[(1,1,''),(1,1,'T')])],
                          seq,name=f'read-{i:04d}')
            reads.append(decode_alignment(a,seq)[0])
        c = Candidate(1,1,'C','T','SNP')
        e = eligible(c, reads)
        full, all_meta = make_tensor(c,e,rows=401,width=9,debug=True,node_walk_counts={1:90,2:3,3:7})
        x, meta = make_tensor(c,list(reversed(e)),rows=200,width=9,debug=True,
                              node_walk_counts={1:90,2:3,3:7})
        indices = [i*400//199 for i in range(200)]
        np.testing.assert_array_equal(x, full[:,indices,:])
        self.assertEqual(meta['selected_grouped_ranks'],indices)
        self.assertEqual([g['path'][0]['node_id'] for g in meta['row_groups']],[3,2])
        self.assertEqual(meta['coverage'],401)
        self.assertEqual(meta['alt_count'],401)
        self.assertEqual(meta['selected_alignments'],200)
        self.assertEqual([r['record_sha256'] for r in meta['rows']],
                         [all_meta['rows'][i]['record_sha256'] for i in indices])

    def test_edits_and_indel_limits(self):
        seq = {1: "A"*110}
        a = alignment([(1, 0, False, [(1, 1, "T"), (0, 50, "C"*50),
            (50, 0, ""), (0, 51, "G"*51), (51, 0, ""), (2, 3, "TGC")])], seq)
        r, unsupported = decode_alignment(a, seq)
        self.assertEqual([(o.candidate.kind, max(len(o.candidate.ref), len(o.candidate.alt)))
                          for o in r.observations], [("SNP", 1), ("INS", 50), ("DEL", 50)])
        self.assertEqual([e["reason"] for e in unsupported],
                         ["indel_exceeds_limit", "indel_exceeds_limit", "complex_replacement_not_supported"])
        for obs in r.observations:
            tensor, meta = make_tensor(obs.candidate, eligible(obs.candidate, [r]))
            self.assertGreaterEqual(meta["candidate_columns"][1]-meta["candidate_columns"][0],
                                    max(len(obs.candidate.ref), len(obs.candidate.alt)))
        # A GAM match is not rediscovered as SNPs from sequence comparison.
        a = alignment([(1, 0, False, [(4, 4, "")])], seq)
        a.sequence = "TTTT"
        self.assertFalse(decode_alignment(a, seq)[0].observations)

    def test_multinode_branches_slices_and_per_row_reference(self):
        seq = {1: "AACG", 2: "TT", 3: "GC", 4: "AT"}
        a = read([(2, 0, False, [(2, 2, "")]), (1, 1, False, [(1, 1, ""), (1, 1, "T"), (1, 1, "")]),
                  (4, 0, False, [(2, 2, "")])], seq)
        b = read([(3, 0, False, [(2, 2, "")]), (1, 1, False, [(3, 3, "")]),
                  (4, 0, False, [(2, 2, "")])], seq)
        c = a.observations[0].candidate
        x, meta = make_tensor(c, eligible(c, [a,b]), rows=3, width=9, debug=True)
        self.assertEqual((meta["coverage"], meta["alt_count"], meta["ref_count"]), (2,1,1))
        self.assertEqual(x[5, 0, 1:8].tolist(), [BASES[v] for v in "TTACGAT"])
        self.assertEqual(x[5, 1, 1:8].tolist(), [BASES[v] for v in "GCACGAT"])
        self.assertFalse(x[:, 2].any())
        self.assertEqual(len(meta["rows"][0]["path"]), 3)
        self.assertEqual([c["node_id"] for c in meta["rows"][0]["columns"][1:8]], [2,2,1,1,1,4,4])

    def test_reverse_equivalence_and_quality_orientation(self):
        seq = {1: "AACG", 2: "TC"}
        forward = read([(2,0,False,[(2,2,"")]), (1,0,False,[(1,1,""),(1,1,"T"),(2,2,"")])], seq)
        a = alignment([(1,0,True,[(2,2,""),(1,1,"A"),(1,1,"")]), (2,0,True,[(2,2,"")])], seq)
        a.quality = bytes([10,11,12,13,14,15])
        reverse = decode_alignment(a, seq)[0]
        c = forward.observations[0].candidate
        self.assertEqual(c, reverse.observations[0].candidate)
        x, meta = make_tensor(c, eligible(c,[reverse]), width=9, debug=True)
        start = meta["candidate_columns"][0]
        self.assertEqual(x[0,0,start], BASES["T"])
        self.assertEqual(x[1,0,start], 12)
        self.assertEqual(x[5,0,start-3:start+3].tolist(), [BASES[v] for v in "TCAACG"])
        self.assertTrue(meta["rows"][0]["reversed_for_candidate"])
        for edits, rev_edits in (([(2,2,""),(0,2,"TG"),(2,2,"")],[(2,2,""),(0,2,"CA"),(2,2,"")]),
                                ([(1,1,""),(2,0,""),(1,1,"")],[(1,1,""),(2,0,""),(1,1,"")])):
            f = read([(1,0,False,edits)],seq)
            r = read([(1,0,True,rev_edits)],seq)
            self.assertEqual(f.observations[0].candidate, r.observations[0].candidate)
            xf, _ = make_tensor(f.observations[0].candidate, eligible(f.observations[0].candidate,[f]), width=9)
            xr, _ = make_tensor(r.observations[0].candidate, eligible(r.observations[0].candidate,[r]), width=9)
            np.testing.assert_array_equal(xf,xr)

    def test_shared_insertion_gap_padding_and_alleles(self):
        seq = {1: "ACGT"}
        a = read([(1,0,False,[(2,2,""),(0,2,"TA"),(2,2,"")])],seq)
        b = read([(1,0,False,[(4,4,"")])],seq)
        d = read([(1,0,False,[(2,2,""),(0,1,"G"),(2,2,"")])],seq)
        terminal = read([(1,0,False,[(2,2,"")])],seq)
        c = a.observations[0].candidate
        x, m = make_tensor(c,eligible(c,[a,b,d,terminal]),width=8,debug=True)
        lo,hi=m["candidate_columns"]
        self.assertEqual((m["coverage"],m["alt_count"],m["ref_count"],m["other_count"]),(4,1,1,2))
        self.assertEqual(hi-lo,2)
        self.assertEqual(x[0,0,lo:hi].tolist(),[BASES["T"],BASES["A"]])
        self.assertTrue(np.all(x[5,:4,lo:hi] == 6))
        self.assertTrue(np.all(x[2,:4,lo:hi] & 2))
        for row, detail in enumerate(m["rows"]):
            if detail["record_sha256"] == b.digest:
                self.assertEqual(x[0,row,lo:hi].tolist(),[6,6])
                self.assertEqual(x[1,row,lo],-1)
            if detail["record_sha256"] == terminal.digest:
                self.assertTrue(np.all(x[0,row,lo:hi] == 0))
            if detail["record_sha256"] == d.digest:
                self.assertEqual(x[0,row,lo:hi].tolist(),[BASES["G"],6])

    def test_position_coverage_deletion_spans_and_row_limits(self):
        seq={1:"ACGTAC"}
        c=Candidate(1,2,"G","T","SNP")
        before=read([(1,0,False,[(2,2,"")])],seq)
        deletion=read([(1,0,False,[(1,1,""),(4,0,""),(1,1,"")])],seq)
        alt=read([(1,0,False,[(2,2,""),(1,1,"T"),(3,3,"")])],seq)
        ref=read([(1,0,False,[(6,6,"")])],seq)
        reads=[before,deletion,alt,alt,ref]
        e=eligible(c,reads)
        self.assertEqual(len(e),4)
        _,small=make_tensor(c,e,rows=1)
        _,large=make_tensor(c,e,rows=10)
        for k in ("coverage","af","alt_count","ref_count","other_count"):
            self.assertEqual(small[k],large[k])
        self.assertEqual(small["af"],.5)
        self.assertEqual(small["selected_alignments"],1)
        self.assertEqual(large["selected_alignments"],4)
        first,_=make_tensor(c,e,rows=3)
        shuffled,_=make_tensor(c,list(reversed(e)),rows=3)
        np.testing.assert_array_equal(first,shuffled)
        ins=Candidate(1,2,"","T","INS")
        self.assertEqual(overlap(deletion,ins,10)[0],"other")
        # A partial affected interval counts, but is insufficient REF evidence.
        dc=Candidate(1,1,"CGTA","","DEL")
        self.assertEqual(overlap(before,dc,10)[0],"other")

    def test_repeated_visits_count_once_and_keep_path(self):
        seq={1:"ACGT",2:"T"}
        r=read([(1,0,False,[(4,4,"")]),(2,0,False,[(1,1,"")]),
                (1,0,False,[(1,1,""),(1,1,"T"),(2,2,"")])],seq)
        c=r.observations[0].candidate
        e=eligible(c,[r,r])
        _,m=make_tensor(c,e,debug=True)
        self.assertEqual(m["coverage"],2)
        self.assertEqual(m["alt_count"],2)
        self.assertEqual(m["rows"][0]["anchor_mapping_index"],2)
        self.assertEqual(len(m["rows"][0]["path"]),3)
        self.assertEqual([n["node_id"] for n in m["row_groups"][0]["path"]],[1,2,1])

    def test_adjacent_edits_limits_and_central_context_insertions(self):
        seq = {1: "A"*60}
        a = alignment([(1,0,False,[(1,1,""),(0,30,"T"*30),(0,21,"G"*21),
                                   (30,0,""),(21,0,"")])],seq)
        r, rejected = decode_alignment(a,seq)
        self.assertFalse(r.observations)
        self.assertEqual([e["event_length"] for e in rejected],[51,51])
        seq = {1: "ACGTAC"}
        alt = read([(1,0,False,[(1,1,""),(4,0,""),(1,1,"")])],seq)
        other = read([(1,0,False,[(2,2,""),(0,2,"TT"),(4,4,"")])],seq)
        c = alt.observations[0].candidate
        x,m = make_tensor(c,eligible(c,[alt,other]),width=12,debug=True)
        lo,hi = m["candidate_columns"]
        self.assertEqual(hi-lo,6)
        self.assertEqual(x[0,1,lo:hi].tolist(),[BASES[b] for b in "CTTGTA"])
        self.assertEqual(x[5,1,lo:hi].tolist(),[BASES[b] for b in "C--GTA"])
        self.assertEqual(m["other_count"],1)

    def test_insertion_node_edges_and_context_branch_insertions(self):
        seq={1:"AC",2:"GT",3:"TA"}
        c=Candidate(1,2,"","T","INS")
        crossing=read([(1,0,False,[(2,2,"")]),(2,0,False,[(2,2,"")])],seq)
        terminal=read([(1,0,False,[(2,2,"")])],seq)
        self.assertEqual(overlap(crossing,c,10)[0],"ref")
        self.assertEqual(overlap(terminal,c,10)[0],"other")
        partial=read([(1,0,False,[(1,1,"")]),(2,0,False,[(2,2,"")])],seq)
        self.assertEqual(overlap(partial,Candidate(1,1,"","T","INS"),10)[0],"other")
        # Start boundary counts, and a preceding mapping can supply REF evidence.
        c0=Candidate(2,0,"","T","INS")
        self.assertEqual(overlap(crossing,c0,10)[0],"ref")
        c=Candidate(1,1,"C","T","SNP")
        a=read([(1,0,False,[(1,1,""),(1,1,"T")]),
                (2,0,False,[(1,1,""),(0,1,"A"),(1,1,"")])],seq)
        b=read([(1,0,False,[(2,2,"")]),(3,0,False,[(2,2,"")])],seq)
        x,m=make_tensor(c,eligible(c,[a,b]),width=9)
        lo=m["candidate_columns"][0]
        self.assertEqual(x[5,0,lo+1:lo+4].tolist(),[BASES[v] for v in "G-T"])
        self.assertEqual(x[5,1,lo+1:lo+3].tolist(),[BASES[v] for v in "TA"])

    def test_low_quality_alt_policy_matches_raw_edit_audit(self):
        from indexed_gam_pipeline.compare_formats import audit_snp
        seq={1:"C"}
        a=alignment([(1,0,True,[(1,1,"T")])],seq)
        a.quality=bytes([3])
        r,_=decode_alignment(a,seq)
        c=r.observations[0].candidate
        self.assertEqual(overlap(r,c,10)[0],"other")
        self.assertEqual(overlap(r,c,0)[0],"alt")
        hit=audit_snp(a,c.metadata(),1,10)
        self.assertEqual((hit["legacy"],hit["candidate_v2"],hit["base_quality"]),("alt","other",3))
        _,m=make_tensor(c,eligible(c,[r]))
        self.assertEqual(m["coverage"],1)
        self.assertEqual(m["other_count"],1)

    def test_node_path_grouping_preserves_selection_counts_and_channels(self):
        seq={1:"AC",2:"GG",3:"TT"}
        rows=[]
        for branch,base in ((3,"T"),(2,"C"),(2,"T"),(3,"C")):
            edits=[(1,1,""),(1,1,"T")] if base=="T" else [(2,2,"")]
            rows.append(read([(branch,0,False,[(2,2,"")]),(1,0,False,edits)],seq))
        c=Candidate(1,1,"C","T","SNP")
        e=eligible(c,rows)
        x,m=make_tensor(c,e,rows=6,width=9,debug=True)
        self.assertEqual([r["support"] for r in m["rows"]],["alt","ref","alt","ref"])
        self.assertEqual([(g["start_row"],g["end_row"]) for g in m["row_groups"]],[(0,2),(2,4)])
        self.assertEqual(set(g["path"][0]["node_id"] for g in m["row_groups"]),{2,3})
        self.assertEqual((m["coverage"],m["alt_count"],m["ref_count"],m["af"]),(4,2,2,.5))
        self.assertFalse(x[:,4:].any())
        # All six channels and debug column maps travel with their record.
        for ri,detail in enumerate(m["rows"]):
            original=next(item for item in e if item[0].digest==detail["record_sha256"])
            single,meta=make_tensor(c,[original],rows=1,width=9,debug=True)
            np.testing.assert_array_equal(x[:,ri],single[:,0])
            self.assertEqual(detail["columns"],meta["rows"][0]["columns"])
        capped,small=make_tensor(c,e,rows=2,width=9,debug=True)
        self.assertEqual(small["selected_counts"],{"alt":1,"ref":1})
        self.assertEqual(small["selected_grouped_ranks"],[0,3])
        self.assertEqual(small["coverage"],4)
        self.assertEqual(small["af"],.5)
        shuffled,_=make_tensor(c,list(reversed(e)),rows=6,width=9)
        np.testing.assert_array_equal(x,shuffled)

    def test_node_group_key_ignores_distant_branches_and_normalizes_strand(self):
        seq={1:"AC",2:"G"*20,3:"TT",4:"AA"}
        a=read([(3,0,False,[(2,2,"")]),(2,0,False,[(20,20,"")]),
                (1,0,False,[(1,1,""),(1,1,"T")])],seq)
        b=read([(1,0,True,[(1,1,"A"),(1,1,"")]),(2,0,True,[(20,20,"")]),
                (4,0,True,[(2,2,"")])],seq)
        c=a.observations[0].candidate
        _,m=make_tensor(c,eligible(c,[a,b]),width=9,debug=True)
        self.assertEqual(len(m["row_groups"]),1)
        self.assertEqual(m["row_groups"][0]["path"],[{"node_id":2,"reverse":False},{"node_id":1,"reverse":False}])

    def test_synthetic_shards_summary_manifest(self):
        with tempfile.TemporaryDirectory() as directory:
            path,_=fixture(directory)
            folder=Path(directory)
            (folder/"nodes.txt").write_text("10\n20\n30\n1000\n")
            (folder/"nodes.json").write_text(json.dumps([dict(node_id=n,sequence="AAAAAA") for n in (10,20,30,1000)]))
            args=argparse.Namespace(command="build",format="candidate-v2",gam=str(path),index=None,
                nodes=str(folder/"nodes.txt"),node_json=str(folder/"nodes.json"),node_sqlite=None,gfa=None,
                output=str(folder/"out"),batch_nodes=2,max_node_span=100,max_batch_segments=100,
                shard_size=2,min_mapq=10,min_af=.05,min_variants=1,min_allele_bq=10.,
                variant_type="all",max_indel_len=50,rows=2,width=100,debug_rows=True)
            with redirect_stdout(io.StringIO()):
                build(args)
            manifest=json.loads((folder/"out/manifest.json").read_text())
            summary=[json.loads(s) for s in (folder/"out/variant_summary.ndjson").read_text().splitlines()]
            self.assertEqual(manifest["tensors"],len(summary))
            self.assertEqual(manifest["shape"],[6,2,100])
            self.assertEqual(len(summary),4)
            for i,meta in enumerate(summary):
                self.assertEqual((meta["shard_index"],meta["index_within_shard"]),(i//2,i%2))
                x=np.load(folder/f"out/shard_{i//2:05d}_data.npy")[i%2]
                self.assertEqual(list(x.shape),manifest["shape"])
                self.assertEqual(str(x.dtype),manifest["dtype"])
                self.assertEqual(meta["selected_alignments"],min(2,meta["coverage"]))
                self.assertEqual(meta["coverage"],meta["alt_count"]+meta["ref_count"]+meta["other_count"])
            # The standalone real-example auditor must also reject corrupted
            # reference channels, not merely accept self-consistent metadata.
            import sqlite3
            from indexed_gam_pipeline.validate_examples import validate
            db=folder/"graph.sqlite"
            with sqlite3.connect(db) as connection:
                connection.execute("CREATE TABLE nodes (node_id TEXT PRIMARY KEY, seq TEXT)")
                connection.executemany("INSERT INTO nodes VALUES (?, ?)",[(str(n),"AAAAAA") for n in (10,20,30,1000)])
            audit=validate(folder/"out",str(path),None,str(db))
            self.assertTrue(audit["passed"])
            self.assertEqual(audit["examples"],4)
            shard=folder/"out/shard_00000_data.npy"
            corrupted=np.load(shard)
            corrupted[0,5,0,49]=BASES["C"]
            np.save(shard,corrupted)
            with self.assertRaisesRegex(ValueError,"graph base"):
                validate(folder/"out",str(path),None,str(db))
