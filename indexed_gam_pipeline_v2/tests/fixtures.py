"""Synthetic GAM / graph-index fixtures shared by the v2 tests. Import this module first."""
import argparse
import json
from pathlib import Path
import sqlite3
import sys

import pysam

ROOT = Path(__file__).resolve().parents[2]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from indexed_gam_pipeline_v2 import vg_pb2  # noqa: E402
from indexed_gam_pipeline_v2.candidates import rc  # noqa: E402
from indexed_gam_pipeline_v2.gam_reader import build_index, encode_varint  # noqa: E402
from indexed_gam_pipeline_v2.graph_index import METRIC, SCHEMA  # noqa: E402


def write_gam(path, alignments, foreign_groups=()):
    """Write one alignment per BGZF block (flush per record) and build its .gai.

    `foreign_groups` lists record indices before which a non-GAM tagged group (like
    giraffe's PARAMS_JSON) is inserted.
    """
    path = Path(path)
    with pysam.BGZFile(str(path), "wb") as stream:
        for i, a in enumerate(alignments):
            if i in foreign_groups:
                payload = b'{"parameters": true}'
                stream.write(encode_varint(2) + encode_varint(11) + b"PARAMS_JSON" + encode_varint(len(payload)) + payload)
                stream.flush()
            raw = a.SerializeToString()
            stream.write(encode_varint(2) + encode_varint(3) + b"GAM" + encode_varint(len(raw)) + raw)
            stream.flush()
    build_index(path, str(path) + ".gai")
    return path


def simple_alignment(nodes=(10,), reverse=False, mapq=60, mutation=True, name="same-name"):
    """Six bases per node: AATAAA with a T>A-site substitution when `mutation`, else AAAAAA."""
    a = vg_pb2.Alignment(name=name, mapping_quality=mapq)
    for nid in nodes:
        m = a.path.mapping.add()
        m.position.node_id = nid
        m.position.is_reverse = reverse
        if mutation:
            m.edit.add(from_length=2, to_length=2)
            m.edit.add(from_length=1, to_length=1, sequence="T")
            m.edit.add(from_length=3, to_length=3)
            a.sequence += "AATAAA"
        else:
            m.edit.add(from_length=6, to_length=6)
            a.sequence += "AAAAAA"
    a.quality = bytes([30] * len(a.sequence))
    return a


def spec_alignment(specs, sequences, name="same-name", mapq=60):
    """specs = [(node, offset, reverse, [(from_length, to_length, replacement), ...]), ...]."""
    a = vg_pb2.Alignment(name=name, mapping_quality=mapq)
    for nid, offset, reverse, edits in specs:
        m = a.path.mapping.add()
        m.position.node_id, m.position.offset, m.position.is_reverse = nid, offset, reverse
        ref = rc(sequences[nid]) if reverse else sequences[nid]
        cursor = offset
        for f, t, seq in edits:
            m.edit.add(from_length=f, to_length=t, sequence=seq)
            a.sequence += seq if seq else ref[cursor:cursor + f] if t else ""
            cursor += f
    a.quality = bytes([30] * len(a.sequence))
    return a


def tiny_gam(directory):
    """Six records over nodes 10/20/30/1000: shared names, exact duplicates, multi-node reads."""
    rows = [simple_alignment((10,)), simple_alignment((10,)), simple_alignment((10, 30)),
            simple_alignment((20,)), simple_alignment((20, 1000)), simple_alignment((30,))]
    return write_gam(Path(directory) / "tiny.gam", rows), rows


def af_gam(directory):
    """Five nodes x 50 reads with 3 or 4 ALT records each: SNP on 10, DEL on 20/30, INS on 40/50."""
    rows = []
    for node, alt_count in [(10, 3), (20, 3), (30, 4), (40, 3), (50, 4)]:
        for i in range(50):
            a = simple_alignment((node,), mutation=(node == 10 and i < alt_count), name=f"{node}-{i}")
            if node != 10 and i < alt_count:
                m = a.path.mapping[0]
                del m.edit[:]
                m.edit.add(from_length=2, to_length=2)
                if node < 40:
                    m.edit.add(from_length=1, to_length=0)
                    m.edit.add(from_length=3, to_length=3)
                    a.sequence = "AAAAA"
                else:
                    m.edit.add(from_length=0, to_length=1, sequence="T")
                    m.edit.add(from_length=4, to_length=4)
                    a.sequence = "AATAAAA"
                a.quality = bytes([30] * len(a.sequence))
            rows.append(a)
    return write_gam(Path(directory) / "af.gam", rows), rows


def graph_fixture(path, nodes):
    """A complete unified graph index from [(node_id, sequence, distinct_path_count), ...]."""
    with sqlite3.connect(path) as db:
        db.execute("CREATE TABLE graph_metadata(value TEXT)")
        db.execute("INSERT INTO graph_metadata VALUES(?)",
                   (json.dumps(dict(schema=SCHEMA, metric=METRIC, status="complete", nodes=len(nodes))),))
        db.execute("CREATE TABLE nodes(node_id INTEGER PRIMARY KEY, seq TEXT, distinct_path_count INTEGER)")
        db.executemany("INSERT INTO nodes VALUES(?,?,?)", nodes)
    return path


def build_args(**overrides):
    """argparse.Namespace with the run.py build defaults, tuned small for synthetic data."""
    args = dict(command="build", gam=None, index=None, nodes=None, graph_index=None, output=None, chromosomes="all", chr_index=None,
                snv_output=None, indel_output=None, snv_min_af=None, indel_min_af=None,
                rows=200, width=101, debug_rows=False, gam_cache_mb=1, batch_nodes=2, max_node_span=100,
                max_batch_alignments=1000, shard_size=2, max_tensors=None, min_mapq=10, min_af=0.05,
                min_variants=1, min_allele_bq=10.0, max_indel_len=50, variant_type="all",
                candidate_unit="site", max_node_reads=800, early_af_filter=True)
    args.update(overrides)
    return argparse.Namespace(**args)
