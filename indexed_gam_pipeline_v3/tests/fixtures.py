"""Synthetic data and helpers shared by the tests and the golden cases (tests/golden.py).

    GAM/GAI writers   write_gam, build_index (GAI v1), encode_varint
    alignments        simple_alignment, spec_alignment, mirror
    datasets          tiny_gam, af_gam, graph_fixture; SEQ, site_rows, mixed_af_rows
    arguments         build_args (a complete builder Namespace), split_args (shared/SNV/INDEL outputs)
    test-only         overlap, make_tensor (single-candidate wrappers around NodeReads / make_site_tensor)
    labels            CHR1, PIECES, CHR2, CHR3, graph_files (GFA + FASTA), write_vcf
"""
import argparse
import gzip
import json
from pathlib import Path
import random
import sqlite3

import pysam

from .. import vg_pb2
from ..candidates import NodeReads, make_site_tensor, rc
from ..gam_reader import decode, group
from ..graph_index import METRIC, SCHEMA


# --- GAM/GAI writers (moved verbatim from v2 gam_reader: only the fixtures write GAIs) ---

def encode_varint(value):
    out = bytearray()
    while value > 127:
        out.append((value & 127) | 128)
        value >>= 7
    out.append(value)
    return bytes(out)


def build_index(gam, output):
    """Create a GAI v1 from complete BGZF GAM groups and their node ID spans."""
    bins = {}
    groups = alignments = 0
    with pysam.BGZFile(str(gam), "rb") as stream:
        while True:
            start = stream.tell()
            messages = group(stream)
            if messages is None:
                break
            end = stream.tell()
            node_ids = [m.position.node_id for raw in messages
                        for m in decode(raw).path.mapping if m.position.node_id > 0]
            alignments += len(messages)
            groups += 1
            if not node_ids:
                continue
            lo, hi = min(node_ids), max(node_ids)
            shift = max(1, (lo ^ hi).bit_length())
            if shift > 64:
                raise ValueError("Node ID exceeds uint64")
            number = (lo >> shift) + (1 << (64 - shift)) - 1
            bins.setdefault(number, []).append((start, end))
    with gzip.open(output, "wb") as stream:
        stream.write(b"GAI!" + encode_varint(1) + encode_varint(len(bins)))
        for number, runs in sorted(bins.items()):
            stream.write(encode_varint(number) + encode_varint(len(runs)))
            for start, end in runs:
                stream.write(encode_varint(start) + encode_varint(end))
        stream.write(encode_varint(0))
    return {"groups": groups, "alignments": alignments, "bins": len(bins)}


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
    """A complete argparse.Namespace with the run.py build defaults, tuned small for synthetic data."""
    args = dict(command="build", gam=None, index=None, nodes=None, graph_index=None, output=None, chromosomes="all", chr_index=None,
                snv_output=None, indel_output=None, snv_min_af=None, indel_min_af=None,
                rows=200, width=101, debug_rows=False, gam_cache_mb=1, batch_nodes=2, max_node_span=100,
                max_batch_alignments=1000, shard_size=2, min_mapq=10, min_af=0.05,
                min_variants=1, min_allele_bq=10.0, max_indel_len=50,
                max_node_reads=800, early_af_filter=True, decoder="auto")
    args.update(overrides)
    return argparse.Namespace(**args)


def split_args(root, name, **overrides):
    """build_args with the shared/SNV/INDEL outputs under root/name and both AF thresholds 0.05."""
    folder = Path(root) / name
    args = dict(output=str(folder / "shared"), snv_output=str(folder / "SNV"), indel_output=str(folder / "INDEL"),
                snv_min_af=0.05, indel_min_af=0.05)
    args.update(overrides)
    return build_args(**args)


# --- test-only single-candidate wrappers (moved verbatim from v2 candidates) ---------

def overlap(read, candidate, min_bq):
    """Classify one record as ('alt'|'ref'|'other', anchor visit), or None if it does not cover.

    See NodeReads.classify for the rule; this single-record form returns the Visit itself.
    """
    hit = NodeReads(candidate.node, [read], len(candidate.ref) + 2).classify(read, candidate, min_bq)
    return None if hit is None else (hit[0], hit[1].visit)


def make_tensor(candidate, eligible, path_counts, rows=200, width=101, debug=False):
    """One-allele site tensor: `eligible` = [(read, support, anchor Visit or VisitView)]."""
    return make_site_tensor([candidate], [eligible], path_counts, rows, width, debug)


def mirror(specs, sequences):
    """The same molecule as the other strand would report it: reversed path, flipped mappings, mirrored edits."""
    mirrored = []
    for node, offset, reverse, edits in reversed(specs):
        span = sum(f for f, _, _ in edits)
        mirrored.append((node, len(sequences[node]) - offset - span, not reverse,
                         [(f, t, rc(seq) if seq else seq) for f, t, seq in reversed(edits)]))
    return mirrored


# --- site fixtures (multi-allelic sites, AF prefilter) ----------------------------

SEQ = {10: "ACGTACGTAC", 20: "TTGCAAGGCT"}


def site_rows():
    """Node 10: position 4 A>C x4, A>G x3, A>T x1; 1-bp DEL at 7 x3 and INS TT at 7 x2 (11 reads;
    neither indel can left-shift: the base before position 7 is G). Node 20: position 3 C>A x3 (5 reads)."""
    snv = ["C"] * 4 + ["G"] * 3 + ["T"] + [None] * 3
    indel = {0: "D", 4: "D", 8: "D", 1: "I", 9: "I"}
    rows = []
    for i, base in enumerate(snv):
        edits = [(4, 4, ""), (1, 1, base or "")]
        edits += {"D": [(2, 2, ""), (1, 0, ""), (2, 2, "")],
                  "I": [(2, 2, ""), (0, 2, "TT"), (3, 3, "")]}.get(indel.get(i), [(5, 5, "")])
        rows.append(spec_alignment([(10, 0, False, edits)], SEQ, name=f"n10_{i}"))
    for i in range(5):
        edits = [(3, 3, ""), (1, 1, "A"), (6, 6, "")] if i < 3 else [(10, 10, "")]
        rows.append(spec_alignment([(20, 0, False, edits)], SEQ, name=f"n20_{i}"))
    return rows


def mixed_af_rows():
    """Node 10 (8 reads): A>C at 4 in 5, G>T at 2 in 1, 1-bp DEL at 6 in 1, INS GG at 8 in 1.
    Node 20 (5 reads): C>A at 3 in 2, G>C at 7 in 1, plus one all-match reverse-strand record."""
    rows = []
    for i in range(8):
        e = [(4, 4, ""), (1, 1, "C" if i < 4 else ""), (5, 5, "")]
        if i == 4:
            e = [(2, 2, ""), (1, 1, "T"), (7, 7, "")]
        if i == 5:
            e = [(4, 4, ""), (1, 1, "C"), (1, 1, ""), (1, 0, ""), (3, 3, "")]
        if i == 6:
            e = [(8, 8, ""), (0, 2, "GG"), (2, 2, "")]
        rows.append(spec_alignment([(10, 0, False, e)], SEQ, name=f"a{i}"))
    for i in range(5):
        e = [(3, 3, ""), (1, 1, "A"), (6, 6, "")] if i < 2 else [(10, 10, "")]
        if i == 2:
            e = [(7, 7, ""), (1, 1, "C"), (2, 2, "")]
        rows.append(spec_alignment([(20, 0, i == 4, e)], SEQ, name=f"b{i}"))
    return rows


# --- label fixtures (reference path, chromosome blocks, truth VCFs) -----------------

rng = random.Random(5)
# chr1: 80 bp with a T run at 3..7 crossing the node1 | node2 junction (node 2 is reverse-oriented).
CHR1 = "ACGTTTTTGCAACACACGTAGGCT" + "".join(rng.choice("ACGT") for _ in range(56))
PIECES = [(1, 5, False), (2, 7, True), (3, 12, False), (4, 20, True), (5, 16, False), (6, 20, False)]
CHR2 = "GATTACAGATTACA"
# chr3: nodes 10 and 11 both reverse-oriented on the reference walk (<10<11).
CHR3 = "".join(rng.choice("ACGT") for _ in range(30))


def graph_files(directory):
    """GFA (chr1 walk over nodes 1-6, chr2 walk 8>9>8, HG1 walk through off-reference node 7) and FASTA."""
    directory = Path(directory)
    segments, walk, cursor = [], "", 0
    for node, length, reverse in PIECES:
        piece = CHR1[cursor:cursor + length]
        segments.append((node, rc(piece) if reverse else piece))
        walk += ("<" if reverse else ">") + str(node)
        cursor += length
    assert cursor == len(CHR1)
    segments += [(7, "GG"), (8, CHR2[:7]), (9, CHR2[7:]), (10, rc(CHR3[:12])), (11, rc(CHR3[12:]))]
    lines = ["H\tVN:Z:1.1"] + [f"S\t{n}\t{s}" for n, s in segments]
    lines += [f"W\tGRCh38\t0\tchr1\t0\t{len(CHR1)}\t{walk}",
              f"W\tGRCh38\t0\tchr2\t0\t{len(CHR2) + 7}\t>8>9>8",
              "W\tHG1\t1\tctg1\t0\t9\t>1>7>3",
              f"W\tGRCh38\t0\tchr3\t0\t{len(CHR3)}\t<10<11"]
    (directory / "g.gfa").write_text("\n".join(lines) + "\n")
    (directory / "g.fa").write_text(f">chr1\n{CHR1}\n>chr2\n{CHR2}{CHR2[:7]}\n>chr3\n{CHR3}\n")
    pysam.faidx(str(directory / "g.fa"))
    return directory / "g.gfa", directory / "g.fa"


def write_vcf(path, records):
    header = ["##fileformat=VCFv4.2", "##contig=<ID=chr1,length=80>", "##contig=<ID=chr2,length=21>",
              '##FILTER=<ID=GAP1,Description="gap">', '##FORMAT=<ID=GT,Number=1,Type=String,Description="GT">',
              "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tS"]
    body = [f"chr1\t{pos}\t.\t{ref}\t{alt}\t.\t{filt}\t.\tGT\t{gt}" for pos, ref, alt, filt, gt in records]
    Path(path).write_text("\n".join(header + body) + "\n")
    return path
