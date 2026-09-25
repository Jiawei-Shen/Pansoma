"""Golden outputs: fingerprints of production-shaped builds and one orchestrated run, recorded from v2.

    python -m indexed_gam_pipeline_v3.tests.golden record --package indexed_gam_pipeline_v2 [--output FILE]
    python -m indexed_gam_pipeline_v3.tests.golden check --package indexed_gam_pipeline_v3
        [--decoder python|native] [--case G1 ...] [--keep DIR]

Cases (all split SNV/INDEL builds through the package's CLI):
    G1   tiny_gam: default (int) --min-allele-bq, --min-variants 1, batch 2, shard 2, cache 1 MiB
    G2   af_gam: --min-allele-bq 10, --debug-rows, shard 1, AF 0.06/0.08
    G3   site_rows: multi-allelic sites, --debug-rows, rows 7, AF 0.2/0.1, --min-variants 2
    G4a  mixed_af_rows: --no-early-af-filter, --max-node-reads 0, AF 0.3/0.3, rows 4, --debug-rows
    G4b  G4a with --early-af-filter
    G5   12 records on one node, --max-node-reads 5, --debug-rows, rows 4
    G6   displaced repeat {1:C, 2:AT, 3:AT, 4:G}: target [3], then the supplement build of [2]
    G7   multi-node deletion {1:GCA, 2:TG, 3:GAT}, --debug-rows
    G8   random_chain_world(20260924): ~40 nodes, ~400 records, --chromosomes autosome (last 10 nodes off)
    O1   mini_world: orchestrate prepare (3 tasks, 2 processes, node stats, autosome selection, supplement
         rounds, merge 4, --keep-sources, reference path, labels), then run.sh (SLURM_CPUS_PER_TASK=2)

Every build is a subprocess (cwd = repository root; O1's run from <root>/source via run.sh), so a
package's native module is the only one loaded in its process and no pipeline package other than
this one is imported here. Inputs are generated here from tests/fixtures.py (identical for every
package; their hashes are recorded) and outputs are compared as tools.compare_runs fingerprints with
the case directory as the <ROOT> anchor. The recorder requires v2 python == v2 native, the expected
decoder in every manifest, G8 independent of PYTHONHASHSEED and O1's coverage (a supplement node, a
removed target, merged tensors in two autosome blocks, somatic/germline/non/ignore labels).
"""
import argparse
from concurrent.futures import ThreadPoolExecutor
import csv
import gzip
import hashlib
import importlib
import json
import os
from pathlib import Path
import random
import shutil
import subprocess
import sys
import tempfile

import pysam

from ..tensor_postprocessing.chr_index import FIELDS
from ..tools.compare_runs import differences, fingerprint
from ..tools.graph_prep import scan
from .fixtures import (SEQ, af_gam, graph_fixture, mirror, mixed_af_rows, site_rows, spec_alignment, tiny_gam,
                       write_gam)

REPO = Path(__file__).resolve().parents[2]
PACKAGE = __package__.split(".")[0]
HASHES = Path(__file__).with_name("golden_hashes.json")
FIXED_MTIME_NS = 1_700_000_000 * 10 ** 9  # inputs get one mtime, so their stamps are reproducible
DEFAULT_AF = ["--snv-min-af", "0.05", "--indel-min-af", "0.05"]
DECODERS = ("python", "native")
HASH_SEEDS = (1, 2)
CONSTANTS_SNIPPET = r"""
import importlib, json, sys
package = sys.argv[1]
names = {
    "candidates": ("FORMAT_VERSION", "SCHEMA_VERSION", "STORAGE_VERSION", "ROW_SELECTION_VERSION",
                   "WINDOW_ENCODING_VERSION", "ROW_ORDER", "CHANNELS", "BASES", "OPS", "STRAND", "COUNT_LINEAR_MAX",
                   "REF_LABEL", "OTHER_LABEL"),
    "build": ("KIND", "SITE_UNIT", "PARAMETERS", "ALLELE_FIELDS", "SITE_DEFINITION", "READ_CAP_RULE"),
    "orchestrate": ("BUILDER_OPTIONS", "LABEL_INPUTS"),
    "tensor_postprocessing.merge_shards": ("LAYOUT", "DEFAULT_SHARD_SIZE", "POSITION", "ADDED", "SHARED_KEYS",
                                           "COPIED_KEYS", "AUDIT_STREAMS"),
    "tensor_postprocessing.truth_labels": ("VERSION", "LABELS", "NEAR_BP"),
}
result = {}
for module, keys in names.items():
    loaded = importlib.import_module(f"{package}.{module}")
    result.update({f"{module}.{key}": getattr(loaded, key) for key in keys})
print(json.dumps(result))
"""


# --- inputs -----------------------------------------------------------------------

def write_nodes(path, nodes):
    Path(path).write_text("".join(f"{n}\n" for n in nodes))
    return Path(path)


def chr_table(path, blocks):
    """A node ID -> chromosome block table from [(chrom, first node, last node, dataset), ...]."""
    with Path(path).open("w", newline="") as stream:
        writer = csv.DictWriter(stream, FIELDS, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows([dict(chrom=c, first_node=a, last_node=b, nodes=b - a + 1, dataset=d, source="golden")
                          for c, a, b, d in blocks])
    return Path(path)


def sources(gam, graph, nodes):
    return ["--gam", str(gam), "--graph-index", str(graph), "--nodes", str(nodes)]


def input_hashes(directory):
    """{relative path: sha256} of every input; a .gai by content (gzip stores its write time)."""
    result = {}
    for path in sorted(Path(directory).rglob("*")):
        if path.is_file():
            data = path.read_bytes()
            result[path.relative_to(directory).as_posix()] = hashlib.sha256(
                gzip.decompress(data) if path.suffix == ".gai" else data).hexdigest()
    return result


def freeze_times(directory):
    for path in Path(directory).rglob("*"):
        if path.is_file():
            os.utime(path, ns=(FIXED_MTIME_NS, FIXED_MTIME_NS))


# --- build cases: world(inputs) -> [(output name, builder arguments)] --------------------

def g1(inputs):
    gam, _ = tiny_gam(inputs)
    graph = graph_fixture(inputs / "graph.sqlite", [(n, "AAAAAA", 86) for n in (10, 20, 30, 1000)])
    nodes = write_nodes(inputs / "nodes.txt", (10, 20, 30, 1000))
    return [("out", sources(gam, graph, nodes) + ["--min-variants", "1", "--batch-nodes", "2", "--shard-size", "2",
                                                   "--gam-cache-mb", "1"])]


def g2(inputs):
    gam, _ = af_gam(inputs)
    graph = graph_fixture(inputs / "graph.sqlite", [(n, "AAAAAA", 331) for n in (10, 20, 30, 40, 50)])
    nodes = write_nodes(inputs / "nodes.txt", (10, 20, 30, 40, 50))
    return [("out", sources(gam, graph, nodes) + ["--min-allele-bq", "10", "--debug-rows", "--shard-size", "1",
                                                   "--snv-min-af", "0.06", "--indel-min-af", "0.08"])]


def site_world(inputs, make_rows):
    gam = write_gam(inputs / f"{make_rows.__name__}.gam", make_rows())
    graph = graph_fixture(inputs / "graph.sqlite", [(n, s, 7) for n, s in SEQ.items()])
    return sources(gam, graph, write_nodes(inputs / "nodes.txt", sorted(SEQ)))


def g3(inputs):
    return [("out", site_world(inputs, site_rows) + ["--snv-min-af", "0.2", "--indel-min-af", "0.1",
                                                      "--min-variants", "2", "--debug-rows", "--rows", "7"])]


def g4(inputs, early):
    return [("out", site_world(inputs, mixed_af_rows) + [
        "--early-af-filter" if early else "--no-early-af-filter", "--max-node-reads", "0", "--snv-min-af", "0.3",
        "--indel-min-af", "0.3", "--rows", "4", "--debug-rows", "--min-variants", "1", "--batch-nodes", "2",
        "--shard-size", "2"])]


def g5(inputs):
    seq = {1: "ACGTACGTAC"}
    rows = [spec_alignment([(1, 0, False, [(4, 4, ""), (1, 1, "C" if i % 2 else ""), (5, 5, "")])], seq,
                           name=f"deep{i}") for i in range(12)]
    gam = write_gam(inputs / "deep.gam", rows)
    graph = graph_fixture(inputs / "graph.sqlite", [(1, seq[1], 3)])
    return [("out", sources(gam, graph, write_nodes(inputs / "nodes.txt", [1])) + [
        "--max-node-reads", "5", "--debug-rows", "--rows", "4", "--min-variants", "1"])]


def displaced_repeat_rows():
    """{1:C, 2:AT, 3:AT, 4:G}: forward records insert AT at the repeat's right end (node 3), reverse
    records at its left end; left-normalization moves every copy to node 2 (6 ALT, 4 REF records)."""
    sequences = {1: "C", 2: "AT", 3: "AT", 4: "G"}
    right = [(1, 0, False, [(1, 1, "")]), (2, 0, False, [(2, 2, "")]), (3, 0, False, [(2, 2, ""), (0, 2, "AT")]),
             (4, 0, False, [(1, 1, "")])]
    left = [(1, 0, False, [(1, 1, "")]), (2, 0, False, [(0, 2, "AT"), (2, 2, "")]), (3, 0, False, [(2, 2, "")]),
            (4, 0, False, [(1, 1, "")])]
    ref = [(n, 0, False, [(len(s), len(s), "")]) for n, s in sequences.items()]
    rows = ([spec_alignment(right, sequences, name=f"f{i}") for i in range(3)]
            + [spec_alignment(mirror(left, sequences), sequences, name=f"r{i}") for i in range(3)]
            + [spec_alignment(ref, sequences, name=f"ref{i}") for i in range(4)])
    return sequences, rows


def g6(inputs):
    sequences, rows = displaced_repeat_rows()
    gam = write_gam(inputs / "g.gam", rows)
    graph = graph_fixture(inputs / "graph.sqlite", [(n, s, 5) for n, s in sequences.items()])
    options = ["--min-variants", "1", "--rows", "12", "--width", "11", "--max-indel-len", "5", "--batch-nodes", "4"]
    return [("main", sources(gam, graph, write_nodes(inputs / "main.txt", [3])) + options),
            ("supplement", sources(gam, graph, write_nodes(inputs / "supplement.txt", [2])) + options)]


def g7(inputs):
    sequences = {1: "GCA", 2: "TG", 3: "GAT"}
    spec = [(1, 0, False, [(2, 2, ""), (1, 0, "")]), (2, 0, False, [(2, 0, "")]),
            (3, 0, False, [(1, 0, ""), (2, 2, "")])]
    ref = [(n, 0, False, [(len(s), len(s), "")]) for n, s in sequences.items()]
    rows = ([spec_alignment(spec, sequences, name=f"a{i}") for i in range(3)]
            + [spec_alignment(mirror(spec, sequences), sequences, name=f"b{i}") for i in range(2)]
            + [spec_alignment(ref, sequences, name=f"r{i}") for i in range(4)])
    gam = write_gam(inputs / "g.gam", rows)
    graph = graph_fixture(inputs / "graph.sqlite", [(n, s, 5) for n, s in sequences.items()])
    return [("out", sources(gam, graph, write_nodes(inputs / "nodes.txt", [1, 2, 3])) + [
        "--min-variants", "1", "--width", "11", "--max-indel-len", "5", "--rows", "12", "--debug-rows"])]


def random_chain_world(seed, nodes=40, records=400):
    """A chain of `nodes` short nodes (random and low-complexity sequences, a rare N) and `records`
    records over windows of it, forward or reverse. Planted events recur at several allele fractions:
    SNVs, short indels (also across node junctions and in repeats), indels longer than
    --max-indel-len 10 and complex replacements; plus random errors and N read bases. MAPQ 0/5/20/60;
    qualities constant, random low or missing. Returns (sequences, records)."""
    rng = random.Random(seed)
    sequences = {}
    for node in range(1, nodes + 1):
        length = rng.randint(3, 24)
        if rng.random() < 0.3:
            seq = (rng.choice(["A", "T", "AT", "CA", "CAG", "GGC"]) * length)[:length]
        else:
            seq = "".join(rng.choice("ACGT") for _ in range(length))
        if rng.random() < 0.05:
            i = rng.randrange(length)
            seq = seq[:i] + "N" + seq[i + 1:]
        sequences[node] = seq
    starts, text = {}, ""
    for node in sorted(sequences):
        starts[node], text = len(text), text + sequences[node]
    node_of = [node for node in sorted(sequences) for _ in sequences[node]]

    def node_end(position):
        return starts[node_of[position]] + len(sequences[node_of[position]])

    def random_bases(k):
        return "".join(rng.choice("ACGT") for _ in range(k))

    events = []
    for _ in range(3 * nodes):
        p = rng.randrange(1, len(text) - 1)
        kind = rng.choices(["SNV", "INS", "DEL", "LONG_INS", "LONG_DEL", "COMPLEX"], [40, 20, 20, 4, 4, 12])[0]
        if kind == "SNV":
            ref, alt = 1, rng.choice([b for b in "ACGT" if b != text[p]])
        elif kind == "INS":
            ref, alt = 0, (text[p:p + rng.randint(1, 4)] if rng.random() < 0.5 else random_bases(rng.randint(1, 4)))
        elif kind == "DEL":
            ref, alt = rng.randint(1, 4), ""
        elif kind == "LONG_INS":
            ref, alt = 0, random_bases(rng.randint(11, 15))
        elif kind == "LONG_DEL":
            ref, alt = rng.randint(11, 14), ""
        else:
            ref = rng.randint(1, 3)
            alt = random_bases(rng.choice([k for k in (1, 2, 3, 4) if k != ref]))
            if p + ref > node_end(p):
                continue  # a complex replacement is one edit of one mapping
        if p + ref < len(text) - 1:
            events.append((p, ref, alt, rng.choice([0.1, 0.2, 0.35, 0.5, 0.8])))
    planted, end = [], 0
    for event in sorted(events):
        if event[0] > end:  # at least one reference base between events
            planted.append(event)
            end = event[0] + event[1]

    rows = []
    for index in range(records):
        s = rng.randrange(0, len(text) - 20)
        e = min(len(text), s + rng.randint(30, 120))
        edits = {}  # node -> [first chain position, edits]

        def add(position, f, t, seq):
            edits.setdefault(node_of[position], [position, []])[1].append((f, t, seq))

        def matches(a, b):
            """Reference bases [a, b), one match edit per node run, with random errors and N read bases."""
            while a < b:
                y, run = min(b, node_end(a)), a
                for position in range(a, y):
                    if text[position] != "N" and rng.random() < 0.012:
                        if position > run:
                            add(run, position - run, position - run, "")
                        add(position, 1, 1, "N" if rng.random() < 0.4 else
                            rng.choice([c for c in "ACGT" if c != text[position]]))
                        run = position + 1
                if y > run:
                    add(run, y - run, y - run, "")
                a = y

        x = s
        for p, ref, alt, af in planted:
            if not (s < p and p + ref < e) or rng.random() >= af:
                continue
            matches(x, p)
            if not ref:
                add(p, 0, len(alt), alt)
            elif alt:
                add(p, ref, len(alt), alt)
            else:  # a deletion is split at node junctions (one deletion edit per mapping)
                cut = p
                while cut < p + ref:
                    y = min(p + ref, node_end(cut))
                    add(cut, y - cut, 0, "")
                    cut = y
            x = p + ref
        matches(x, e)
        specs = [(node, first - starts[node], False, node_edits) for node, (first, node_edits) in sorted(edits.items())]
        if rng.random() < 0.5:
            specs = mirror(specs, sequences)
        a = spec_alignment(specs, sequences, name=f"r{index:03d}", mapq=rng.choice([0, 5, 20, 60, 60, 60]))
        mode = rng.random()
        if mode < 0.1:
            a.quality = b""
        elif mode < 0.4:
            a.quality = bytes(rng.randint(0, 45) for _ in a.sequence)
        rows.append(a)
    return sequences, rows


def g8(inputs):
    sequences, rows = random_chain_world(20260924)
    rng = random.Random(7)
    gam = write_gam(inputs / "chain.gam", rows)
    graph = graph_fixture(inputs / "graph.sqlite", [(n, s, rng.choice([1, 2, 5, 50, 99, 100, 101, 150, 1000]))
                                                    for n, s in sequences.items()])
    last = max(sequences)
    table = chr_table(inputs / "chr.tsv", [("chr1", 1, last - 10, "autosome"),
                                           ("chrX", last - 9, last, "non_autosomal")])
    return [("out", sources(gam, graph, write_nodes(inputs / "nodes.txt", sorted(sequences))) + [
        "--rows", "20", "--width", "21", "--max-node-reads", "30", "--batch-nodes", "7", "--max-node-span", "12",
        "--shard-size", "3", "--max-indel-len", "10", "--chromosomes", "autosome", "--chr-index", str(table),
        "--debug-rows"])]


BUILD_CASES = {"G1": g1, "G2": g2, "G3": g3, "G4a": lambda i: g4(i, False), "G4b": lambda i: g4(i, True),
               "G5": g5, "G6": g6, "G7": g7, "G8": g8}
CASES = (*BUILD_CASES, "O1")


# --- O1: an orchestrated run ----------------------------------------------------------

def write_vcf(path, contigs, records):
    """records = [(chrom, pos, ref, alt, filter, gt)], chromosome- and position-sorted."""
    header = ["##fileformat=VCFv4.2", *(f"##contig=<ID={c},length={len(s)}>" for c, s in contigs.items()),
              '##FILTER=<ID=GAP1,Description="gap">', '##FORMAT=<ID=GT,Number=1,Type=String,Description="GT">',
              "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tS"]
    body = [f"{c}\t{pos}\t.\t{ref}\t{alt}\t.\t{filt}\t.\tGT\t{gt}" for c, pos, ref, alt, filt, gt in records]
    Path(path).write_text("\n".join(header + body) + "\n")
    return Path(path)


def mini_world(inputs, root):
    """Inputs of O1 and its prepare arguments.

    chr1 = nodes 1-6 (12 bp each, forward), off-reference node 7 (walk HG1 >1>7>3); chr2 = nodes 8-12
    with the displaced repeat 8(..C) 9(AT) 10(AT) 11(G..); chrX = nodes 13-15. Targets: every node but 9,
    so the chr2 insertion (placed at node 10 by forward records) is displaced onto node 9 and built by
    a supplement task; the chrX targets are removed by --chromosomes autosome. Truth: germline SNV
    (node 2) and DEL (node 3), somatic SNV (node 4), the chr2 insertion and an SNV on node 12; SNVs on
    nodes 1, 5 and 6 are non-truth (label non) and the node-7 SNV is off GRCh38 (label ignore).
    """
    rng = random.Random(20260925)

    def bases(k):
        return "".join(rng.choice("ACGT") for _ in range(k))

    seq = {n: bases(12) for n in range(1, 7)}
    seq[3] = seq[3][:3] + "A" + "C" + "G" + seq[3][6:]  # the node-3 DEL (offset 4) cannot shift
    seq.update({7: bases(6), 8: bases(9) + "C", 9: "AT", 10: "AT", 11: "G" + bases(9), 12: bases(10)})
    seq.update({n: bases(10) for n in (13, 14, 15)})
    contigs = {"chr1": "".join(seq[n] for n in range(1, 7)), "chr2": "".join(seq[n] for n in range(8, 13)),
               "chrX": "".join(seq[n] for n in (13, 14, 15))}

    def alt_of(base):
        return {"A": "C", "C": "G", "G": "T", "T": "A"}[base]

    def full(node, edits=None):
        return (node, 0, False, edits or [(len(seq[node]), len(seq[node]), "")])

    def snv(node, offset):
        rest = len(seq[node]) - offset - 1
        return full(node, [(offset, offset, ""), (1, 1, alt_of(seq[node][offset])), (rest, rest, "")])

    rows = []
    chr1 = {"A": (1, range(0, 4)), "B": (2, range(0, 10)), "D": (4, range(10, 15)), "E": (5, range(12, 16)),
            "F": (6, range(15, 19))}
    for i in range(20):
        specs = []
        for node in range(1, 7):
            carried = [k for k, (n, reads) in chr1.items() if n == node and i in reads]
            if carried:
                specs.append(snv(node, 4))
            elif node == 3 and 5 <= i < 10:
                specs.append(full(3, [(4, 4, ""), (1, 0, ""), (7, 7, "")]))
            else:
                specs.append(full(node))
        rows.append(spec_alignment(mirror(specs, seq) if i >= 16 else specs, seq, name=f"chr1_{i:02d}"))
    for i in range(5):
        rows.append(spec_alignment([full(1), snv(7, 2) if i < 4 else full(7), full(3)], seq, name=f"alt_{i}"))
    right = [full(8), full(9), full(10, [(2, 2, ""), (0, 2, "AT")]), full(11), full(12)]
    left = [full(8), full(9, [(0, 2, "AT"), (2, 2, "")]), full(10), full(11), full(12)]
    for i in range(10):
        spec = right if i < 3 else left if i < 6 else [full(n) for n in range(8, 13)]
        if i in (0, 1, 6, 7):
            spec = spec[:4] + [snv(12, 3)]
        rows.append(spec_alignment(mirror(spec, seq) if 3 <= i < 6 else spec, seq, name=f"chr2_{i}"))
    for i in range(5):
        rows.append(spec_alignment([full(13), snv(14, 5) if i < 4 else full(14), full(15)], seq, name=f"chrX_{i}"))

    gam = write_gam(inputs / "mini.gam", rows)
    graph = graph_fixture(inputs / "graph.sqlite", [(n, s, 1 + n % 4) for n, s in seq.items()])
    nodes = write_nodes(inputs / "nodes.txt", [n for n in seq if n != 9])
    stats = inputs / "node_stats.json"
    stats.write_text(json.dumps({str(n): dict(perfect=10, not_perfect=n, max_read_length=72) for n in seq},
                                indent=2) + "\n")
    table = chr_table(inputs / "chr.tsv", [("chr1", 1, 7, "autosome"), ("chr2", 8, 12, "autosome"),
                                           ("chrX", 13, 15, "non_autosomal")])
    lines = ["H\tVN:Z:1.1"] + [f"S\t{n}\t{s}" for n, s in seq.items()]
    lines += [f"W\tGRCh38\t0\tchr1\t0\t{len(contigs['chr1'])}\t>1>2>3>4>5>6",
              f"W\tGRCh38\t0\tchr2\t0\t{len(contigs['chr2'])}\t>8>9>10>11>12",
              f"W\tGRCh38\t0\tchrX\t0\t{len(contigs['chrX'])}\t>13>14>15",
              f"W\tHG1\t1\tctg1\t0\t{len(seq[1]) + len(seq[7]) + len(seq[3])}\t>1>7>3"]
    (inputs / "mini.gfa").write_text("\n".join(lines) + "\n")
    scan(inputs / "mini.gfa", inputs / "rp")
    meta = json.loads((inputs / "rp" / "meta.json").read_text())
    meta.update(source=dict(path="mini.gfa", size=meta["source"]["size"], mtime_ns=0), scan_seconds=0.0)
    (inputs / "rp" / "meta.json").write_text(json.dumps(meta, indent=2) + "\n")
    fasta = inputs / "mini.fa"
    fasta.write_text("".join(f">{c}\n{s}\n" for c, s in contigs.items()))
    pysam.faidx(str(fasta))
    c1, c2 = contigs["chr1"], contigs["chr2"]
    germline = write_vcf(inputs / "germline.vcf", contigs, [
        ("chr1", 17, c1[16], alt_of(c1[16]), "PASS", "0|1"), ("chr1", 28, c1[27:29], c1[27], "PASS", "0|1")])
    somatic = write_vcf(inputs / "somatic.vcf", contigs, [
        ("chr1", 41, c1[40], alt_of(c1[40]), "PASS", "0|1"), ("chr2", 10, c2[9], c2[9] + "AT", "PASS", "0|1"),
        ("chr2", 28, c2[27], alt_of(c2[27]), "PASS", "0|1")])
    bed = "".join(f"{c}\t0\t{len(s)}\n" for c, s in contigs.items() if c != "chrX")
    (inputs / "somatic.bed").write_text(bed)
    (inputs / "germline.bed").write_text(bed)
    return ["--root", str(root), "--gam", str(gam), "--nodes", str(nodes), "--node-stats", str(stats),
            "--graph-index", str(graph), "--tasks", "3", "--processes", "2", "--supplement-rounds", "3",
            "--supplement-min-records", "1", "--snv-min-af", "0.06", "--indel-min-af", "0.08", "--gam-cache-mb", "1",
            "--batch-nodes", "2", "--shard-size", "2", "--chromosomes", "autosome", "--chr-index", str(table),
            "--merge-shard-size", "4", "--keep-sources", "--reference-path", str(inputs / "rp"),
            "--somatic-vcf", str(somatic), "--somatic-bed", str(inputs / "somatic.bed"),
            "--germline-vcf", str(germline), "--germline-bed", str(inputs / "germline.bed"),
            "--reference-fasta", str(fasta), "--truth-dir", str(root / "truth"), "--decoder", "auto"]


def o1_coverage(root):
    """The O1 properties the goldens rely on; returns a list of the missing ones."""
    config = json.loads((root / "config.json").read_text())
    missing = []
    if not config["supplement"]["rounds"] or not config["supplement"]["rounds"][0]["nodes"]:
        missing.append("a supplement node")
    if not config["chromosome_selection"].get("removed"):
        missing.append("a removed target")
    blocks, labels = set(), set()
    for kind in config["variant_outputs"]:
        merged = json.loads((root / "tensors" / kind / "manifest.json").read_text())
        blocks |= {c for c, info in merged["chromosomes"].items() if info["tensors"]}
        labels |= set(json.loads((root / "tensors" / kind / "labels.manifest.json").read_text())["totals"])
    if len(blocks) < 2:
        missing.append("merged tensors in two autosome blocks")
    missing += [f"a {name} label" for name in ("somatic", "germline", "non", "ignore") if name not in labels]
    return missing


# --- running --------------------------------------------------------------------------

def environment(**extra):
    env = {k: v for k, v in os.environ.items() if k not in ("PANSOMA_DECODER", "PYTHONHASHSEED", "SLURM_CPUS_PER_TASK")}
    env.update(PYTHONDONTWRITEBYTECODE="1", PYTHONUNBUFFERED="1", OMP_NUM_THREADS="1", OPENBLAS_NUM_THREADS="1",
               MKL_NUM_THREADS="1", NUMEXPR_NUM_THREADS="1")
    env.update({k: str(v) for k, v in extra.items()})
    return env


def execute(command, cwd, env):
    result = subprocess.run(command, cwd=cwd, env=env, capture_output=True, text=True)
    if result.returncode:
        raise RuntimeError(f"{' '.join(map(str, command))} exited with {result.returncode}:\n"
                           f"{result.stdout[-3000:]}\n{result.stderr[-3000:]}")
    return result


def run_job(package, case, decoder, work, hashseed=None):
    """Generate a case's inputs under work/<label>/inputs, build them with `package`, fingerprint the outputs."""
    label = f"{case}-{decoder}" + (f"-seed{hashseed}" if hashseed is not None else "")
    directory = (Path(work) / label).resolve()
    inputs = directory / "inputs"
    inputs.mkdir(parents=True)
    env = environment(**({"PYTHONHASHSEED": hashseed} if hashseed is not None else {}))
    python = sys.executable
    if case == "O1":
        tree = directory / "run"
        prepare = mini_world(inputs, tree) + ["--max-batch-alignments", "20000"]  # v2's default (v3's is 200000)
        freeze_times(inputs)
        env.update(SLURM_CPUS_PER_TASK="2", **({"PANSOMA_DECODER": "python"} if decoder == "python" else {}))
        execute([python, "-m", f"{package}.orchestrate", "prepare", *prepare], REPO, env)
        execute(["bash", str(tree / "run.sh")], REPO, env)
        manifests = sorted((tree / "tensors" / "shared").glob("task_*/manifest.json"))
    else:
        tree = directory / "out"
        builds = BUILD_CASES[case](inputs)
        freeze_times(inputs)
        for name, arguments in builds:
            if "--snv-min-af" not in arguments:
                arguments = arguments + DEFAULT_AF
            if "--batch-nodes" not in arguments:  # v2's default, which the goldens were recorded with
                arguments = arguments + ["--batch-nodes", "512"]
            if "--max-batch-alignments" not in arguments:  # v2's default (v3's is 200000)
                arguments = arguments + ["--max-batch-alignments", "20000"]
            out = tree / name
            execute([python, "-m", f"{package}.run", "build", *arguments, "--output", str(out / "shared"),
                     "--snv-output", str(out / "SNV"), "--indel-output", str(out / "INDEL"), "--decoder", decoder],
                    REPO, env)
        manifests = sorted(tree.glob("*/shared/manifest.json"))
    used = sorted({json.loads(p.read_text())["decoder"]["used"] for p in manifests})
    return dict(label=label, case=case, decoder=decoder, hashseed=hashseed, directory=str(directory),
                inputs=input_hashes(inputs), outputs=fingerprint(tree, anchor=directory), used=used,
                o1_missing=o1_coverage(tree) if case == "O1" else [])


def jobs(cases=CASES, decoders=DECODERS, hash_seeds=True):
    """(case, decoder, hashseed) of a full check: every case per decoder, plus G8 per PYTHONHASHSEED."""
    result = [(case, decoder, None) for case in cases for decoder in decoders]
    if hash_seeds and "G8" in cases:
        result += [("G8", decoder, seed) for decoder in decoders for seed in HASH_SEEDS]
    return result


def run_jobs(package, selected, work, workers=None):
    """{label: result dict or the exception it raised}, run in parallel (each job is a subprocess chain)."""
    cpus = len(os.sched_getaffinity(0)) if hasattr(os, "sched_getaffinity") else (os.cpu_count() or 2)
    workers = workers or max(1, min(len(selected), 16, cpus // 2))

    def one(job):
        case, decoder, seed = job
        label = f"{case}-{decoder}" + (f"-seed{seed}" if seed is not None else "")
        try:
            return label, run_job(package, case, decoder, work, seed)
        except Exception as exc:  # noqa: BLE001 -- reported per job
            return label, exc

    with ThreadPoolExecutor(max_workers=workers) as pool:
        return dict(pool.map(one, selected))


def native_module(package):
    """(True, path) when a compiled _fastdecode sits in the package, else (False, reason)."""
    found = sorted((REPO / package).glob("_fastdecode*.so"))
    return (True, str(found[0])) if found else (False, f"no compiled _fastdecode*.so in {package}")


def constants(package):
    """The package's output constants, dumped by a subprocess (never imported here)."""
    return json.loads(execute([sys.executable, "-c", CONSTANTS_SNIPPET, package], REPO, environment()).stdout)


def compare_constants(recorded, current):
    """Differences of the output constants; BUILDER_OPTIONS must be v2's minus candidate_unit (v3 drops it)."""
    problems = []
    for key in sorted(recorded.keys() | current.keys()):
        a, b = recorded.get(key), current.get(key)
        if key == "orchestrate.BUILDER_OPTIONS" and a:
            a = [k for k in a if k != "candidate_unit"]
        if a != b:
            problems.append(f"constant {key}: recorded {a!r}, now {b!r}")
    return problems


def versions():
    import google.protobuf
    import numpy
    result = dict(python=sys.version.split()[0], executable=sys.executable, numpy=numpy.__version__,
                  protobuf=google.protobuf.__version__, pysam=pysam.__version__)
    for name in ("scipy", "pybind11"):
        try:
            result[name] = importlib.import_module(name).__version__
        except ImportError:
            result[name] = None
    return result


def git(*arguments):
    result = subprocess.run(["git", "-C", str(REPO), *arguments], capture_output=True, text=True, check=True)
    return result.stdout.strip()


def provenance(package):
    folder = REPO / package
    if git("status", "--porcelain", "--", package):
        raise ValueError(f"{package} has uncommitted changes; record goldens only from a committed package")
    files = sorted(p for p in folder.rglob("*") if p.is_file() and p.suffix in (".py", ".cpp")
                   and not {"__pycache__", "tests", "examples"} & set(p.relative_to(folder).parts))
    return dict(package=package, commit=git("log", "-1", "--format=%H", "--", package),
                source_sha256={p.relative_to(folder).as_posix(): hashlib.sha256(p.read_bytes()).hexdigest()
                               for p in files}, **versions())


def load(path=HASHES):
    return json.loads(Path(path).read_text())


def check_results(recorded, results):
    """Problems of each job's result against the recorded goldens ({label: [problem, ...]})."""
    problems = {}
    for label, result in results.items():
        if isinstance(result, Exception):
            problems[label] = [f"failed: {result}"]
            continue
        golden = recorded["cases"][result["case"]]
        found = []
        if result["inputs"] != golden["inputs"]:
            found.append("fixture inputs differ from the recorded ones (" + "; ".join(
                differences(golden["inputs"], result["inputs"])) + "): the generators changed; re-record the "
                "goldens from v2 (python -m " + PACKAGE + ".tests.golden record --package indexed_gam_pipeline_v2)")
        if result["used"] != [result["decoder"]]:
            found.append(f"decoder used {result['used']}, expected {result['decoder']}")
        found += differences(golden["outputs"], result["outputs"])
        found += [f"O1 lacks {m}" for m in result["o1_missing"]]
        if found:
            problems[label] = found
    return problems


def environment_drift(recorded):
    now = versions()
    return [f"{k}: recorded {v}, now {now.get(k)}" for k, v in recorded["recorded_with"].items()
            if k in now and now.get(k) != v]


def record(package, output=HASHES, workers=None):
    """Build every case with `package` under both decoders and write its fingerprints to `output`."""
    ok, reason = native_module(package)
    if not ok:
        raise ValueError(f"recording needs the native decoder of {package}: {reason}")
    with tempfile.TemporaryDirectory(prefix="golden-record-") as work:
        results = run_jobs(package, jobs(), work, workers)
        failed = {label: str(r) for label, r in results.items() if isinstance(r, Exception)}
        if failed:
            raise RuntimeError("recording failed:\n" + "\n".join(f"{k}: {v}" for k, v in failed.items()))
        problems = []
        for label, result in results.items():
            if result["used"] != [result["decoder"]]:
                problems.append(f"{label}: decoder used {result['used']}")
            reference = results[f"{result['case']}-python"]
            for key in ("inputs", "outputs"):
                found = differences(reference[key], result[key])
                if found:
                    problems.append(f"{label} {key} differ from {result['case']}-python: {found}")
            problems += [f"{label}: O1 lacks {m}" for m in result["o1_missing"]]
        if problems:
            raise RuntimeError("inconsistent recording:\n" + "\n".join(problems))
    cases = {case: dict(inputs=results[f"{case}-python"]["inputs"], outputs=results[f"{case}-python"]["outputs"])
             for case in CASES}
    document = dict(format="golden-hashes-v1", recorded_with=provenance(package), constants=constants(package),
                    cases=cases)
    Path(output).write_text(json.dumps(document, indent=1) + "\n")
    return document


def check(package, decoders=DECODERS, cases=CASES, keep=None, workers=None, recorded=None):
    """Problems of `package` against the recorded goldens: {'constants' or job label: [problem, ...]}."""
    recorded = recorded or load()
    problems = {}
    found = compare_constants(recorded["constants"], constants(package))
    if found:
        problems["constants"] = found
    work = Path(keep) if keep else Path(tempfile.mkdtemp(prefix="golden-check-"))
    work.mkdir(parents=True, exist_ok=True)
    if any(work.iterdir()):
        raise ValueError(f"--keep directory must be new or empty: {work}")
    try:
        results = run_jobs(package, jobs(cases, decoders), work, workers)
        problems.update(check_results(recorded, results))
    finally:
        if not keep:
            shutil.rmtree(work, ignore_errors=True)
    return problems, results


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    commands = parser.add_subparsers(dest="command", required=True)
    r = commands.add_parser("record", help="record the goldens from a package (v2, before any trimming)")
    r.add_argument("--package", required=True)
    r.add_argument("--output", default=str(HASHES))
    c = commands.add_parser("check", help="rebuild every case with a package and compare with the goldens")
    c.add_argument("--package", default=PACKAGE)
    c.add_argument("--decoder", choices=DECODERS, action="append", help="default: both")
    c.add_argument("--case", nargs="+", choices=CASES, default=list(CASES))
    c.add_argument("--keep", help="build under this new or empty directory and keep it")
    for p in (r, c):
        p.add_argument("--workers", type=int, help="parallel jobs (default: min(16, allocated CPUs / 2))")
    args = parser.parse_args(argv)
    if args.command == "record":
        document = record(args.package, args.output, args.workers)
        print(f"recorded {len(document['cases'])} cases from {args.package} "
              f"({document['recorded_with']['commit'][:7]}) into {args.output}")
        return
    decoders = args.decoder or list(DECODERS)
    ok, reason = native_module(args.package)
    if "native" in decoders and not ok:
        print(f"native cases skipped: {reason}")
        decoders = [d for d in decoders if d != "native"]
    recorded = load()
    problems, results = check(args.package, decoders, args.case, args.keep, args.workers, recorded)
    for label in sorted(results):
        result = results[label]
        files = len(result["outputs"]) if isinstance(result, dict) else 0
        print(f"{label:16s} {'MISMATCH' if label in problems else 'ok'} ({files} files)")
    for label, found in sorted(problems.items()):
        print(f"--- {label}")
        for line in found:
            print("  " + line)
    if problems:
        for line in environment_drift(recorded):
            print("environment differs from the recording: " + line)
        sys.exit(1)
    print(f"all {len(results)} golden jobs match ({args.package}, decoders {', '.join(decoders)})")


if __name__ == "__main__":
    main()
