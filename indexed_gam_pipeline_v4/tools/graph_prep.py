"""One-time graph preparation: the GFA, the reference-path directory, the chromosome block table and the graph
index audit (the graph index itself: tools.graph_index_build; all steps as one job: tools/jobs/graph_prep.sh).

    python -m indexed_gam_pipeline_v4.tools.graph_prep gfa --gbz G.gbz --output G.gfa [--threads 16] [--vg VG]
    python -m indexed_gam_pipeline_v4.tools.graph_prep ref-path-scan --gfa G.gfa --output DIR [--reference-sample GRCh38]
    python -m indexed_gam_pipeline_v4.tools.graph_prep ref-path-check --path DIR --graph-index DB --fasta FA [--samples 100000]
    python -m indexed_gam_pipeline_v4.tools.graph_prep components --gbz G.gbz --reference-path DIR --output DIR [--vg VG]
    python -m indexed_gam_pipeline_v4.tools.graph_prep chr-index --components-dir D --reference-path DIR --output PREFIX [--graph-index DB]
    python -m indexed_gam_pipeline_v4.tools.graph_prep audit --graph-index DB --gfa G.gfa --chr-index TSV --output JSON [--processes 16]

The formats are documented with their run-time readers: tensor_postprocessing.reference_path
(gfa-reference-path directory, ReferencePath) and tensor_postprocessing.chr_index (chr-node-ranges
TSV plus JSON with tsv_sha256, ChrIndex). scan writes the reference-path FORMAT, build the
CHR_INDEX_FORMAT; nothing here is imported at run time. `gfa` and `components` run vg: --vg, else $PANSOMA_VG
(scripts/use_vg.sh), else vg on PATH.
"""
import argparse
import contextlib
import csv
import hashlib
import json
import os
from pathlib import Path
import random
import subprocess
import time

import numpy as np

from ..common import read_json, stamp, write_json
from ..tensor_postprocessing.chr_index import AUTOSOMES, FIELDS, ChrIndex
from ..tensor_postprocessing.reference_path import FORMAT, ReferencePath, rc

# --- GFA from the GBZ ------------------------------------------------------------------------------

def vg_command(vg=None):
    """--vg, else $PANSOMA_VG, else vg on PATH."""
    return vg or os.environ.get("PANSOMA_VG") or "vg"


def vg_version(vg):
    return subprocess.run([vg, "version"], capture_output=True, text=True, check=True).stdout.splitlines()[0]


def to_gfa(gbz, output, vg=None, threads=1):
    """GBZ -> GFA (vg convert -f --no-translation): segment names are the GBZ node IDs, the IDs in GAM alignments
    and in the graph index, not the original segment names a translation would restore. Published atomically."""
    vg, output, started = vg_command(vg), Path(output), time.perf_counter()
    if output.exists():
        raise ValueError(f"Output exists: {output}")
    work = output.with_name(output.name + ".tmp")
    with work.open("wb") as stream:
        subprocess.run([vg, "convert", "-f", "--no-translation", "-t", str(threads), str(gbz)], stdout=stream,
                       check=True)
    with work.open("rb") as stream:
        if stream.read(2) != b"H\t":
            raise ValueError(f"vg convert wrote no GFA header to {work}")
    work.rename(output)
    return dict(gbz=stamp(gbz), gfa=stamp(output), vg=vg_version(vg), seconds=time.perf_counter() - started)


# --- reference path: one pass over a GFA --------------------------------------------------------

SEPARATORS = bytes.maketrans(b"<>", b"  ")
PATH_SEPARATORS = bytes.maketrans(b"+-,", b"   ")
# GBWT's sample name of generic paths (plain names such as chr1, written as P lines by vg convert -f): the reference of
# a GBZ built from a FASTA alone, e.g. a GRCh38-only vg autoindex graph. ref-path-scan --reference-sample _gbwt_ref
# takes the GFA's P lines as the reference walks.
GENERIC_SAMPLE = "_gbwt_ref"
# S lines -> "id<TAB>length" into LENGTHS (LN:i: when the sequence is '*'); W lines pass through, P lines too when
# GENERIC is 1.
AWK = r'''BEGIN { FS = OFS = "\t" }
$1 == "S" { n = length($3)
            if ($3 == "*") { n = -1; for (i = 4; i <= NF; i++) if ($i ~ /^LN:i:/) n = substr($i, 6) + 0 }
            print $2, n > LENGTHS; next }
$1 == "W" { print }
$1 == "P" && GENERIC { print }'''


def parse_walk(walk):
    """b'>12<3' -> (ids int64, reverse bool)."""
    raw = np.frombuffer(walk, dtype=np.uint8)
    reverse = raw[(raw == 60) | (raw == 62)] == 60
    ids = np.array(walk.translate(SEPARATORS).split(), dtype=np.int64)
    if ids.size != reverse.size or not ids.size:
        raise ValueError("Malformed GFA walk")
    return ids, reverse


def parse_path(segments):
    """b'12+,3-' (the segment field of a GFA P line) -> (ids int64, reverse bool)."""
    raw = np.frombuffer(segments, dtype=np.uint8)
    reverse = raw[(raw == 43) | (raw == 45)] == 45
    ids = np.array(segments.translate(PATH_SEPARATORS).split(), dtype=np.int64)
    if ids.size != reverse.size or not ids.size:
        raise ValueError("Malformed GFA path")
    return ids, reverse


def reference_path_name(meta, contig):
    """The GBZ path name of a reference contig: <sample>#<hap>#<contig>, or the contig itself for generic paths."""
    if meta["reference_sample"] == GENERIC_SAMPLE:
        return contig
    hap = next(c["hap"] for c in meta["contigs"] if c["name"] == contig)
    return f"{meta['reference_sample']}#{hap}#{contig}"


def scan(gfa, output, reference_sample="GRCh38"):
    """Build the output directory from one streaming pass (awk prefilter) over `gfa`. The reference walks are the W
    lines of `reference_sample`, or with GENERIC_SAMPLE the P lines (generic paths: hap 0, start 0)."""
    started = time.perf_counter()
    output = Path(output)
    if output.exists():
        raise ValueError(f"Output exists: {output}")
    work = output.with_name(output.name + ".tmp")
    work.mkdir(parents=True)
    lengths_txt = work / "segment_lengths.txt"
    # Walk records are streamed to walks.ndjson, never kept: a graph can have hundreds of millions of W lines.
    reference, n_walks, samples = [], 0, set()
    generic, paths = reference_sample == GENERIC_SAMPLE, []  # paths: P lines, logged once their lengths are known
    process = subprocess.Popen(["awk", "-v", f"LENGTHS={lengths_txt}", "-v", f"GENERIC={int(generic)}", AWK, str(gfa)],
                               stdout=subprocess.PIPE)
    with (work / "walks.ndjson").open("w") as log:
        for line in process.stdout:
            fields = line.rstrip(b"\n").split(b"\t")
            if fields[0] == b"P":
                if len(fields) < 3:
                    raise ValueError("Malformed GFA P line")
                ids, reverse = parse_path(fields[2])
                paths.append((dict(sample=GENERIC_SAMPLE, hap="0", contig=fields[1].decode(), start=0, end=None,
                                   nodes=int(ids.size), min=int(ids.min()), max=int(ids.max())), ids, reverse))
                continue
            if len(fields) != 7:
                raise ValueError("Malformed GFA W line")
            sample, hap, contig = (f.decode() for f in fields[1:4])
            if sample == reference_sample:
                ids, reverse = parse_walk(fields[6])
            else:  # only the node count and ID range are recorded
                ids = np.fromstring(fields[6].translate(SEPARATORS), dtype=np.int64, sep=" ")
                if not ids.size:
                    raise ValueError("Malformed GFA walk")
            record = dict(sample=sample, hap=hap, contig=contig, start=int(fields[4]), end=int(fields[5]),
                          nodes=int(ids.size), min=int(ids.min()), max=int(ids.max()))
            log.write(json.dumps(record) + "\n")
            n_walks += 1
            samples.add(sample)
            if sample == reference_sample:
                reference.append((record, ids, reverse))
    process.stdout.close()
    if process.wait():
        raise RuntimeError(f"awk failed on {gfa}")
    if not reference and not paths:
        raise ValueError(f"No {'P' if generic else 'W'} lines for sample {reference_sample} in {gfa}")
    pairs = np.fromfile(lengths_txt, sep=" ", dtype=np.int64).reshape(-1, 2)
    if (pairs[:, 1] < 0).any() or (pairs[:, 0] <= 0).any():
        raise ValueError("Segments need a positive integer ID and a sequence or LN tag")
    size = int(pairs[:, 0].max()) + 1
    lengths = np.zeros(size, dtype=np.int32)
    lengths[pairs[:, 0]] = pairs[:, 1]
    if np.count_nonzero(lengths) != len(pairs):
        raise ValueError("Duplicate or empty segments in the GFA")
    with (work / "walks.ndjson").open("a") as log:  # a generic path spans its whole contig: end = its length
        for record, ids, reverse in paths:
            if ids.max() >= size or not lengths[ids].all():
                raise ValueError(f"Reference path {record['contig']} uses nodes without segments")
            record["end"] = int(lengths[ids].astype(np.int64).sum())
            log.write(json.dumps(record) + "\n")
            n_walks += 1
            samples.add(GENERIC_SAMPLE)
            reference.append((record, ids, reverse))
    chrom = np.full(size, -1, dtype=np.int16)
    start0 = np.zeros(size, dtype=np.int64)
    reverse_of = np.zeros(size, dtype=bool)
    visits = np.zeros(size, dtype=np.int64)
    contigs, offset, order = [], 0, []
    for k, (record, ids, reverse) in enumerate(reference):
        if ids.max() >= size or not lengths[ids].all():
            raise ValueError(f"Reference walk {record['contig']} uses nodes without segments")
        span = lengths[ids].astype(np.int64)
        starts = record["start"] + np.cumsum(span) - span
        if record["start"] + int(span.sum()) != record["end"]:
            raise ValueError(f"Reference walk {record['contig']}: node lengths disagree with the W line interval")
        unique, first, counts = np.unique(ids, return_index=True, return_counts=True)
        visits[unique] += counts
        new = unique[chrom[unique] == -1]
        firsts = first[chrom[unique] == -1]
        chrom[new], start0[new], reverse_of[new] = k, starts[firsts], reverse[firsts]
        contigs.append(dict(name=record["contig"], hap=record["hap"], start=record["start"], end=record["end"],
                            nodes=int(ids.size), min_node=record["min"], max_node=record["max"], offset=offset))
        offset += ids.size
        order.append((ids, starts, reverse))
    arrays = dict(lengths=lengths, chrom=chrom, start0=start0, reverse=reverse_of,
                  visits=np.minimum(visits, np.iinfo(np.uint16).max).astype(np.uint16),
                  path_nodes=np.concatenate([o[0] for o in order]), path_starts=np.concatenate([o[1] for o in order]),
                  path_reverse=np.concatenate([o[2] for o in order]))
    for name, values in arrays.items():
        np.save(work / f"{name}.npy", values)
    lengths_txt.unlink()
    meta = dict(format=FORMAT, source=stamp(gfa), reference_sample=reference_sample, max_node=size - 1,
                segments=int(len(pairs)), walks=n_walks, samples=len(samples),
                contigs=contigs, reference_nodes=int(np.count_nonzero(visits)),
                ambiguous_reference_nodes=int(np.count_nonzero(visits > 1)),
                scan_seconds=time.perf_counter() - started, checks={})
    write_json(work / "meta.json", meta)
    work.rename(output)
    return meta


def check(directory, graph_index, fasta, samples=100000, seed=20260923):
    """GFA lengths vs the GBZ graph index, and reference node sequences vs the FASTA; recorded in meta.json."""
    from ..graph_index import GraphIndex
    import pysam
    path = ReferencePath(directory)
    rng = random.Random(seed)
    lengths = np.asarray(path.lengths)
    present = np.flatnonzero(lengths)
    any_nodes = [int(present[rng.randrange(len(present))]) for _ in range(samples)]
    reference = np.flatnonzero(np.asarray(path.visits) == 1)
    ref_nodes = [int(reference[rng.randrange(len(reference))]) for _ in range(samples)]
    report = dict(graph_index=str(graph_index), fasta=str(fasta), samples=samples, seed=seed)
    with GraphIndex(graph_index) as graph, pysam.FastaFile(str(fasta)) as genome:
        report["graph_index_nodes"] = graph.metadata["nodes"]
        report["gfa_segments"] = path.meta["segments"]
        report["gfa_max_node"] = path.meta["max_node"]  # above the node count when IDs have gaps (filtered graphs)
        records = graph.get_nodes(set(any_nodes) | set(ref_nodes))
        report["length_mismatches"] = int(sum(len(records[n]["sequence"]) != int(lengths[n]) for n in set(any_nodes)))
        mismatches, skipped = 0, 0
        for n in sorted(set(ref_nodes)):
            contig = path.contigs[path.chrom[n]]
            if contig not in genome.references:
                skipped += 1
                continue
            start = int(path.start0[n])
            linear = genome.fetch(contig, start, start + int(lengths[n])).upper()
            node = records[n]["sequence"].upper()
            mismatches += int((rc(node) if path.reverse[n] else node) != linear)
        report.update(fasta_sequence_mismatches=mismatches, fasta_nodes_checked=len(set(ref_nodes)) - skipped,
                      fasta_contig_missing=skipped)
    report["passed"] = (report["length_mismatches"] == 0 and mismatches == 0
                        and report["graph_index_nodes"] == report["gfa_segments"])
    meta = read_json(Path(directory) / "meta.json")
    meta["checks"]["graph_index_and_fasta"] = report
    write_json(Path(directory) / "meta.json", meta)
    if not report["passed"]:
        raise ValueError(f"Reference path check failed: {report}")
    return report


# --- chromosome components (the chr-index input) ------------------------------------------------

def components(gbz, reference_path, output, vg=None, threads=1, autosomes=AUTOSOMES):
    """<output>/chrN/chrN.component.nodes.raw.txt (one node ID per line) for chr-index: the connected component of
    each autosome's reference path, `vg chunk -C -p <sample>#<hap>#<contig>` (a generic path: `-p <contig>`) on the
    GBZ, off-reference nodes
    included; plus <output>/summary.json (per chromosome: node count, ID range, nodes on the reference path).
    The path names come from the ref-path-scan directory. Published atomically."""
    vg, output, started = vg_command(vg), Path(output), time.perf_counter()
    if output.exists():
        raise ValueError(f"Output exists: {output}")
    path = ReferencePath(reference_path)
    contigs = {c["name"]: c for c in path.meta["contigs"]}
    visits = np.asarray(path.visits)
    work = output.with_name(output.name + ".tmp")
    work.mkdir(parents=True)
    summary = {}
    for chrom in autosomes:
        if chrom not in contigs:
            raise ValueError(f"No reference walk for {chrom} in {reference_path}")
        name = reference_path_name(path.meta, chrom)
        chunk = subprocess.Popen([vg, "chunk", "-x", str(gbz), "-p", name, "-C", "-t", str(threads)],
                                 stdout=subprocess.PIPE)
        convert = subprocess.Popen([vg, "convert", "-f", "-"], stdin=chunk.stdout, stdout=subprocess.PIPE)
        chunk.stdout.close()
        with convert.stdout:
            ids = [int(line.split(b"\t", 2)[1]) for line in convert.stdout if line.startswith(b"S\t")]
        if convert.wait() or chunk.wait():
            raise RuntimeError(f"vg chunk -C -p {name} failed on {gbz}")
        ids = np.array(ids, dtype=np.int64)
        if not ids.size:
            raise ValueError(f"Empty component for {name}")
        (work / chrom).mkdir()
        (work / chrom / f"{chrom}.component.nodes.raw.txt").write_text("".join(f"{n}\n" for n in ids.tolist()))
        on_path = ids[ids < visits.size]
        reference_nodes = int(np.count_nonzero(visits[on_path] > 0))
        summary[chrom] = dict(path=name, nodes=int(ids.size), first_node=int(ids.min()), last_node=int(ids.max()),
                              reference_nodes=reference_nodes, off_reference_nodes=int(ids.size) - reference_nodes)
    meta = dict(gbz=stamp(gbz), reference_path=str(Path(reference_path).resolve()), vg=vg_version(vg),
                chromosomes=summary, seconds=time.perf_counter() - started)
    write_json(work / "summary.json", meta)
    work.rename(output)
    return meta


# --- chromosome block table ---------------------------------------------------------------------

CHR_INDEX_FORMAT = "chr-node-ranges"
NAMED_GROUPS = ("chrX", "chrY", "chrM", "chrEBV")
UNPLACED = "unplaced"


def group_of(contig):
    """Non-autosomal group of a reference contig name."""
    return contig if contig in NAMED_GROUPS else UNPLACED


def component_block(path, graph=None):
    """(first, last, count) of one component node list. Its node-ID interval [first, last] may have gaps (IDs of a
    filtered graph, e.g. HPRC v2.1 d46) only when `graph` (a GraphIndex) shows that every graph node in the interval
    is in the component; without a graph index it must be gap-free."""
    ids = np.fromfile(path, sep=" ", dtype=np.int64)
    if not ids.size:
        raise ValueError(f"Empty component node list: {path}")
    first, last = int(ids.min()), int(ids.max())
    if np.unique(ids).size != ids.size:
        raise ValueError(f"Duplicate node IDs in {path}")
    if ids.size != last - first + 1:
        if graph is None:
            raise ValueError(f"Component {path} is not one contiguous node-ID interval "
                             f"({ids.size} nodes in [{first}, {last}]); give the graph index to allow ID gaps")
        inside = graph_nodes_between(graph, first, last)
        if inside != ids.size:
            raise ValueError(f"Component {path}: {inside - ids.size} graph nodes of other components inside its "
                             f"node-ID interval [{first}, {last}]")
    return first, last, int(ids.size)


def graph_nodes_between(graph, first, last):
    """Number of graph index nodes with first <= ID <= last."""
    return graph.db.execute("SELECT count(*) FROM nodes WHERE node_id BETWEEN ? AND ?", (first, last)).fetchone()[0]


def build(components_dir, reference_path, output, graph_index=None, autosomes=AUTOSOMES):
    """Write <output>.tsv + <output>.json; raises on any violated invariant. With the graph index, node IDs may have
    gaps (filtered graphs): block node counts are counted in the index, and a component interval may have gaps."""
    from ..graph_index import GraphIndex
    with GraphIndex(graph_index) if graph_index else contextlib.nullcontext() as graph:
        return _build(components_dir, reference_path, output, graph_index, graph, autosomes)


def _build(components_dir, reference_path, output, graph_index, graph, autosomes):
    from ..tensor_postprocessing.reference_path import ReferencePath
    components_dir, output = Path(components_dir), Path(output)
    path = ReferencePath(reference_path)
    max_node = path.meta["max_node"]
    blocks, sources = [], {}
    for chrom in autosomes:
        source = components_dir / chrom / f"{chrom}.component.nodes.raw.txt"
        first, last, count = component_block(source, graph)
        blocks.append(dict(chrom=chrom, first_node=first, last_node=last, nodes=count, dataset="autosome",
                           source=f"vg chunk -C connected component ({source.name})"))
        sources[chrom] = stamp(source)
    blocks.sort(key=lambda b: b["first_node"])
    for a, b in zip(blocks, blocks[1:]):
        if a["last_node"] >= b["first_node"]:
            raise ValueError(f"Overlapping autosome blocks: {a['chrom']} and {b['chrom']}")
    # Non-autosomal groups: every reference contig starts an interval at its smallest node, which runs to the next
    # start (another contig or an autosome block); adjacent intervals of one group are joined. A group whose contigs
    # are consecutive in node-ID order (HPRC graphs) is one interval, an interleaved one (the GRCh38-only vg autoindex
    # graph) one interval per run.
    autosome_starts = [b["first_node"] for b in blocks]
    contigs = sorted(((c["min_node"], group_of(c["name"])) for c in path.meta["contigs"] if c["name"] not in autosomes))
    intervals = []
    for i, (first, group) in enumerate(contigs):
        if any(b["first_node"] <= first <= b["last_node"] for b in blocks):
            raise ValueError(f"Reference contigs of {group} start inside an autosome block")
        following = [s for s in autosome_starts if s > first] + [s for s, _ in contigs[i + 1:] if s > first]
        last = min(following) - 1 if following else max_node
        if intervals and intervals[-1][0] == group and intervals[-1][2] + 1 == first:
            intervals[-1][2] = last
        elif not intervals or intervals[-1][1] != first:  # contigs sharing a smallest node start one interval
            intervals.append([group, first, last])
    for group, first, last in intervals:
        nodes = graph_nodes_between(graph, first, last) if graph else last - first + 1
        blocks.append(dict(chrom=group, first_node=first, last_node=last, nodes=nodes,
                           dataset="non_autosomal", source="reference contigs' smallest node .. next block"))
    blocks.sort(key=lambda b: b["first_node"])
    index = ChrIndex.from_blocks(blocks)
    # Invariants: every reference contig inside its own block, every walk inside one block.
    for contig in path.meta["contigs"]:
        expected = contig["name"] if contig["name"] in autosomes else group_of(contig["name"])
        got = index.names_of([contig["min_node"], contig["max_node"]])
        if list(got) != [expected, expected]:
            raise ValueError(f"Reference contig {contig['name']} is not inside block {expected}: {got}")
    walks = dict(total=0, crossing_autosome=0, crossing_non_autosomal=0, unassigned=0, examples=[])
    for walk in path.walks():
        walks["total"] += 1
        a, b = index.lookup([walk["min"], walk["max"]])
        if a == b and a >= 0:
            continue
        if a < 0 or b < 0:
            walks["unassigned"] += 1
        elif "autosome" in (index.dataset[a], index.dataset[b]):
            walks["crossing_autosome"] += 1
        else:
            walks["crossing_non_autosomal"] += 1
        if len(walks["examples"]) < 20:
            walks["examples"].append(walk)
    if walks["crossing_autosome"] or walks["unassigned"]:
        raise ValueError(f"Walks leave their chromosome block: {walks}")
    covered = sum(b["nodes"] for b in blocks)
    total = graph.metadata["nodes"] if graph else max_node  # the same when node IDs are 1..N
    meta = dict(format=CHR_INDEX_FORMAT, max_node=max_node, covered_nodes=covered, uncovered_nodes=total - covered,
                autosome_components=sources, reference_path=dict(path=str(Path(reference_path).resolve()),
                                                                 source=path.meta["source"]),
                walk_check=walks)
    if graph:
        meta["graph_index"] = dict(path=str(Path(graph_index).resolve()), nodes=graph.metadata["nodes"],
                                   gbz=graph.metadata["source"])
        if graph.metadata["nodes"] != path.meta["segments"]:
            raise ValueError("Graph index and GFA disagree on the node count")
    output.parent.mkdir(parents=True, exist_ok=True)
    table = output.with_name(output.name + ".tsv")
    with table.open("w", newline="") as stream:
        writer = csv.DictWriter(stream, FIELDS, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(blocks)
    meta["tsv_sha256"] = hashlib.sha256(table.read_bytes()).hexdigest()
    write_json(output.with_name(output.name + ".json"), meta)
    return blocks, meta


# --- graph index audit ---------------------------------------------------------------------------

WALK_SEPARATORS = bytes.maketrans(b"<>", b"  ")      # W line walk: >12<13
PATH_SEPARATORS = bytes.maketrans(b",+-", b"   ")    # P line segments: 12+,13-


def audit_part(gfa, start, end, size, targets):
    """One byte range of the GFA (the lines that start in [start, end)): every node's number of W/P lines visiting it
    (numpy fancy-index increments count a node once per line), the S lines of `targets`, and totals."""
    counts = np.zeros(size, dtype=np.int32)
    sequences, totals = {}, dict(s_lines=0, paths=0, path_visits=0, beyond_index=0)
    with open(gfa, "rb", buffering=1 << 24) as stream:
        stream.seek(max(start - 1, 0))
        if start:
            stream.readline()  # the rest of the line in progress at start - 1 belongs to the previous range
        position = stream.tell()
        while position < end:
            line = stream.readline()
            if not line:
                break
            position += len(line)
            kind = line[:2]
            if kind == b"S\t":
                totals["s_lines"] += 1
                tab = line.index(b"\t", 2)
                node = int(line[2:tab])
                if node in targets:
                    sequences[node] = line[tab + 1:].split(b"\t", 1)[0].rstrip(b"\n").decode()
                continue
            if kind == b"W\t":
                walk = line.split(b"\t", 7)[6].rstrip(b"\n")  # field 7; optional tags would follow
                ids = np.fromstring(walk.translate(WALK_SEPARATORS), dtype=np.int64, sep=" ")
            elif kind == b"P\t":
                ids = np.fromstring(line.split(b"\t", 3)[2].translate(PATH_SEPARATORS), dtype=np.int64, sep=" ")
            else:
                continue
            totals["paths"] += 1
            totals["path_visits"] += int(ids.size)
            if ids.size and (ids.max() >= size or ids.min() < 0):
                totals["beyond_index"] += int(np.count_nonzero((ids >= size) | (ids < 0)))
                ids = ids[(ids < size) & (ids >= 0)]
            counts[ids] += 1
    return counts, sequences, totals


def _audit_part(task):
    return audit_part(*task)


def audit(graph_index, gfa, chr_index, output, random_nodes=24, seed=20260930, processes=None):
    """Independent check of the graph index against the GFA of the same GBZ. distinct_path_count of every index node
    must equal the number of GFA W/P lines (paths) visiting it, recounted in one pass over the GFA (byte ranges in
    `processes` processes; default $SLURM_CPUS_PER_TASK, else 1); every path node must be an index node, and the GFA
    must have one S line per index node. The sequence of the first node of every chr-index block and of
    `random_nodes` random nodes must equal its S line. Writes the report to `output` (JSON); on any difference to
    <output>.failed instead, and raises."""
    from multiprocessing import get_context
    from ..graph_index import GraphIndex
    started = time.perf_counter()
    processes = processes or int(os.environ.get("SLURM_CPUS_PER_TASK") or 1)
    with Path(chr_index).open() as stream:
        blocks = [int(row["first_node"]) for row in csv.DictReader(stream, delimiter="\t")]
    with GraphIndex(graph_index) as index:
        metadata = index.metadata
        low, high = index.db.execute("SELECT min(node_id), max(node_id) FROM nodes").fetchone()
        rng = random.Random(seed)
        picks = {index.db.execute("SELECT node_id FROM nodes WHERE node_id >= ? ORDER BY node_id LIMIT 1",
                                  (rng.randint(low, high),)).fetchone()[0] for _ in range(random_nodes)}
        records = index.get_nodes(set(blocks) | picks)
    total = Path(gfa).stat().st_size
    cuts = [total * k // processes for k in range(processes + 1)]
    tasks = [(str(gfa), a, b, high + 1, set(records)) for a, b in zip(cuts, cuts[1:])]
    counts = np.zeros(high + 1, dtype=np.int64)
    sequences, totals = {}, dict(s_lines=0, paths=0, path_visits=0, beyond_index=0)
    with get_context("fork").Pool(processes) as pool:
        for part_counts, part_sequences, part_totals in pool.imap_unordered(_audit_part, tasks):
            counts += part_counts
            sequences.update(part_sequences)
            for key, value in part_totals.items():
                totals[key] += value
    scan_seconds = time.perf_counter() - started
    present = np.zeros(high + 1, dtype=bool)
    compared, mismatched, examples = 0, 0, []
    with GraphIndex(graph_index) as index:
        cursor = index.db.execute("SELECT node_id, distinct_path_count FROM nodes")
        while True:
            rows = cursor.fetchmany(5_000_000)
            if not rows:
                break
            ids, expected = (np.array(column, dtype=np.int64) for column in zip(*rows))
            present[ids] = True
            wrong = np.flatnonzero(counts[ids] != expected)
            compared += len(ids)
            mismatched += len(wrong)
            examples += [dict(node=int(ids[k]), index_count=int(expected[k]), gfa_paths=int(counts[ids[k]]))
                         for k in wrong[:20 - len(examples)]]
    outside = int(np.count_nonzero(counts[~present]))
    sequence_mismatches = [dict(node=n, in_gfa=n in sequences) for n in sorted(records)
                           if records[n]["sequence"] != sequences.get(n)]
    checks = dict(path_count_mismatches=mismatched, path_nodes_outside_index=outside + totals["beyond_index"],
                  s_lines_minus_index_nodes=totals["s_lines"] - metadata["nodes"],
                  sequence_mismatches=len(sequence_mismatches))
    report = dict(passed=not any(checks.values()), checks=checks, nodes_compared=compared,
                  path_count_examples=examples, sequence_examples=sequence_mismatches[:20],
                  sequence_nodes={str(n): records[n]["distinct_path_count"] for n in sorted(records)},
                  block_first_nodes=blocks, random_nodes=random_nodes, seed=seed, gfa_totals=totals,
                  index_totals={k: metadata.get(k) for k in ("nodes", "logical_paths", "path_visits")},
                  processes=processes, scan_seconds=scan_seconds, seconds=time.perf_counter() - started,
                  graph_index=str(Path(graph_index).resolve()), gfa=stamp(gfa),
                  scope="distinct_path_count of every index node vs the GFA W/P lines visiting it; GFA S lines vs "
                        "index nodes; sequences of the chr-index block starts and random nodes vs their S lines")
    output = Path(output)
    if not report["passed"]:  # the report goes to <output>.failed, so `output` exists only for a passed audit
        write_json(output.with_name(output.name + ".failed"), report)
        raise ValueError(f"Graph index audit failed: {checks}; {examples[:3]} {sequence_mismatches[:3]}")
    write_json(output, report)
    return report


def main(argv=None):
    parser = argparse.ArgumentParser(prog="python -m indexed_gam_pipeline_v4.tools.graph_prep",
                                     description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    commands = parser.add_subparsers(dest="command", required=True)

    p = commands.add_parser("gfa", help="GBZ -> GFA whose segment names are the GBZ node IDs (vg convert)")
    p.add_argument("--gbz", required=True)
    p.add_argument("--output", required=True, help="new .gfa")
    p.add_argument("--threads", type=int, default=1)
    p.add_argument("--vg", help="vg executable (default: $PANSOMA_VG, else vg on PATH)")

    p = commands.add_parser("ref-path-scan", help="one GFA pass: node lengths, GRCh38 path coordinates, walk ranges")
    p.add_argument("--gfa", required=True)
    p.add_argument("--output", required=True, help="new directory")
    p.add_argument("--reference-sample", default="GRCh38",
                   help=f"sample of the reference W lines; {GENERIC_SAMPLE}: the generic paths, the GFA's P lines (a GBZ "
                        "built from a FASTA alone, e.g. a GRCh38-only vg autoindex graph)")

    p = commands.add_parser("ref-path-check", help="check a ref-path directory against the graph index and a FASTA")
    p.add_argument("--path", required=True)
    p.add_argument("--graph-index", required=True)
    p.add_argument("--fasta", required=True)
    p.add_argument("--samples", type=int, default=100000)

    p = commands.add_parser("components", help="connected component node lists of chr1-22 (vg chunk -C)")
    p.add_argument("--gbz", required=True)
    p.add_argument("--reference-path", required=True, help="ref-path-scan output directory (reference path names)")
    p.add_argument("--output", required=True, help="new directory; chr-index's --components-dir")
    p.add_argument("--threads", type=int, default=1)
    p.add_argument("--vg", help="vg executable (default: $PANSOMA_VG, else vg on PATH)")

    p = commands.add_parser("chr-index", help="node ID -> chromosome block table (chr1-22 + non-autosomal groups)")
    p.add_argument("--components-dir", required=True, help="directory with chrN/chrN.component.nodes.raw.txt")
    p.add_argument("--reference-path", required=True, help="ref-path-scan output directory")
    p.add_argument("--output", required=True, help="output prefix; writes <prefix>.tsv and <prefix>.json")
    p.add_argument("--graph-index", help="graph index SQLite (records the GBZ fingerprint, checks the node count)")

    p = commands.add_parser("audit", help="graph index vs the GFA on known nodes (sequence, distinct path count)")
    p.add_argument("--graph-index", required=True)
    p.add_argument("--gfa", required=True, help="the GFA of the index's GBZ (the gfa command's output)")
    p.add_argument("--chr-index", required=True, help="chr-index .tsv (its block starts are audited)")
    p.add_argument("--output", required=True, help="report JSON (e.g. graph_audit.json)")
    p.add_argument("--random-nodes", type=int, default=24, help="sequence check: random nodes besides the block starts")
    p.add_argument("--processes", type=int, help="GFA byte ranges read in parallel (default: $SLURM_CPUS_PER_TASK, else 1)")

    args = parser.parse_args(argv)
    if args.command == "gfa":
        result = to_gfa(args.gbz, args.output, args.vg, args.threads)
    elif args.command == "components":
        result = components(args.gbz, args.reference_path, args.output, args.vg, args.threads)
    elif args.command == "audit":
        result = audit(args.graph_index, args.gfa, args.chr_index, args.output, args.random_nodes,
                       processes=args.processes)
        result = {k: v for k, v in result.items() if k != "sequence_nodes"}
    elif args.command == "ref-path-scan":
        result = scan(args.gfa, args.output, args.reference_sample)
        result = {k: v for k, v in result.items() if k != "contigs"}
    elif args.command == "ref-path-check":
        result = check(args.path, args.graph_index, args.fasta, args.samples)
    elif args.command == "chr-index":
        blocks, meta = build(args.components_dir, args.reference_path, args.output, args.graph_index)
        result = dict(blocks=blocks, walk_check=meta["walk_check"], uncovered_nodes=meta["uncovered_nodes"])
    print(json.dumps(result, indent=2))


if __name__ == "__main__":
    main()
