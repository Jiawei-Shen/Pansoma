"""One-time graph preparation: the GFA, the reference-path directory, the chromosome block table and the graph
index audit (the graph index itself: tools.graph_index_build; all steps as one job: tools/jobs/graph_prep.sh).

    python -m indexed_gam_pipeline_v4.tools.graph_prep gfa --gbz G.gbz --output G.gfa [--threads 16] [--vg VG]
    python -m indexed_gam_pipeline_v4.tools.graph_prep ref-path-scan --gfa G.gfa --output DIR [--reference-sample GRCh38]
    python -m indexed_gam_pipeline_v4.tools.graph_prep ref-path-check --path DIR --graph-index DB --fasta FA [--samples 100000]
    python -m indexed_gam_pipeline_v4.tools.graph_prep components --gbz G.gbz --reference-path DIR --output DIR [--vg VG]
    python -m indexed_gam_pipeline_v4.tools.graph_prep chr-index --components-dir D --reference-path DIR --output PREFIX [--graph-index DB]
    python -m indexed_gam_pipeline_v4.tools.graph_prep audit --graph-index DB --gfa G.gfa --chr-index TSV --output JSON

The formats are documented with their run-time readers: tensor_postprocessing.reference_path
(gfa-reference-path directory, ReferencePath) and tensor_postprocessing.chr_index (chr-node-ranges
TSV plus JSON with tsv_sha256, ChrIndex). scan writes the reference-path FORMAT, build the
CHR_INDEX_FORMAT; nothing here is imported at run time. `gfa` and `components` run vg: --vg, else $PANSOMA_VG
(scripts/use_vg.sh), else vg on PATH.
"""
import argparse
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
# S lines -> "id<TAB>length" into LENGTHS (LN:i: when the sequence is '*'); W lines pass through.
AWK = r'''BEGIN { FS = OFS = "\t" }
$1 == "S" { n = length($3)
            if ($3 == "*") { n = -1; for (i = 4; i <= NF; i++) if ($i ~ /^LN:i:/) n = substr($i, 6) + 0 }
            print $2, n > LENGTHS; next }
$1 == "W" { print }'''


def parse_walk(walk):
    """b'>12<3' -> (ids int64, reverse bool)."""
    raw = np.frombuffer(walk, dtype=np.uint8)
    reverse = raw[(raw == 60) | (raw == 62)] == 60
    ids = np.array(walk.translate(SEPARATORS).split(), dtype=np.int64)
    if ids.size != reverse.size or not ids.size:
        raise ValueError("Malformed GFA walk")
    return ids, reverse


def scan(gfa, output, reference_sample="GRCh38"):
    """Build the output directory from one streaming pass (awk prefilter) over `gfa`."""
    started = time.perf_counter()
    output = Path(output)
    if output.exists():
        raise ValueError(f"Output exists: {output}")
    work = output.with_name(output.name + ".tmp")
    work.mkdir(parents=True)
    lengths_txt = work / "segment_lengths.txt"
    reference, walks = [], []
    process = subprocess.Popen(["awk", "-v", f"LENGTHS={lengths_txt}", AWK, str(gfa)], stdout=subprocess.PIPE)
    with (work / "walks.ndjson").open("w") as log:
        for line in process.stdout:
            fields = line.rstrip(b"\n").split(b"\t")
            if len(fields) != 7:
                raise ValueError("Malformed GFA W line")
            ids, reverse = parse_walk(fields[6])
            sample, hap, contig = (f.decode() for f in fields[1:4])
            record = dict(sample=sample, hap=hap, contig=contig, start=int(fields[4]), end=int(fields[5]),
                          nodes=int(ids.size), min=int(ids.min()), max=int(ids.max()))
            log.write(json.dumps(record) + "\n")
            walks.append(record)
            if sample == reference_sample:
                reference.append((record, ids, reverse))
    process.stdout.close()
    if process.wait():
        raise RuntimeError(f"awk failed on {gfa}")
    if not reference:
        raise ValueError(f"No W lines for sample {reference_sample} in {gfa}")
    pairs = np.fromfile(lengths_txt, sep=" ", dtype=np.int64).reshape(-1, 2)
    if (pairs[:, 1] < 0).any() or (pairs[:, 0] <= 0).any():
        raise ValueError("Segments need a positive integer ID and a sequence or LN tag")
    size = int(pairs[:, 0].max()) + 1
    lengths = np.zeros(size, dtype=np.int32)
    lengths[pairs[:, 0]] = pairs[:, 1]
    if np.count_nonzero(lengths) != len(pairs):
        raise ValueError("Duplicate or empty segments in the GFA")
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
                segments=int(len(pairs)), walks=len(walks), samples=len({w["sample"] for w in walks}),
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
        report["gfa_max_node"] = path.meta["max_node"]
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
                        and report["graph_index_nodes"] == report["gfa_max_node"])
    meta = read_json(Path(directory) / "meta.json")
    meta["checks"]["graph_index_and_fasta"] = report
    write_json(Path(directory) / "meta.json", meta)
    if not report["passed"]:
        raise ValueError(f"Reference path check failed: {report}")
    return report


# --- chromosome components (the chr-index input) ------------------------------------------------

def components(gbz, reference_path, output, vg=None, threads=1, autosomes=AUTOSOMES):
    """<output>/chrN/chrN.component.nodes.raw.txt (one node ID per line) for chr-index: the connected component of
    each autosome's reference path, `vg chunk -C -p <sample>#<hap>#<contig>` on the GBZ, off-reference nodes
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
        name = f"{path.meta['reference_sample']}#{contigs[chrom]['hap']}#{chrom}"
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


def component_block(path):
    """(first, last, count) of one component node list; it must be one gap-free interval."""
    ids = np.fromfile(path, sep=" ", dtype=np.int64)
    if not ids.size:
        raise ValueError(f"Empty component node list: {path}")
    first, last = int(ids.min()), int(ids.max())
    if np.unique(ids).size != ids.size:
        raise ValueError(f"Duplicate node IDs in {path}")
    if ids.size != last - first + 1:
        raise ValueError(f"Component {path} is not one contiguous node-ID interval "
                         f"({ids.size} nodes in [{first}, {last}])")
    return first, last, int(ids.size)


def build(components_dir, reference_path, output, graph_index=None, autosomes=AUTOSOMES):
    """Write <output>.tsv + <output>.json; raises on any violated invariant."""
    from ..tensor_postprocessing.reference_path import ReferencePath
    components_dir, output = Path(components_dir), Path(output)
    path = ReferencePath(reference_path)
    max_node = path.meta["max_node"]
    blocks, sources = [], {}
    for chrom in autosomes:
        source = components_dir / chrom / f"{chrom}.component.nodes.raw.txt"
        first, last, count = component_block(source)
        blocks.append(dict(chrom=chrom, first_node=first, last_node=last, nodes=count, dataset="autosome",
                           source=f"vg chunk -C connected component ({source.name})"))
        sources[chrom] = stamp(source)
    blocks.sort(key=lambda b: b["first_node"])
    for a, b in zip(blocks, blocks[1:]):
        if a["last_node"] >= b["first_node"]:
            raise ValueError(f"Overlapping autosome blocks: {a['chrom']} and {b['chrom']}")
    # Non-autosomal groups start at their reference contigs' smallest node.
    starts = {}
    for contig in path.meta["contigs"]:
        if contig["name"] in autosomes:
            continue
        group = group_of(contig["name"])
        starts[group] = min(starts.get(group, contig["min_node"]), contig["min_node"])
    autosome_starts = [b["first_node"] for b in blocks]
    ordered = sorted(starts.items(), key=lambda kv: kv[1])
    for i, (group, first) in enumerate(ordered):
        if any(b["first_node"] <= first <= b["last_node"] for b in blocks):
            raise ValueError(f"Reference contigs of {group} start inside an autosome block")
        following = [s for s in autosome_starts if s > first] + [s for _, s in ordered[i + 1:]]
        last = min(following) - 1 if following else max_node
        blocks.append(dict(chrom=group, first_node=first, last_node=last, nodes=last - first + 1,
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
    meta = dict(format=CHR_INDEX_FORMAT, max_node=max_node, covered_nodes=covered, uncovered_nodes=max_node - covered,
                autosome_components=sources, reference_path=dict(path=str(Path(reference_path).resolve()),
                                                                 source=path.meta["source"]),
                walk_check=walks)
    if graph_index:
        from ..graph_index import GraphIndex
        with GraphIndex(graph_index) as graph:
            meta["graph_index"] = dict(path=str(Path(graph_index).resolve()), nodes=graph.metadata["nodes"],
                                       gbz=graph.metadata["source"])
            if graph.metadata["nodes"] != max_node:
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

# Target S lines -> "S id seq"; every W/P line visiting targets -> one "V id" per visited target (revisits once).
AUDIT_AWK = r'''BEGIN { FS = "\t"; while ((getline id < NODES) > 0) t[id] = 1 }
$1 == "S" { if ($2 in t) print "S\t" $2 "\t" $3; next }
$1 == "W" || $1 == "P" { w = ($1 == "W") ? $7 : $3; n = split(w, a, /[<>,+-]+/); split("", seen)
                         for (i = 1; i <= n; i++) if (a[i] in t) seen[a[i]] = 1
                         for (k in seen) print "V\t" k }'''


def audit(graph_index, gfa, chr_index, output, random_nodes=24, seed=20260930):
    """Independent check of the graph index against the GFA of the same GBZ, on the first node of every chr-index
    block and `random_nodes` random index nodes: the index sequence must equal the S line, and distinct_path_count
    the number of W/P lines (paths) that visit the node. Writes the report to `output` (JSON); on any difference to
    <output>.failed instead, and raises."""
    from ..graph_index import GraphIndex
    with Path(chr_index).open() as stream:
        blocks = [int(row["first_node"]) for row in csv.DictReader(stream, delimiter="\t")]
    with GraphIndex(graph_index) as index:
        low, high = index.db.execute("SELECT min(node_id), max(node_id) FROM nodes").fetchone()
        rng = random.Random(seed)
        picks = {index.db.execute("SELECT node_id FROM nodes WHERE node_id >= ? ORDER BY node_id LIMIT 1",
                                  (rng.randint(low, high),)).fetchone()[0] for _ in range(random_nodes)}
        records = index.get_nodes(set(blocks) | picks)
    nodes = sorted(records)
    output = Path(output)
    listed = output.with_name(output.name + ".nodes.tmp")
    listed.write_text("".join(f"{n}\n" for n in nodes))
    sequences, visits = {}, dict.fromkeys(nodes, 0)
    with subprocess.Popen(["awk", "-v", f"NODES={listed}", AUDIT_AWK, str(gfa)], stdout=subprocess.PIPE,
                          text=True) as process:
        for line in process.stdout:
            kind, node, *rest = line.rstrip("\n").split("\t")
            if kind == "S":
                sequences[int(node)] = rest[0]
            else:
                visits[int(node)] += 1
    if process.returncode:
        raise RuntimeError(f"awk failed on {gfa}")
    listed.unlink()
    mismatches = [dict(node=n, index_count=records[n]["distinct_path_count"], gfa_paths=visits[n],
                       sequence_equal=records[n]["sequence"] == sequences.get(n))
                  for n in nodes if records[n]["distinct_path_count"] != visits[n]
                  or records[n]["sequence"] != sequences.get(n)]
    report = dict(passed=not mismatches, nodes={str(n): records[n]["distinct_path_count"] for n in nodes},
                  block_first_nodes=blocks, random_nodes=random_nodes, seed=seed, mismatches=mismatches,
                  graph_index=str(Path(graph_index).resolve()), gfa=stamp(gfa),
                  scope="index sequence vs GFA S line, distinct_path_count vs GFA W/P lines visiting the node, "
                        "on the chr-index block starts and random nodes")
    if mismatches:  # the report goes to <output>.failed, so `output` exists only for a passed audit
        write_json(output.with_name(output.name + ".failed"), report)
        raise ValueError(f"Graph index audit failed: {mismatches[:5]}")
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
    p.add_argument("--reference-sample", default="GRCh38")

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
    p.add_argument("--random-nodes", type=int, default=24)

    args = parser.parse_args(argv)
    if args.command == "gfa":
        result = to_gfa(args.gbz, args.output, args.vg, args.threads)
    elif args.command == "components":
        result = components(args.gbz, args.reference_path, args.output, args.vg, args.threads)
    elif args.command == "audit":
        result = audit(args.graph_index, args.gfa, args.chr_index, args.output, args.random_nodes)
        result = {k: v for k, v in result.items() if k != "nodes"}
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
