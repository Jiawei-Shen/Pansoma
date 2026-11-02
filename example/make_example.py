"""Build the Pansoma minimal example (HG008-T Illumina, GRCh38 chr1:242,500,000-243,000,000, HPRC v1.1 d9 graph).

Run once by the developer; the reviewer only downloads its output directory (example/data). Steps:
  gam      every record of the genome sorted GAM that visits a graph node of the window's node-ID range (+2 kb
           flanks), written as an unsorted GAM; `vg gamsort -i` then makes the example's sorted GAM + GAI
  graph    graph-index SQLite of the nodes those records visit plus the whole node-ID range (same schema and metadata)
  refpath  GRCh38 reference-path directory cut to those node IDs (node-ID arrays end at the largest one; values
           outside the kept IDs are 'no such node'); the chr1 reference-walk arrays cut to the window
  fasta    GRCh38 chr1 (no-alt analysis set) with every base outside chr1:242.40-243.10 Mb set to N, bgzipped + faidx
  vcfs     the four PoN VCFs and the HG008-T somatic truth VCF / BED cut to the window (+100 kb), bgzip + tabix
  ckpt     the two checkpoints without optimizer / scheduler state

Usage: python make_example.py OUT_DIR STEP [STEP ...]   (developer python; imports the pipeline from the clone)
"""
import gzip
import hashlib
import json
from pathlib import Path
import shutil
import sqlite3
import subprocess
import sys

S = Path("/scratch/jshen/Github/Pansoma/tmp/assessment_20261030/r2code")
sys.path.insert(0, str(S / "repo"))

import numpy as np  # noqa: E402
import pysam  # noqa: E402

from indexed_gam_pipeline_v4.gam_reader import IndexedGam, decode  # noqa: E402
from indexed_gam_pipeline_v4.tensor_postprocessing.reference_path import ReferencePath, ARRAYS, PATH_ARRAYS  # noqa: E402

GAM = Path("/scratch/jshen/data/HG008_GIAB/vg_giraffe_GAM/Liss_lab_BCM_Illumina-WGS_20240313/"
           "HG008-T_p23_BCM_IlluminaWGS_20240313.sorted.gam")
G = Path("/scratch/jshen/data/pansoma_v2_tensors/graph_index")
GRAPH = G / "hprc-v1.1-mc-grch38.d9.graph_index.sqlite"
REFPATH = G / "hprc-v1.1-mc-grch38.d9.grch38_path"
CHR_INDEX = G / "hprc-v1.1-mc-grch38.d9.chr_node_ranges.tsv"
FASTA = Path("/scratch/jshen/data/HapMap/GCA_000001405.15_GRCh38_no_alt_analysis_set.fasta")
PON = Path("/scratch/jshen/data/Pansoma/panel_of_normal_VCFs")
PONS = ["af-only-gnomad.hg38.vcf.gz", "Homo_sapiens_assembly38.dbsnp138.vcf.gz", "1000g_pon.hg38.vcf.gz",
        "CoLoRSdb.GRCh38.v1.1.0.deepvariant.glnexus.vcf.gz"]
TRUTH = Path("/scratch/jshen/data/HG008_GIAB/draft_v02_benchmark")
TRUTH_VCF = TRUTH / "HG008-T_somatic_smvar_benchmark_v0.2_tumorvariants.vcf.gz"
TRUTH_BED = TRUTH / "HG008-T_somatic_smvar_benchmark_v0.2_all.bed"
RIS = Path("/scratch/jshen/data/pansoma_net_v2_runs/ris_7fc7396_20261012")
CHECKPOINTS = {"SNV": RIS / "HG008_Illumina_SNV_d9_normbq" / "best.pth",
               "INDEL": RIS / "HG008_Illumina_INDEL_d9_normbq" / "best.pth"}
VG = "/scratch/jshen/bin/vg_v1.77.0"

CHROM = "chr1"
WINDOW = (242_500_000, 243_000_000)          # the evaluation window
READS = (WINDOW[0] - 2_000, WINDOW[1] + 2_000)  # reads of every node starting in here
FASTA_KEEP = (242_400_000, 243_100_000)      # FASTA bases kept (the rest N)
VCF_KEEP = (242_400_000, 243_100_000)        # PoN / truth records kept
NAME = "HG008T_Illumina.chr1_242.5-243.0Mb"


def varint(value):
    out = bytearray()
    while value > 127:
        out.append((value & 127) | 128)
        value >>= 7
    out.append(value)
    return bytes(out)


def node_range():
    path = ReferencePath(REFPATH)
    part = path.contig_slice(CHROM)
    nodes, starts = np.asarray(path.path_nodes[part]), np.asarray(path.path_starts[part])
    keep = (starts >= READS[0]) & (starts < READS[1])
    return int(nodes[keep].min()), int(nodes[keep].max())


def step_gam(out):
    lo, hi = node_range()
    reader = IndexedGam(GAM, cache_bytes=1 << 30)
    raws, groups, counts, metrics = [], [], {}, {}
    from collections import defaultdict
    counts = defaultdict(int)
    reader._collect(range(lo, hi + 1), metrics, raws, groups, counts)
    visited = set()
    for raw in raws:
        visited.update(m.position.node_id for m in decode(raw).path.mapping)
    unsorted = out / "work" / f"{NAME}.unsorted.gam"
    unsorted.parent.mkdir(parents=True, exist_ok=True)
    with pysam.BGZFile(str(unsorted), "wb") as stream:  # vg's framing: groups of <= 1000 records, tag "GAM"
        for i in range(0, len(raws), 1000):
            chunk = raws[i:i + 1000]
            stream.write(varint(len(chunk) + 1) + varint(3) + b"GAM" +
                         b"".join(varint(len(r)) + r for r in chunk))
            stream.flush()
    (out / "work" / "visited_nodes.txt").write_text("\n".join(map(str, sorted(visited))) + "\n")
    gam_dir = out / "data" / "gam"
    gam_dir.mkdir(parents=True, exist_ok=True)
    sorted_gam = gam_dir / f"{NAME}.sorted.gam"
    with open(sorted_gam, "wb") as stream:
        subprocess.run([VG, "gamsort", "-t", "4", "-i", str(sorted_gam) + ".gai", str(unsorted)], stdout=stream,
                       check=True)
    stats = {}
    for label, path in (("unsorted", unsorted), ("sorted", sorted_gam)):
        stats[label] = subprocess.run([VG, "stats", "-a", str(path)], capture_output=True, text=True,
                                      check=True).stdout
    info = dict(source_gam=str(GAM), node_range=[lo, hi], window=dict(chrom=CHROM, reads=READS),
                records=len(raws), visited_nodes=len(visited), min_visited=min(visited), max_visited=max(visited),
                vg=subprocess.run([VG, "version"], capture_output=True, text=True).stdout.splitlines()[0],
                vg_stats=stats, fetch_metrics=metrics)
    (out / "work" / "gam_slice.json").write_text(json.dumps(info, indent=2) + "\n")
    print(json.dumps({k: v for k, v in info.items() if k != "vg_stats"}, indent=2))
    print(stats["sorted"])


def kept_nodes(out):
    lo, hi = json.loads((out / "work" / "gam_slice.json").read_text())["node_range"]
    visited = {int(x) for x in (out / "work" / "visited_nodes.txt").read_text().split()}
    return sorted(visited | set(range(lo, hi + 1)))


def step_graph(out):
    nodes = kept_nodes(out)
    target = out / "data" / "graph" / "hprc-v1.1-mc-grch38.d9.chr1_242.5Mb.graph_index.sqlite"
    target.parent.mkdir(parents=True, exist_ok=True)
    if target.exists():
        target.unlink()
    src = sqlite3.connect(GRAPH.as_uri() + "?mode=ro", uri=True)
    (meta,), = src.execute("SELECT value FROM graph_metadata").fetchall()
    meta = json.loads(meta)
    rows = []
    for i in range(0, len(nodes), 900):
        chunk = nodes[i:i + 900]
        rows += src.execute("SELECT node_id, seq, distinct_path_count FROM nodes WHERE node_id IN (" +
                            ",".join("?" * len(chunk)) + ")", chunk).fetchall()
    src.close()
    meta["example_slice"] = dict(source=str(GRAPH), nodes_in_source=meta.get("nodes"), nodes=len(rows),
                                 min_node=nodes[0], max_node=nodes[-1],
                                 note="nodes visited by the example GAM's records plus its node-ID range")
    with sqlite3.connect(target) as db:
        db.execute("CREATE TABLE graph_metadata(value TEXT)")
        db.execute("INSERT INTO graph_metadata VALUES(?)", (json.dumps(meta),))
        db.execute("CREATE TABLE nodes(node_id INTEGER PRIMARY KEY, seq TEXT, distinct_path_count INTEGER)")
        db.executemany("INSERT INTO nodes VALUES(?,?,?)", rows)
    print(f"graph index slice: {len(rows)} of {len(nodes)} wanted nodes, {target.stat().st_size} bytes")
    if len(rows) != len(nodes):  # IDs absent from the (frequency-filtered) graph
        print(f"  {len(nodes) - len(rows)} node IDs of the range are not in the graph")


def step_refpath(out):
    nodes = np.asarray(kept_nodes(out), dtype=np.int64)
    path = ReferencePath(REFPATH)
    target = out / "data" / "graph" / "hprc-v1.1-mc-grch38.d9.chr1_242.5Mb.grch38_path"
    target.mkdir(parents=True, exist_ok=True)
    size = int(nodes.max()) + 1
    fill = dict(lengths=0, chrom=-1, start0=0, reverse=False, visits=0)
    for name in ARRAYS:
        source = getattr(path, name)
        values = np.full(size, fill[name], dtype=source.dtype)
        values[nodes] = np.asarray(source[nodes])
        np.save(target / f"{name}.npy", values)
    part = path.contig_slice(CHROM)
    starts = np.asarray(path.path_starts[part])
    keep = np.flatnonzero((starts >= READS[0]) & (starts < READS[1]))
    first, last = int(keep[0]), int(keep[-1]) + 1
    for name in PATH_ARRAYS:
        np.save(target / f"{name}.npy", np.asarray(getattr(path, name)[part][first:last]))
    meta = json.loads((REFPATH / "meta.json").read_text())
    contig = next(c for c in meta["contigs"] if c["name"] == CHROM)
    sub = np.asarray(path.path_nodes[part][first:last])
    contig = dict(contig, offset=0, nodes=int(last - first), min_node=int(sub.min()), max_node=int(sub.max()))
    # chrom.npy indexes meta["contigs"]: chr1 is index 0 in the source too (asserted), so its values stay valid
    assert meta["contigs"][0]["name"] == CHROM
    meta.update(contigs=[contig], max_node=size - 1,
                example_slice=dict(source=str(REFPATH), node_ids=[int(nodes.min()), int(nodes.max())],
                                   kept_node_ids=int(len(nodes)), walk=dict(contig=CHROM, start=READS[0], end=READS[1]),
                                   note="node-ID arrays end at max_node; IDs not kept read as 'no such node'; "
                                        "path arrays hold the chr1 GRCh38 walk inside the window only"))
    (target / "meta.json").write_text(json.dumps(meta, indent=2) + "\n")
    with open(REFPATH / "walks.ndjson") as src, open(target / "walks.ndjson", "w") as dst:
        for line in src:
            w = json.loads(line)
            if w["sample"] == meta["reference_sample"] and w["contig"] == CHROM:
                dst.write(line)
    print(f"reference path slice: {size} node-ID entries, {last - first} walk entries")


def step_fasta(out):
    target = out / "data" / "reference" / "GRCh38_no_alt.chr1_masked_outside_242.4-243.1Mb.fa"
    target.parent.mkdir(parents=True, exist_ok=True)
    with pysam.FastaFile(str(FASTA)) as fa:
        seq = fa.fetch(CHROM)
    masked = "N" * FASTA_KEEP[0] + seq[FASTA_KEEP[0]:FASTA_KEEP[1]] + "N" * (len(seq) - FASTA_KEEP[1])
    assert len(masked) == len(seq)
    with open(target, "w") as stream:
        stream.write(f">{CHROM}\n")
        for i in range(0, len(masked), 60):
            stream.write(masked[i:i + 60] + "\n")
    pysam.tabix_compress(str(target), str(target) + ".gz", force=True)
    target.unlink()
    pysam.faidx(str(target) + ".gz")
    print(f"FASTA: {target}.gz {Path(str(target) + '.gz').stat().st_size} bytes")


def cut_vcf(source, target, region):
    with pysam.TabixFile(str(source)) as tbx, open(target, "w") as stream:
        for line in tbx.header:
            stream.write(line + "\n")
        n = 0
        contig = CHROM if CHROM in tbx.contigs else CHROM.removeprefix("chr")
        for line in tbx.fetch(contig, *region):
            stream.write(line + "\n")
            n += 1
    pysam.tabix_index(str(target), preset="vcf", force=True)
    return n


def step_vcfs(out):
    pon = out / "data" / "pon"
    pon.mkdir(parents=True, exist_ok=True)
    for name in PONS:
        n = cut_vcf(PON / name, pon / name.removesuffix(".vcf.gz").__add__(".chr1_242.4-243.1Mb.vcf"), VCF_KEEP)
        print(f"PoN {name}: {n} records")
    truth = out / "data" / "truth"
    truth.mkdir(parents=True, exist_ok=True)
    n = cut_vcf(TRUTH_VCF, truth / "HG008-T_somatic_smvar_benchmark_v0.2_tumorvariants.chr1_242.4-243.1Mb.vcf", VCF_KEEP)
    print(f"truth: {n} records")
    bed = truth / "HG008-T_somatic_smvar_benchmark_v0.2_all.chr1_242.5-243.0Mb.bed"
    rows = []
    for line in TRUTH_BED.read_text().splitlines():
        c, s, e = line.split("\t")[:3]
        s, e = int(s), int(e)
        if c == CHROM and e > WINDOW[0] and s < WINDOW[1]:
            rows.append(f"{c}\t{max(s, WINDOW[0])}\t{min(e, WINDOW[1])}")
    bed.write_text("\n".join(rows) + "\n")
    print(f"benchmark BED: {len(rows)} intervals, {sum(int(r.split()[2]) - int(r.split()[1]) for r in rows)} bp")


def step_ckpt(out):
    import torch
    target = out / "data" / "models"
    target.mkdir(parents=True, exist_ok=True)
    for kind, src in CHECKPOINTS.items():
        c = torch.load(src, map_location="cpu", weights_only=False)
        dropped = [k for k in ("optimizer_state_dict", "scheduler_state_dict") if k in c]
        for k in dropped:
            del c[k]
        c["release_note"] = dict(source=src.parent.name, source_sha256=sha256(src), dropped=dropped)
        dst = target / f"{src.parent.name}.pth"
        torch.save(c, dst)
        print(f"{kind}: {dst} {dst.stat().st_size} bytes (dropped {dropped}; epoch {c['epoch']}, "
              f"somatic_threshold {c.get('somatic_threshold')})")


def sha256(path):
    h = hashlib.sha256()
    with open(path, "rb") as stream:
        for block in iter(lambda: stream.read(1 << 22), b""):
            h.update(block)
    return h.hexdigest()


def step_chrindex(out):
    target = out / "data" / "graph"
    target.mkdir(parents=True, exist_ok=True)
    for src in (CHR_INDEX, CHR_INDEX.with_suffix(".json")):  # .tsv and its metadata (tsv_sha256)
        shutil.copy2(src, target / src.name)
    print("chr index copied")


STEPS = dict(gam=step_gam, graph=step_graph, refpath=step_refpath, fasta=step_fasta, vcfs=step_vcfs,
             ckpt=step_ckpt, chrindex=step_chrindex)

if __name__ == "__main__":
    out = Path(sys.argv[1])
    for step in sys.argv[2:]:
        print(f"== {step}", flush=True)
        STEPS[step](out)
