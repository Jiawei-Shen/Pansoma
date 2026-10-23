"""E20 step 4: how the COLO829BL (normal) and COLO829T (tumor) ONT reads align at the COLO829T absorbed truth alleles.

Reviewer 2, point 3 (R2.3): for false negatives at existing graph nodes, do the normal sample's reads preferentially
align to the alternative path? Per locus (loci_<graph>.tsv of loci.py) and per read, two independent labels:

  * sequence: the read's own bases between two GRCh38 anchors, against GRCh38 (REF hap) and GRCh38 + truth (ALT hap):
    ref_exact / alt_exact / ref_like / alt_like (closer by edit distance) / unclear (equally close);
  * path: the graph bases of the nodes the alignment walks between the same anchors:
    grch38    the site segment stays on GRCh38 nodes with consecutive GRCh38 positions (edits allowed);
    alt       it leaves GRCh38 and the walked graph sequence equals the ALT hap (the path that spells the truth);
    alt_like  it leaves GRCh38 and the walked sequence is closer to the ALT hap than to GRCh38 (e.g. the truth plus a
              nearby germline allele, or another length of the same repeat closer to the ALT);
    other     it leaves GRCh38 otherwise (another allele, or a branch that spells GRCh38).
A read on the ALT path whose own bases are REF-like makes an edit there: that is the misalignment the reviewer asks
about. Read selection, anchors and decoding are those of E5c (reads.py): records touching the GRCh38 nodes of
[span - 30 - 300, span + 30 + 300), MAPQ > 10, matched GRCh38 bases within 2 kb outside [span - 30, span + 30) on
both sides in one orientation; decoded with the pipeline's frozen decoder (COLO829T_ONT_d46 run source), trimmed to
the window plus 600 read bases.
Modes: STREAM=1 (default for normal / tumor) reads the GAM once in file order and hands each record to every locus
window it touches (the same membership rule as the GAI fetch); the GAI fetch (STREAM=0, and tumorgam) re-reads the large
bins of long reads for every locus and was too slow on these GAMs (2026-10-21: ~3 GB read per worker before the first
5 loci).

Inputs (environment): SAMPLE = normal | tumor | tumorgam, GRAPH = d9 | d46; SUBSET = n (a seeded random n loci).
  normal / tumor: $D/gam/<SAMPLE>_<GRAPH>.sorted.gam (extract.sh + giraffe.sh: reads of the GRCh38 BAMs at the loci,
  realigned with the GAMs' giraffe settings); tumorgam: the whole-genome COLO829T ONT GAM Pansoma used on that graph.
Output: $D/classify/<SAMPLE>_<GRAPH>.jsonl, one record per locus.
usage: SAMPLE=normal GRAPH=d9 /wanglab/jshen/anaconda3/bin/python classify.py
"""
import csv
import gzip
import json
import os
import sys
from bisect import bisect_left, bisect_right
from collections import Counter
from concurrent.futures import ProcessPoolExecutor
from pathlib import Path

import numpy as np
import pysam

SOURCE = Path("/scratch/jshen/data/pansoma_v2_tensors/COLO829T_ONT_d46/run/source")
sys.path.insert(0, str(SOURCE))
from indexed_gam_pipeline_v4 import candidates as dec  # noqa: E402
from indexed_gam_pipeline_v4.gam_reader import IndexedGam, decode, group  # noqa: E402
from indexed_gam_pipeline_v4.graph_index import GraphIndex  # noqa: E402
from indexed_gam_pipeline_v4.tensor_postprocessing.reference_path import ReferencePath, rc  # noqa: E402
from indexed_gam_pipeline_v4.tensor_postprocessing.truth_labels import placements  # noqa: E402

D = Path(os.environ.get("TEST_DIR", "/scratch/jshen/data/reviewer_experiments_20261013/E20_normal_reads_colo829bl"))
FASTA = "/scratch/jshen/data/HapMap/GCA_000001405.15_GRCh38_no_alt_analysis_set.fasta"
GRAPHS = {
    "d9": ("/scratch/jshen/data/pansoma_v2_tensors/graph_index/hprc-v1.1-mc-grch38.d9.graph_index.sqlite",
           "/scratch/jshen/data/pansoma_v2_tensors/graph_index/hprc-v1.1-mc-grch38.d9.grch38_path",
           "/scratch/jshen/data/COLO829T/COLO829T_ONT/COLO829T_ONT_merged.sorted.gam"),
    "d46": ("/scratch/jshen/data/AF-Filtered_VG_Indexes_HPRC_V2/graph_index/hprc-v2.1-mc-grch38.d46.graph_index.sqlite",
            "/scratch/jshen/data/AF-Filtered_VG_Indexes_HPRC_V2/graph_index/hprc-v2.1-mc-grch38.d46.grch38_path",
            "/scratch/jshen/data/COLO829T/COLO829T_ONT/vg_hprc_v2.1_d46/COLO829T_ONT_hprc-v2.1-d46_vg.sorted.gam"),
}
SAMPLE, GRAPH = os.environ.get("SAMPLE", "normal"), os.environ.get("GRAPH", "d9")
GRAPH_INDEX, REFPATH, TUMOR_GAM = GRAPHS[GRAPH]
GAM = TUMOR_GAM if SAMPLE == "tumorgam" else str(D / "gam" / f"{SAMPLE}_{GRAPH}.sorted.gam")
A, MIN_MAPQ, FAR, FETCH, TRIM_BP = 30, 10, 2000, 300, 600
csv.field_size_limit(sys.maxsize)


def levenshtein(a, b):
    if not a or not b:
        return len(a) + len(b)
    x = np.frombuffer(a.encode(), dtype=np.uint8)
    y = np.frombuffer(b.encode(), dtype=np.uint8)
    idx = np.arange(len(y) + 1)
    row = idx.copy()
    for i, ch in enumerate(x, 1):
        tmp = np.empty_like(row)
        tmp[0] = i
        tmp[1:] = np.minimum(row[1:] + 1, row[:-1] + (y != ch))
        row = np.minimum.accumulate(tmp - idx) + idx
    return int(row[-1])


def trim(alignment, nodes, vg_pb2):
    """E5c's trim: the mappings on `nodes` plus >= TRIM_BP read bases on each side."""
    maps = alignment.path.mapping
    hits = [i for i, m in enumerate(maps) if m.position.node_id in nodes]
    if not hits:
        return None
    lengths = [sum(e.to_length for e in m.edit) for m in maps]
    i0, got = hits[0], 0
    while i0 > 0 and got < TRIM_BP:
        i0 -= 1
        got += lengths[i0]
    i1, got = hits[-1] + 1, 0
    while i1 < len(maps) and got < TRIM_BP:
        got += lengths[i1]
        i1 += 1
    offsets = [0]
    for n in lengths:
        offsets.append(offsets[-1] + n)
    sub = vg_pb2.Alignment(name=alignment.name, mapping_quality=alignment.mapping_quality)
    for m in maps[i0:i1]:
        sub.path.mapping.add().CopyFrom(m)
    sub.sequence = alignment.sequence[offsets[i0]:offsets[i1]]
    if alignment.quality:
        sub.quality = alignment.quality[offsets[i0]:offsets[i1]]
    return sub


def seq_class(seq, ref_hap, alt_hap):
    if seq == ref_hap:
        return "ref_exact"
    if seq == alt_hap:
        return "alt_exact"
    dr, da = levenshtein(seq, ref_hap), levenshtein(seq, alt_hap)
    return "ref_like" if dr < da else "alt_like" if da < dr else "unclear"


class Context:
    """Per-process reference path arrays, graph index, FASTA."""

    def __init__(self):
        self.path = ReferencePath(REFPATH)
        self.visits, self.chrom_of, self.start0, self.rev_of, self.lengths = (
            np.asarray(self.path.visits), np.asarray(self.path.chrom), np.asarray(self.path.start0),
            np.asarray(self.path.reverse), np.asarray(self.path.lengths))
        self.fasta = pysam.FastaFile(FASTA)
        self.graph = GraphIndex(GRAPH_INDEX)
        self.graph.__enter__()
        from indexed_gam_pipeline_v4 import vg_pb2
        self.vg_pb2 = vg_pb2


def locus_params(m, genome, ctx):
    """Allele span, anchor window, the GRCh38 nodes to fetch, and the GRCh38 sequence around it (E5c's rule)."""
    contig = ctx.path.contig_index[m["chrom"]]
    part = ctx.path.contig_slice(m["chrom"])
    wnodes, wstarts = np.asarray(ctx.path.path_nodes[part]), np.asarray(ctx.path.path_starts[part])
    pos0 = int(m["pos0"])
    places = placements(genome, pos0, m["ref"], m["alt"], m["kind"])
    lo = min(p[0] for p in places)
    hi = max(p[0] + len(p[1]) for p in places)
    left, right = lo - A, hi + A
    a = max(0, int(np.searchsorted(wstarts, left - FETCH, side="right")) - 1)
    b = int(np.searchsorted(wstarts, right + FETCH, side="right"))
    ws = max(0, left - FAR - 10)
    return dict(truth_id=m["truth_id"], kind2=m["kind2"], chrom=m["chrom"], vcf_pos=int(m["vcf_pos"]), contig=contig,
                pos0=pos0, ref=m["ref"], alt=m["alt"], lo=lo, hi=hi, left=left, right=right,
                fetch_nodes=frozenset(int(n) for n in wnodes[a:b]), ws=ws, win=genome[ws:right + FAR + 10])


def classify_read(x, P, ctx):
    """None (does not span the allele) or (path label, sequence label, example)."""
    need = {mp.position.node_id for mp in x.path.mapping}
    sequences = {n: r["sequence"] for n, r in ctx.graph.get_nodes(need).items()}
    read, _ = dec.decode_alignment(x, sequences, 50, None)
    visits, chrom_of, start0, rev_of, lengths = ctx.visits, ctx.chrom_of, ctx.start0, ctx.rev_of, ctx.lengths
    contig, lo, hi, left, right = P["contig"], P["lo"], P["hi"], P["left"], P["right"]
    lin_col = {}
    for i, c in enumerate(read.columns):
        if c.boundary or c.op != "M" or c.read != c.ref:
            continue
        n = c.node
        if visits[n] != 1 or chrom_of[n] != contig:
            continue
        lin = int(start0[n] + (lengths[n] - 1 - c.pos if rev_of[n] else c.pos))
        lin_col.setdefault(lin, (i, c.reverse != bool(rev_of[n])))
    keys_sorted = sorted(lin_col)
    k = bisect_right(keys_sorted, left) - 1
    j = bisect_left(keys_sorted, right)
    if k < 0 or j >= len(keys_sorted) or left - keys_sorted[k] > FAR or keys_sorted[j] - right > FAR:
        return None
    pl, pr = keys_sorted[k], keys_sorted[j]
    (il, flip_l), (ir, flip_r) = lin_col[pl], lin_col[pr]
    if flip_l != flip_r:
        return None
    i0, i1 = (il, ir) if il < ir else (ir, il)
    between = read.columns[i0 + 1:i1]
    seq = "".join(c.read for c in between if c.read != "-")
    walked = "".join(c.ref for c in between if not c.boundary and c.ref != "-")
    if flip_l:
        seq, walked = rc(seq), rc(walked)
    ws, win, pos0 = P["ws"], P["win"], P["pos0"]
    ref_hap = win[pl + 1 - ws:pr - ws]
    alt_hap = win[pl + 1 - ws:pos0 - ws] + P["alt"] + win[pos0 + len(P["ref"]) - ws:pr - ws]
    # path: does the site segment stay on GRCh38 nodes with consecutive positions (E5c's site_on_path)?
    ordered = between[::-1] if flip_l else between
    lins = []
    for c in ordered:
        if not (visits[c.node] == 1 and chrom_of[c.node] == contig):
            lins.append(None)
        elif c.boundary:
            lins.append(("b", 0))
        else:
            lins.append(int(start0[c.node] + (lengths[c.node] - 1 - c.pos if rev_of[c.node] else c.pos)))
    before = max((i for i, v in enumerate(lins) if isinstance(v, int) and v < lo), default=-1)
    after = min((i for i, v in enumerate(lins) if isinstance(v, int) and v >= hi and i > before), default=len(lins))
    segment = lins[before + 1:after]
    bases = [v for v in segment if isinstance(v, int)]
    offpath = any(v is None for v in segment)
    expected = lins[before] + 1 if before >= 0 else lo
    contiguous = bases == list(range(expected, expected + len(bases)))
    reaches = after < len(lins) and isinstance(lins[after], int) and lins[after] == expected + len(bases)
    if not offpath and contiguous and reaches:
        pcls = "grch38"
    elif walked == alt_hap:
        pcls = "alt"
    else:
        dr, da = levenshtein(walked, ref_hap), levenshtein(walked, alt_hap)
        pcls = "alt_like" if da < dr else "other"
    scls = seq_class(seq, ref_hap, alt_hap)
    example = (dict(read=read.name, path=pcls, seq=scls, walked_len=len(walked), ref_len=len(ref_hap),
                    alt_len=len(alt_hap)) if pcls != "grch38" and scls.startswith("ref") else None)
    return pcls, scls, example


def new_info(P):
    return dict(truth_id=P["truth_id"], kind2=P["kind2"], chrom=P["chrom"], vcf_pos=P["vcf_pos"], span=[P["lo"], P["hi"]],
                records=0, low_mapq=0, not_spanning=0, spanning=0, table=Counter(), path=Counter(), seq=Counter(),
                examples=[])


def add(info, result):
    if result is None:
        info["not_spanning"] += 1
        return
    pcls, scls, example = result
    info["spanning"] += 1
    info["table"][f"{pcls}|{scls}"] += 1
    info["path"][pcls] += 1
    info["seq"][scls] += 1
    if example and len(info["examples"]) < 2:
        info["examples"].append(example)


def finish(info):
    for key in ("table", "path", "seq"):
        info[key] = dict(info[key])
    return info


def analyse(chunk):
    """Fetch mode: per locus, the records the GAI returns for its GRCh38 window nodes."""
    ctx = Context()
    reader = IndexedGam(GAM, cache_bytes=512 << 20)
    results, genome_chrom = [], None
    for m in chunk:
        if m["chrom"] != genome_chrom:
            genome_chrom, genome = m["chrom"], ctx.fasta.fetch(m["chrom"]).upper()
        P = locus_params(m, genome, ctx)
        info = new_info(P)
        for x in reader.fetch(P["fetch_nodes"], {}):
            info["records"] += 1
            if x.mapping_quality <= MIN_MAPQ:
                info["low_mapq"] += 1
                continue
            add(info, classify_read(trim(x, P["fetch_nodes"], ctx.vg_pb2), P, ctx))
        results.append(finish(info))
    return results


# --- stream mode: one sequential pass over a GAM that holds only reads near the loci ----------------------------
_CTX, _PARAMS, _NODE_LOCI = None, None, None


def stream_init(params, node_loci):
    global _CTX, _PARAMS, _NODE_LOCI
    _CTX, _PARAMS, _NODE_LOCI = Context(), params, node_loci


def stream_work(raws):
    """Records of one batch -> [(locus index, 'low_mapq' | classify_read result)] for every locus window they touch
    (the membership rule of IndexedGam.fetch: the record's path visits a node of the window)."""
    out = []
    for raw in raws:
        x = decode(raw)
        hit = set()
        for mp in x.path.mapping:
            hit.update(_NODE_LOCI.get(mp.position.node_id, ()))
        for i in sorted(hit):
            if x.mapping_quality <= MIN_MAPQ:
                out.append((i, "low_mapq"))
            else:
                out.append((i, classify_read(trim(x, _PARAMS[i]["fetch_nodes"], _CTX.vg_pb2), _PARAMS[i], _CTX)))
    return out


def stream(work, n, output):
    ctx = Context()
    params, genome_chrom = [], None
    for m in sorted(work, key=lambda m: (m["chrom"], int(m["pos0"]))):
        if m["chrom"] != genome_chrom:
            genome_chrom, genome = m["chrom"], ctx.fasta.fetch(m["chrom"]).upper()
        params.append(locus_params(m, genome, ctx))
    del genome
    node_loci = {}
    for i, P in enumerate(params):
        for node in P["fetch_nodes"]:
            node_loci.setdefault(node, []).append(i)
    infos = [new_info(P) for P in params]

    def batches():
        with gzip.open(GAM, "rb") as fh:
            batch = []
            while True:
                messages = group(fh)
                if messages is None:
                    break
                batch.extend(messages)
                if len(batch) >= 200:
                    yield batch
                    batch = []
            if batch:
                yield batch

    seen = 0
    with ProcessPoolExecutor(n, initializer=stream_init, initargs=(params, node_loci)) as pool:
        for part in pool.map(stream_work, batches()):
            for i, result in part:
                infos[i]["records"] += 1
                if result == "low_mapq":
                    infos[i]["low_mapq"] += 1
                else:
                    add(infos[i], result)
            seen += 1
            if seen % 100 == 0:
                print(f"{seen} batches", flush=True)
    with output.open("a") as out:
        for info in infos:
            out.write(json.dumps(finish(info)) + "\n")


def main():
    with open(D / f"loci_{GRAPH}.tsv") as f:
        work = list(csv.DictReader(f, delimiter="\t"))
    if os.environ.get("LIMIT"):
        work = work[: int(os.environ["LIMIT"])]
    if os.environ.get("SUBSET"):  # tumorgam check: a seeded random subset (the whole-genome GAM is slow to fetch)
        import random
        work = random.Random(1).sample(sorted(work, key=lambda m: int(m["truth_id"])), int(os.environ["SUBSET"]))
    out_dir = D / "classify"
    out_dir.mkdir(parents=True, exist_ok=True)
    output = out_dir / f"{SAMPLE}_{GRAPH}{os.environ.get('SUFFIX', '')}.jsonl"
    done = set()
    if output.exists():
        with output.open() as f:
            done = {json.loads(line)["truth_id"] for line in f if line.strip()}
    work = sorted((m for m in work if m["truth_id"] not in done), key=lambda m: (m["chrom"], int(m["pos0"])))
    print(f"{SAMPLE} {GRAPH} {GAM}: {len(done)} done, {len(work)} to go", flush=True)
    n = int(os.environ.get("SLURM_CPUS_PER_TASK", 4))
    if os.environ.get("STREAM", "1") == "1" and SAMPLE != "tumorgam":
        if done:
            raise SystemExit(f"stream mode writes every locus at the end; move the partial {output} away first")
        stream(work, n, output)
        return
    chunks = [work[i:i + 5] for i in range(0, len(work), 5)]
    finished = 0
    with ProcessPoolExecutor(n) as pool, output.open("a") as out:
        for part in pool.map(analyse, chunks):
            for r in part:
                out.write(json.dumps(r) + "\n")
            out.flush()
            finished += len(part)
            if finished % 100 < 5:
                print(f"{finished}/{len(work)}", flush=True)


if __name__ == "__main__":
    main()
