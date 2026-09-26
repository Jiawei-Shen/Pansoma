"""Read-level view of somatic truth alleles: how do the ALT reads represent the variant (v6 tensors)?

For each locus (every no_candidate miss plus controls that did get a tensor):
  * fetch every GAM record touching the GRCh38 nodes of [span - A, span + A) (A = 30 bp),
  * decode it with the FROZEN v6 decoder (v6_run/source, the code that built the tensors), restricted to the
    mappings around the window plus >= TRIM_BP (600) read bases on each side (v6 normalizes indels across
    nodes; validate_trim.py checks this against whole-record decoding),
  * keep records that are anchored on GRCh38 on both sides (a matched base within 10 bp outside
    the window), rebuild the read's own bases between the anchors in GRCh38 orientation and
    compare with the REF and the ALT haplotype (edit distance),
  * for ALT-like records (closer to ALT than to REF) record how the alignment represents it:
    the v6 candidates it produced there, whether one is a truth key, whether the path leaves
    GRCh38 (alternative graph branch), complex / over-length edits.
Output: reads.json (one record per locus) and reads_summary.json.
"""
import csv, importlib.util, json, os, random, sys
from collections import Counter, defaultdict
from bisect import bisect_left, bisect_right
from concurrent.futures import ProcessPoolExecutor
from pathlib import Path

import numpy as np
import pysam

sys.path.insert(0, "/scratch/jshen/data/pansoma_v2_tensors/Liss_lab_PacBio_Revio_20240125/v6_run/source")  # the run's frozen package (indexed_gam_pipeline_v2 was removed from the repo)
from indexed_gam_pipeline_v2 import vg_pb2
from indexed_gam_pipeline_v2.gam_reader import IndexedGam
from indexed_gam_pipeline_v2.graph_index import GraphIndex
from indexed_gam_pipeline_v2.tensor_postprocessing.reference_path import ReferencePath, rc
from indexed_gam_pipeline_v2.tensor_postprocessing.truth_labels import placements

csv.field_size_limit(sys.maxsize)
RUN = Path("/scratch/jshen/data/pansoma_v2_tensors/Liss_lab_PacBio_Revio_20240125")
GAM = "/scratch/jshen/data/HG008_GIAB/raw_sequencing_data/Liss_lab_PacBio_Revio_20240125/HG008-T_PacBio-HiFi-Revio_20240125_116x.AF-HPRC.sorted.gam"
GRAPH = "/scratch/jshen/data/pansoma_v2_tensors/graph_index/hprc-v1.1-mc-grch38.d9.graph_index.sqlite"
REFPATH = "/scratch/jshen/data/pansoma_v2_tensors/graph_index/hprc-v1.1-mc-grch38.d9.grch38_path"
FASTA = "/scratch/jshen/data/HapMap/GCA_000001405.15_GRCh38_no_alt_analysis_set.fasta"
TRUTH = "/scratch/jshen/data/pansoma_v2_tensors/truth/somatic.graph.tsv"
OUT = Path(__file__).resolve().parent
A, SLACK, MIN_MAPQ, FAR, FETCH = 30, 10, 10, 2000, 300
TRIM_BP = int(os.environ.get("TRIM_BP", 600))  # 0 = decode whole records

spec = importlib.util.spec_from_file_location("v6_candidates", RUN / "v6_run/source/indexed_gam_pipeline_v2/candidates.py")
v6 = importlib.util.module_from_spec(spec)
spec.loader.exec_module(v6)


def levenshtein(a, b):
    """Edit distance, one numpy row at a time (the left dependency via a running minimum)."""
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


def trim(alignment, nodes, margin=TRIM_BP):
    """The alignment restricted to its mappings on `nodes` plus flanking mappings holding >= `margin` read bases
    on each side, with the matching slice of sequence and quality. v6 left-normalizes indels across nodes, so the
    margin must exceed any repeat an indel could slide through (checked against whole-record decoding)."""
    maps = alignment.path.mapping
    hits = [i for i, m in enumerate(maps) if m.position.node_id in nodes]
    if not hits:
        return None
    lengths = [sum(e.to_length for e in m.edit) for m in maps]
    i0, got = hits[0], 0
    while i0 > 0 and got < margin:
        i0 -= 1
        got += lengths[i0]
    i1, got = hits[-1] + 1, 0
    while i1 < len(maps) and got < margin:
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


def loci():
    status = {}
    with open(RUN / "v6_tensors/somatic.recall.tsv") as f:
        for row in csv.DictReader(f, delimiter="\t"):
            status[row["truth_id"]] = row["status"]
    rows = []
    with open(TRUTH) as f:
        for row in csv.DictReader(f, delimiter="\t"):
            row.update(status=status[row["truth_id"]], pos0=int(row["pos0"]),
                       keys=row["keys"].split(",") if row["keys"] else [])
            rows.append(row)
    # INDEL only: the SNV tensors (and SNV misses) of v6 equal v5's, whose SNV analysis is copied.
    misses = [r for r in rows if r["status"] == "no_candidate" and r["kind"] != "SNP"]
    rng = random.Random(1)
    hits = [r for r in rows if r["status"] == "tensor_representative"]
    controls = (rng.sample([r for r in hits if r["kind"] != "SNP" and max(len(r["ref"]), len(r["alt"])) >= 2], 150)
                + rng.sample([r for r in hits if r["kind"] != "SNP" and max(len(r["ref"]), len(r["alt"])) == 1], 50))
    return misses + controls


def analyse(chunk):
    path = ReferencePath(REFPATH)
    fasta = pysam.FastaFile(FASTA)
    reader = IndexedGam(GAM, cache_bytes=512 << 20)
    targets = TARGETS
    visits, chrom_of, start0, rev_of, lengths = (np.asarray(path.visits), np.asarray(path.chrom),
                                                 np.asarray(path.start0), np.asarray(path.reverse), np.asarray(path.lengths))
    results = []
    genome_chrom = None
    with GraphIndex(GRAPH) as graph:
        for m in chunk:
            if m["chrom"] != genome_chrom:
                genome_chrom, genome = m["chrom"], fasta.fetch(m["chrom"]).upper()
                contig = path.contig_index[m["chrom"]]
                part = path.contig_slice(m["chrom"])
                wnodes, wstarts = np.asarray(path.path_nodes[part]), np.asarray(path.path_starts[part])
            places = placements(genome, m["pos0"], m["ref"], m["alt"], m["kind"])
            lo = min(p[0] for p in places)
            hi = max(p[0] + len(p[1]) for p in places)
            left, right = lo - A, hi + A          # anchors: a matched base at <= left and >= right
            a = max(0, int(np.searchsorted(wstarts, left - FETCH, side="right")) - 1)
            b = int(np.searchsorted(wstarts, right + FETCH, side="right"))
            fetch_nodes = {int(n) for n in wnodes[a:b]}
            key_nodes = sorted({int(k.split(":")[0]) for k in m["keys"]})
            ref_hap = None  # built per anchor pair
            info = dict(truth_id=m["truth_id"], status=m["status"], kind=m["kind"], chrom=m["chrom"],
                        vcf_pos=int(m["vcf_pos"]), vcf_ref=m["vcf_ref"], vcf_alt=m["vcf_alt"], gt=m["gt"],
                        length=max(len(m["ref"]), len(m["alt"])), placements=len(places), span=[lo, hi],
                        key_nodes_target=[bool(n in targets) for n in key_nodes],
                        window_nodes=len(fetch_nodes), window_nodes_target=sum(n in targets for n in fetch_nodes),
                        records=0, low_mapq=0, not_spanning=0, spanning=0, ref_exact=0, alt_exact=0,
                        ref_like=0, alt_like=0, unclear=0, alt_representation=Counter(), alt_candidates=Counter(), anchor_distance=[],
                        alt_examples=[], context=genome[max(0, lo - 10):hi + 10])
            alignments = list(reader.fetch(fetch_nodes, {}))
            info["records"] = len(alignments)
            good = [x for x in alignments if x.mapping_quality > MIN_MAPQ]
            info["low_mapq"] = len(alignments) - len(good)
            if TRIM_BP:
                good = [t for t in (trim(x, fetch_nodes) for x in good) if t is not None]
            need = {mp.position.node_id for x in good for mp in x.path.mapping}
            sequences = {n: r["sequence"] for n, r in graph.get_nodes(need).items()} if need else {}
            keys = set(m["keys"])
            for x in good:
                read, unsupported = v6.decode_alignment(x, sequences, 50, None)
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
                    info["not_spanning"] += 1
                    continue
                pl, pr = keys_sorted[k], keys_sorted[j]
                (il, flip_l), (ir, flip_r) = lin_col[pl], lin_col[pr]
                if flip_l != flip_r:
                    info["not_spanning"] += 1
                    continue
                i0, i1 = (il, ir) if il < ir else (ir, il)
                between = read.columns[i0 + 1:i1]
                seq = "".join(c.read for c in between if c.read != "-")
                if flip_l:
                    seq = rc(seq)
                ref_hap = genome[pl + 1:pr]
                alt_hap = genome[pl + 1:m["pos0"]] + m["alt"] + genome[m["pos0"] + len(m["ref"]):pr]
                info["spanning"] += 1
                info["anchor_distance"].append(max(left - pl, pr - right))
                if seq == ref_hap:
                    info["ref_exact"] += 1
                    continue
                exact = seq == alt_hap
                dr, da = (1, 0) if exact else (levenshtein(seq, ref_hap), levenshtein(seq, alt_hap))
                if not exact and da >= dr:
                    info["ref_like" if dr < da else "unclear"] += 1
                    continue
                info["alt_exact" if exact else "alt_like"] += 1
                # How does this ALT record represent the variant?
                # Walk the columns in GRCh38 orientation; the site segment runs from the last GRCh38 base
                # before the variant span to the first GRCh38 base after it. The read "expresses the site on
                # GRCh38" iff that segment holds only GRCh38 bases with consecutive positions (edits allowed).
                ordered = between[::-1] if flip_l else between
                lins = []
                for c in ordered:
                    on = visits[c.node] == 1 and chrom_of[c.node] == contig
                    if not on:
                        lins.append(None)
                    elif c.boundary:
                        lins.append(("b", int(start0[c.node] + (lengths[c.node] - c.pos if rev_of[c.node] else c.pos))))
                    else:
                        lins.append(int(start0[c.node] + (lengths[c.node] - 1 - c.pos if rev_of[c.node] else c.pos)))
                before = max((i for i, x in enumerate(lins) if isinstance(x, int) and x < lo), default=-1)
                after = min((i for i, x in enumerate(lins) if isinstance(x, int) and x >= hi and i > before), default=len(lins))
                segment = lins[before + 1:after]
                bases = [x for x in segment if isinstance(x, int)]
                offpath = any(x is None for x in segment)
                expected_start = lins[before] + 1 if before >= 0 else lo
                contiguous = bases == list(range(expected_start, expected_start + len(bases)))
                reaches = (after < len(lins) and isinstance(lins[after], int)
                           and lins[after] == expected_start + len(bases))
                site_on_path = not offpath and contiguous and reaches
                local_edit = any(ordered[i].op in ("I", "D", "X", "C") for i in range(before + 1, after)) or any(
                    ordered[i].op in ("I", "D", "X", "C") and isinstance(x, int) and lo - 5 <= x <= hi + 5
                    for i, x in enumerate(lins))
                maps = {c.visit for c in between}
                obs = [o for o in read.observations if o.visit in maps]
                local_obs = []
                for o in obs:
                    lin = path.linear(o.candidate.node, o.candidate.start, o.candidate.ref, o.candidate.alt, o.candidate.kind)
                    if lin is None or (lin["chrom"] == m["chrom"] and lo - 5 <= lin["pos0"] <= hi + 5):
                        local_obs.append(o)
                complex_edit = any(e.get("reason") == "complex_replacement_not_supported" and e["mapping_index"] in maps
                                   for e in unsupported)
                too_long = any(e.get("reason", "").startswith("indel_exceeds_limit") and e["mapping_index"] in maps
                               for e in unsupported)
                key_hit = [o for o in obs if o.candidate.metadata()["candidate_id"] in keys]
                if key_hit:
                    kind = "key_observed_on_target_node" if any(o.candidate.node in targets for o in key_hit) \
                        else "key_observed_on_non_target_node"
                    if kind.endswith("_target_node") and max(o.quality for o in key_hit) < 10:
                        kind = "key_observed_low_base_quality"
                elif too_long:
                    kind = "indel_longer_than_50"
                elif complex_edit:
                    kind = "complex_replacement_edit"
                elif site_on_path:
                    kind = "site_on_grch38:" + ("other_edit" if local_obs or local_edit else "no_local_edit")
                else:
                    kind = ("site_bypassed:" + ("branch" if offpath else "skip_edge") + ":"
                            + ("edit_on_path_taken" if local_obs or local_edit else "no_edit"))
                info["alt_representation"][kind + (":exact_read" if exact else ":read_with_errors")] += 1
                for o in local_obs:
                    info["alt_candidates"][o.candidate.metadata()["candidate_id"]] += 1
                if len(info["alt_examples"]) < 3:
                    info["alt_examples"].append(dict(read=read.name, exact=exact, kind=kind,
                                                     local=sorted({o.candidate.metadata()["candidate_id"] for o in local_obs}),
                                                     ops=dict(Counter(c.op for c in between)), offpath=offpath))
            info["alt_representation"] = dict(info["alt_representation"])
            d = sorted(info.pop("anchor_distance"))
            info["anchor_distance_median"] = d[len(d) // 2] if d else None
            info["alt_candidates"] = {c: dict(reads=n, tensor_label=LABEL_OF.get(c))
                                      for c, n in info["alt_candidates"].most_common(8)}
            results.append(info)
    return results


class NodeSet:
    """Membership in a sorted int64 array: ~8 bytes per node and read-only pages that forked workers share,
    instead of a Python set (~80 bytes per node, and every worker ends up with its own copy)."""

    def __init__(self, nodes):
        self.nodes = np.unique(np.asarray(nodes, dtype=np.int64))

    def __contains__(self, node):
        i = int(np.searchsorted(self.nodes, node))
        return i < len(self.nodes) and int(self.nodes[i]) == node

    def __len__(self):
        return len(self.nodes)


def init(targets, labels):
    global TARGETS, LABEL_OF
    TARGETS, LABEL_OF = targets, labels


def tensor_labels():
    """{candidate_id of every allele of every labelled tensor: (label name, reason, representative?)}."""
    import re
    result = {}
    ids = re.compile(r'"candidate_id": "([^"]+)"')
    for kind in ("SNV", "INDEL"):
        for f in sorted((RUN / "v6_tensors" / kind).glob("chr*_labels.ndjson")):
            summary = f.with_name(f.name.replace("_labels.ndjson", "_variant_summary.ndjson"))
            with f.open() as lab, summary.open() as summ:
                for a, b in zip(lab, summ):
                    r = json.loads(a)
                    alleles = ids.findall(b.split('"alleles": ', 1)[1]) if '"alleles": ' in b else [r["candidate_id"]]
                    for c in alleles:
                        result[c] = (r["label_name"], r["reason"], c == r["candidate_id"])
    return result


def main():
    work = loci()
    if os.environ.get("LIMIT"):  # quick test: a few misses of each kind and a few controls
        rng = random.Random(2)
        work = (rng.sample([m for m in work if m["status"] == "no_candidate"], int(os.environ["LIMIT"]))
                + rng.sample([m for m in work if m["status"] != "no_candidate"], 3))
    output = OUT / os.environ.get("OUTPUT", "reads.jsonl")
    done = set()
    for f in [output] + sorted(OUT.glob("reads*.jsonl")):
        if f.exists() and (f == output or not os.environ.get("LIMIT")):
            with f.open() as stream:
                done |= {json.loads(line)["truth_id"] for line in stream if line.strip()}
    if os.environ.get("PART"):  # "k/N": this job takes every N-th locus of the full sorted list, starting at k
        k, n_parts = map(int, os.environ["PART"].split("/"))
        work = sorted(work, key=lambda m: (m["chrom"], m["pos0"]))[k::n_parts]
    work = [m for m in work if m["truth_id"] not in done]
    print(f"{len(done)} loci already done, {len(work)} to go", flush=True)
    # v6 targets: the autosome node parts plus the supplement rounds (nodes that received normalized indels)
    targets = NodeSet(np.concatenate([np.fromfile(f, sep=" ", dtype=np.int64)
                               for f in sorted((RUN / "v6_run/parts").glob("nodes_*.txt")) + sorted((RUN / "v6_run/parts").glob("supplement_0[0-9]/nodes_*.txt"))]))
    print(f"{len(targets)} target nodes", flush=True)
    work.sort(key=lambda m: (m["chrom"], m["pos0"]))
    n = int(os.environ.get("SLURM_CPUS_PER_TASK", 8))
    chunks = [work[i:i + 5] for i in range(0, len(work), 5)]
    labels = {}  # tensor labels are looked up by report.py; the workers do not need them
    finished = 0
    with ProcessPoolExecutor(n, initializer=init, initargs=(targets, labels)) as pool, output.open("a") as out:
        for part in pool.map(analyse, chunks):
            for r in part:
                out.write(json.dumps(r) + "\n")
            out.flush()
            finished += len(part)
            if finished % 100 < 5:
                print(f"{finished}/{len(work)} loci", flush=True)


if __name__ == "__main__":
    main()
