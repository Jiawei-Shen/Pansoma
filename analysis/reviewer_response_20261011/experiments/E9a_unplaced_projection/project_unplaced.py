"""E9a + E9b: place the graph-coordinate ("unplaced") calls of one run's vcf_chr1 on GRCh38 through the GRCh38 nodes
that flank their bubble, and on T2T-CHM13 directly or through CHM13 anchors.

    python project_unplaced.py --vcf-dir <run>/vcf_chr1 --run RUN --output <scratch>/<RUN>/<vcf_dir_name>

Input: <vcf-dir>/<name>.graph.vcf.gz (its header gives the tensors and graph index), <name>.linear.vcf.gz and
<name>.linear.unplaced.vcf.gz (linear_vcf: records without a GRCh38 position, in graph coordinates). Graph data
(prep_graph.py): <scratch>/graph/<graph>/edges (CSR over oriented node sides) and chm13_path (CHM13 walk coordinates);
the graph's GRCh38 reference path and graph index (node sequences); the GRCh38 FASTA of the set's labels.

Per unplaced record (event = INFO KIND, START, EVREF, EVALT, PATH on node CHROM, node-forward):
1. Anchors. A Dijkstra search (bp distance) from the end of the event's last node through its successors, and from its
   start through its predecessors, stops at the first nodes with exactly one GRCh38 visit (the anchors; it never walks
   through one) and gives up after --max-bp bases or --max-nodes nodes. The pair of a left and a right anchor on the same
   contig, both traversed in the same direction relative to GRCh38 and in order, with the smallest GRCh38 span plus walk
   length, is the bubble; the walk between them through the event node is the branch.
2. Alignment. The branch (oriented to the GRCh38 forward strand), with up to 32 anchor bases on each side, is aligned
   globally (affine gaps: match 2, mismatch -4, gap 4 + 2 per base) to the GRCh38 interval between the anchors with the
   same flanks.
3. Classes (proj_class / proj_detail):
   i   written as a GRCh38 variant. SNV: the call's base is aligned to a GRCh38 base, the ALT differs from that base, and
       the alignment has no gap within 10 bases; REF = the FASTA base (detail exact; branch_base_differs when the branch
       base itself differs from GRCh38, i.e. the call sits on a germline difference). INDEL: the event's sliding range
       on the branch (+-10 bases) is aligned without a gap; the event is moved to GRCh38 and left-aligned on the FASTA
       (linear_vcf.Genome; detail exact, or composite when the sliding range +-2 has mismatches).
   ii  a GRCh38 position or interval only: the event lies in an insertion relative to GRCh38 (in_insertion), within 10
       bases of a branch INDEL so its GRCh38 representation is ambiguous (near_branch_indel), its ALT is the GRCh38 base
       (alt_is_grch38_base: the read shows GRCh38's own base on a branch node, i.e. no GRCh38 variant at all), or the
       bubble is too large to align (too_long). Interval = the GRCh38 bases between the aligned neighbours of the
       event (+-10 bases for near_branch_indel), or between the anchors.
   iii no consistent anchor pair within the limits (no_left_anchor, no_right_anchor, inconsistent_anchors).
4. CHM13 (E9b): direct when the event's node(s) have one CHM13 visit (ReferencePath.linear on the CHM13 walks, the same
   projection linear_vcf uses for GRCh38); otherwise the CHM13 anchor interval found the same way (anchor_interval);
   otherwise none.

Outputs in --output:
    projected.tsv              one row per unplaced record (columns below)
    projected.grch38.vcf.gz    the class-i records as GRCh38 records (linear_vcf header; INFO PROJ, NODE, LEFTSHIFT)
    complete.linear.vcf.gz     linear.vcf.gz + the class-i records (same (CHROM, POS, REF, ALT): one record kept, PASS
                               first, then p_somatic; INFO MERGED), sorted; a linear_vcf output for vcfeval
    complete.linear.unplaced.vcf.gz   the class ii / iii records, in graph coordinates
    complete.linear.stats.json unplaced_pass = PASS records still without a GRCh38 variant (vcfeval reads it)
    projection.json            counts per class, filter and detail, CHM13 counts, timing
"""
import argparse
import gzip
import heapq
import json
import sys
import time
from collections import Counter
from pathlib import Path

import numpy as np

REPO = Path(__file__).resolve().parents[4]
sys.path.insert(0, str(REPO))
sys.path.insert(0, str(REPO / "machine_learning"))
from indexed_gam_pipeline_v4.graph_index import GraphIndex  # noqa: E402
from indexed_gam_pipeline_v4.tensor_postprocessing.reference_path import ReferencePath, rc  # noqa: E402
from pansoma_net_v2.linear_vcf import Genome, a1_reads  # noqa: E402
from pansoma_net_v2.vcfio import format_info, meta, parse_info, read_vcf, write_vcf  # noqa: E402

SCRATCH = Path("/scratch/jshen/data/reviewer_experiments_20261013/E9a_unplaced_projection")
GRAPHS = {  # graph index file name -> prep directory name
    "hprc-v1.1-mc-grch38.d9.graph_index.sqlite": "d9",
    "hprc-v2.1-mc-grch38.d46.graph_index.sqlite": "d46",
    "hg38_linear.giraffe.graph_index.sqlite": "linear",
}
MATCH, MISMATCH, GAP_OPEN, GAP_EXT = 2, -4, 4, 2
FLANK = 32
CONTEXT = 2  # bases beyond an INDEL's sliding range checked for mismatches (detail composite)
NEAR_GAP = 10  # class i: no alignment gap within this many bases of the SNV / the INDEL's sliding range
TSV_COLUMNS = ["run", "set", "kind", "chrom_graph_id", "node", "offset", "ref", "alt", "p_somatic", "proj_class",
               "grch38_chrom", "grch38_pos", "grch38_ref", "grch38_alt", "chm13_chrom", "chm13_pos",
               # extra columns
               "candidate_id", "block", "filter", "qual", "event_kind", "evref", "evalt", "path", "af", "alt_reads",
               "proj_detail", "near_gap", "grch38_start0", "grch38_end0", "anchor_left", "anchor_right",
               "anchor_bp_left", "anchor_bp_right", "branch_len", "grch38_span", "chm13_class", "chm13_end",
               "dup_of_placed"]
EXTRA_INFO = """##INFO=<ID=PROJ,Number=1,Type=String,Description="E9a: placed on GRCh38 through the bubble anchors (detail)">
##INFO=<ID=PROJ_NEARGAP,Number=0,Type=Flag,Description="E9a: an alignment gap within 10 bases of the event">"""


def side_rc(seq, side):
    return rc(seq) if side & 1 else seq


class Graph:
    def __init__(self, name, graph_index, grch38_path):
        d = SCRATCH / "graph" / name
        self.name = name
        self.ref = ReferencePath(grch38_path)
        self.chm13 = ReferencePath(d / "chm13_path") if (d / "chm13_path" / "meta.json").exists() else None
        self.off = np.load(d / "edges" / "succ_offsets.npy", mmap_mode="r")
        self.tgt = np.load(d / "edges" / "succ_targets.npy", mmap_mode="r")
        self.lengths = self.ref.lengths
        self.index = GraphIndex(graph_index)
        self.seq_cache = {}

    def succ(self, side):
        return [int(t) for t in self.tgt[int(self.off[side]):int(self.off[side + 1])]]

    def sequences(self, nodes):
        need = [n for n in set(nodes) if n not in self.seq_cache]
        if need:
            for n, v in self.index.get_nodes(need).items():
                self.seq_cache[int(n)] = v["sequence"].upper()
        return {n: self.seq_cache[n] for n in nodes}

    def search(self, start_side, path, max_bp, max_nodes, slack=2000, max_anchors=16):
        """Anchors reachable from the end of start_side: [(side, interior bp)], parents."""
        heap, dist, parent, anchors, expanded = [(0, start_side)], {start_side: 0}, {start_side: None}, [], 0
        best = None
        while heap:
            d, s = heapq.heappop(heap)
            if d > dist[s]:
                continue
            if best is not None and d > best + slack:
                break
            node = s >> 1
            if s != start_side and path.visits[node] == 1:
                anchors.append((s, d))
                best = d if best is None else best
                if len(anchors) >= max_anchors:
                    break
                continue
            if expanded >= max_nodes:
                break
            expanded += 1
            step = 0 if s == start_side else int(self.lengths[node])
            for t in self.succ(s):
                nd = d + step
                if nd > max_bp:
                    continue
                if t not in dist or nd < dist[t]:
                    dist[t], parent[t] = nd, s
                    heapq.heappush(heap, (nd, t))
        return anchors, parent


def chain(parent, side):
    out = []
    while side is not None:
        out.append(side)
        side = parent[side]
    return out[::-1]  # start side first


def best_pair(path, left, right, max_span):
    """The left/right anchor pair of one bubble: same contig and direction, in order; smallest span + walk."""
    best = None
    for sl, dl in left:
        ml, tl = sl >> 1, (sl ^ 1) & 1  # forward-walk orientation of the left anchor
        for sr, dr in right:
            mr, tr = sr >> 1, sr & 1
            if path.contigs[int(path.chrom[ml])] != path.contigs[int(path.chrom[mr])]:
                continue  # by name: a CHM13 contig can be several W lines
            plus_l, plus_r = tl == int(path.reverse[ml]), tr == int(path.reverse[mr])
            if plus_l != plus_r:
                continue
            if plus_l:
                g0, g1 = int(path.start0[ml]) + int(path.lengths[ml]), int(path.start0[mr])
            else:
                g0, g1 = int(path.start0[mr]) + int(path.lengths[mr]), int(path.start0[ml])
            if g1 < g0 or g1 - g0 > max_span:
                continue
            score = (g1 - g0) + dl + dr
            if best is None or score < best[0]:
                best = (score, sl, dl, sr, dr, plus_l, g0, g1)
    return best


def global_align(a, b):
    """Affine-gap global alignment; for each index of a, the aligned index of b or -1."""
    m, n = len(a), len(b)
    if m == 0:
        return np.zeros(0, dtype=np.int64)
    if n == 0:
        return np.full(m, -1, dtype=np.int64)
    if m == n:
        mism = sum(x != y for x, y in zip(a, b))
        if mism <= max(2, m // 10):
            return np.arange(m, dtype=np.int64)
    A = np.frombuffer(a.encode(), dtype=np.uint8)
    Bv = np.frombuffer(b.encode(), dtype=np.uint8)
    NEG = -(1 << 40)
    o, e = GAP_OPEN, GAP_EXT
    ar = np.arange(n + 1, dtype=np.int64)
    H_prev = -(o + e * ar)
    H_prev[0] = 0
    F_prev = np.full(n + 1, NEG, dtype=np.int64)
    ptrH = np.zeros((m + 1, n + 1), dtype=np.uint8)
    ptrF = np.zeros((m + 1, n + 1), dtype=np.uint8)
    ptrE = np.zeros((m + 1, n + 1), dtype=np.uint8)
    ptrH[0, 1:] = 2
    ptrE[0, 2:] = 1
    D = np.full(n + 1, NEG, dtype=np.int64)
    E = np.full(n + 1, NEG, dtype=np.int64)
    for i in range(1, m + 1):
        f_open, f_ext = H_prev - o - e, F_prev - e
        F = np.maximum(f_open, f_ext)
        ptrF[i] = f_ext > f_open
        D[1:] = H_prev[:-1] + np.where(Bv == A[i - 1], MATCH, MISMATCH)
        Hp = np.maximum(D, F)
        pm = np.maximum.accumulate(Hp + e * ar)
        E[1:] = pm[:-1] - o - e * ar[1:]
        H = np.maximum(Hp, E)
        ptr = np.where(D >= F, 0, 1).astype(np.uint8)
        ptr[E > Hp] = 2
        ptr[0] = 1
        ptrH[i] = ptr
        ptrE[i, 1:] = (E[:-1] - e) > (H[:-1] - o - e)
        H_prev, F_prev = H, F
    amap = np.full(m, -1, dtype=np.int64)
    i, j, state = m, n, 0
    while i > 0 or j > 0:
        if state == 0:
            if i == 0:
                state = 2
                continue
            if j == 0:
                state = 1
                continue
            p = ptrH[i, j]
            if p == 0:
                amap[i - 1] = j - 1
                i, j = i - 1, j - 1
            else:
                state = int(p)
        elif state == 1:
            ext = ptrF[i, j]
            i -= 1
            if not ext:
                state = 0
        else:
            ext = ptrE[i, j]
            j -= 1
            if not ext:
                state = 0
    return amap


def slide_range(seq, start, length, kind, alt):
    """Range [lo, hi) of seq over which an INDEL (DEL of seq[start:start+length] / INS of alt at start) can slide."""
    if kind == "DEL":
        lo, hi = start, start + length
        while lo > 0 and seq[lo - 1] == seq[hi - 1]:
            lo, hi = lo - 1, hi - 1
        lo2, hi2 = start, start + length
        while hi2 < len(seq) and seq[lo2] == seq[hi2]:
            lo2, hi2 = lo2 + 1, hi2 + 1
        return lo, hi2
    # INS of alt before seq[start]
    lo, rot = start, alt
    while lo > 0 and seq[lo - 1] == rot[-1]:
        rot, lo = rot[-1] + rot[:-1], lo - 1
    hi, rot = start, alt
    while hi < len(seq) and seq[hi] == rot[0]:
        rot, hi = rot[1:] + rot[0], hi + 1
    return lo, hi


class Projector:
    def __init__(self, graph, genome, max_bp, max_nodes, max_span, max_cells):
        self.g, self.genome = graph, genome
        self.max_bp, self.max_nodes, self.max_span, self.max_cells = max_bp, max_nodes, max_span, max_cells

    def bubble(self, path, node, last):
        left, lparent = self.g.search(2 * node + 1, path, self.max_bp, self.max_nodes)
        right, rparent = self.g.search(2 * last, path, self.max_bp, self.max_nodes)
        if not left:
            return None, "no_left_anchor"
        if not right:
            return None, "no_right_anchor"
        pair = best_pair(path, left, right, self.max_span)
        if pair is None:
            return None, "inconsistent_anchors"
        _, sl, dl, sr, dr, plus, g0, g1 = pair
        lwalk = [s ^ 1 for s in chain(lparent, sl)[::-1]]  # forward orientation: left anchor ... node
        rwalk = chain(rparent, sr)  # last ... right anchor
        return dict(left=lwalk[0], right=rwalk[-1], before=lwalk[1:-1], after=rwalk[1:-1], dl=dl, dr=dr,
                    plus=plus, g0=g0, g1=g1, contig=path.contigs[int(path.chrom[lwalk[0] >> 1])]), None

    def project(self, rec):
        """dict of the projection fields for one unplaced graph record."""
        info = parse_info(rec[7])
        node, start, kind = int(rec[0]), int(info["START"]), info["KIND"]
        evref, evalt = (info.get("EVREF") or "").upper(), (info.get("EVALT") or "").upper()
        path_steps = [tuple(int(v) for v in s.split(":")) for s in info["PATH"].split(",")] if info.get("PATH") else []
        out = dict(event_kind=kind, evref=evref, evalt=evalt, path=info.get("PATH") or "", node=node, offset=start,
                   proj_detail="", near_gap="", grch38_chrom="", grch38_pos="", grch38_ref="", grch38_alt="",
                   grch38_start0="", grch38_end0="", anchor_left="", anchor_right="", anchor_bp_left="",
                   anchor_bp_right="", branch_len="", grch38_span="", chm13_chrom="", chm13_pos="", chm13_class="none",
                   chm13_end="")
        out.update(self.chm13(node, start, evref, evalt, kind, path_steps))
        last = path_steps[-1][0] if path_steps else node
        b, why = self.bubble(self.g.ref, node, last)
        if b is None:
            out.update(proj_class="iii", proj_detail=why)
            return out
        out.update(anchor_left=f"{b['left'] >> 1}{'+-'[b['left'] & 1]}", anchor_right=f"{b['right'] >> 1}{'+-'[b['right'] & 1]}",
                   anchor_bp_left=b["dl"], anchor_bp_right=b["dr"], grch38_span=b["g1"] - b["g0"])
        interior = b["before"] + [2 * node] + [2 * n for n, _, _ in path_steps] + b["after"]
        seqs = self.g.sequences([s >> 1 for s in interior])
        parts = [side_rc(seqs[s >> 1], s) for s in interior]
        branch = "".join(parts)
        out["branch_len"] = len(branch)
        p_node = sum(len(x) for x in parts[:len(b["before"])])  # branch index of the node's first base
        contig, g0, g1 = b["contig"], b["g0"], b["g1"]
        # event on the branch (forward walk): SNP/DEL [p, p + k), INS boundary p
        if kind == "SNP":
            p, k = p_node + start, 1
        elif kind == "DEL":
            p, k = p_node + start, len(evref)
        else:
            p, k = p_node + start, 0
        ref_b, alt_b = evref, evalt
        if kind != "INS" and branch[p:p + k] != evref:
            raise ValueError(f"{rec[2]}: branch bases {branch[p:p + k]!r} at {p} are not EVREF {evref!r}")
        if not b["plus"]:
            branch = rc(branch)
            p = len(branch) - p - k
            ref_b, alt_b = rc(evref), rc(evalt)
        fl = min(FLANK, g0)
        gseq = self.genome.get(contig, g0 - fl, g1 + FLANK)
        fr = len(gseq) - fl - (g1 - g0)
        if len(gseq) != fl + (g1 - g0) + fr or fr < 0:
            raise ValueError(f"{rec[2]}: FASTA interval {contig}:{g0 - fl}-{g1 + FLANK}")
        bext = gseq[:fl] + branch + gseq[len(gseq) - fr:]
        pe = p + fl
        if (len(bext) + 1) * (len(gseq) + 1) > self.max_cells:
            out.update(proj_class="ii", proj_detail="too_long", grch38_chrom=contig, grch38_start0=g0, grch38_end0=g1,
                       grch38_pos=g0 + 1)
            return out
        amap = global_align(bext, gseq)
        gpos = lambda j: g0 - fl + int(j)  # noqa: E731  (G_ext index -> GRCh38 0-based)
        gap_cols = np.flatnonzero(amap < 0)
        # also gaps in G (deleted GRCh38 bases): consecutive mapped b indices with a jump in G
        mapped = np.flatnonzero(amap >= 0)
        jumps = mapped[1:][np.diff(amap[mapped]) > 1] if len(mapped) > 1 else np.zeros(0, dtype=np.int64)
        gap_idx = np.concatenate([gap_cols, jumps])
        near = bool(np.any(np.abs(gap_idx - pe) <= NEAR_GAP)) if len(gap_idx) else False
        out["near_gap"] = int(near)

        def interval_of(lo, hi):
            """GRCh38 [start0, end0) between the aligned neighbours of bext[lo:hi]."""
            left = [j for j in amap[:lo] if j >= 0]
            right = [j for j in amap[hi:] if j >= 0]
            s0 = gpos(left[-1]) + 1 if left else g0
            e0 = gpos(right[0]) if right else g1
            return min(s0, e0), max(s0, e0)

        def as_ii(detail, lo, hi):
            s0, e0 = interval_of(lo, hi)
            out.update(proj_class="ii", proj_detail=detail, grch38_chrom=contig, grch38_start0=s0, grch38_end0=e0,
                       grch38_pos=s0 + 1)
            return out

        def gap_free(lo, hi):
            js = amap[max(0, lo):min(len(amap), hi)]
            return bool((js >= 0).all() and (len(js) < 2 or (np.diff(js) == 1).all()))

        if kind == "SNP":
            j = amap[pe]
            if j < 0:
                return as_ii("in_insertion", pe, pe + 1)
            pos0 = gpos(j)
            fasta = gseq[j]
            if alt_b == fasta:
                out.update(proj_class="ii", proj_detail="alt_is_grch38_base", grch38_chrom=contig, grch38_start0=pos0,
                           grch38_end0=pos0 + 1, grch38_pos=pos0 + 1)
                return out
            if not gap_free(pe - NEAR_GAP, pe + 1 + NEAR_GAP):
                return as_ii("near_branch_indel", pe - NEAR_GAP, pe + 1 + NEAR_GAP)
            detail = "exact" if bext[pe] == fasta else "branch_base_differs"
            out.update(proj_class="i", proj_detail=detail, grch38_chrom=contig, grch38_pos=pos0 + 1, grch38_ref=fasta,
                       grch38_alt=alt_b, grch38_start0=pos0, grch38_end0=pos0 + 1)
            return out
        # INDEL: the sliding range on the branch, +-NEAR_GAP, must be aligned without gaps
        lo, hi = slide_range(bext, pe, k, kind, alt_b)
        if kind == "DEL":
            in_ins = bool((amap[pe:pe + k] < 0).any())
        else:
            in_ins = bool(amap[pe - 1] < 0 or (pe < len(amap) and amap[pe] < 0))
        if in_ins:
            return as_ii("in_insertion", lo, hi)
        if not gap_free(lo - NEAR_GAP, hi + NEAR_GAP):
            return as_ii("near_branch_indel", lo - NEAR_GAP, hi + NEAR_GAP)
        clo, chi = max(0, lo - CONTEXT), min(len(bext), hi + CONTEXT)
        js = amap[clo:chi]
        mism = sum(x != y for x, y in zip(bext[clo:chi], gseq[js[0]:js[-1] + 1]))
        if kind == "DEL":
            pos0 = gpos(amap[pe])
            gref = gseq[amap[pe]:amap[pe] + k]
            gpos0, gr, ga, moved = self.genome.left_align(contig, pos0, gref, "", "DEL")
        else:
            pos0 = gpos(amap[pe]) if pe < len(amap) else gpos(amap[pe - 1]) + 1
            gpos0, gr, ga, moved = self.genome.left_align(contig, pos0, "", alt_b, "INS")
        vpos, vref, valt = self.genome.vcf_allele(contig, gpos0, gr, ga, kind)
        out.update(proj_class="i", proj_detail="exact" if mism == 0 else "composite", grch38_chrom=contig,
                   grch38_pos=vpos, grch38_ref=vref, grch38_alt=valt, grch38_start0=gpos0,
                   grch38_end0=gpos0 + len(gr), leftshift=moved)
        return out

    def chm13(self, node, start, evref, evalt, kind, path_steps):
        c = self.g.chm13
        if c is None:
            return dict(chm13_class="not_built")
        lin = c.linear(node, start, evref, evalt, kind, path_steps or None)
        if lin is not None:
            return dict(chm13_class="direct", chm13_chrom=lin["chrom"], chm13_pos=lin["pos0"] + 1,
                        chm13_end=lin["pos0"] + max(len(lin["ref"]), 1))
        last = path_steps[-1][0] if path_steps else node
        b, why = self.bubble(c, node, last)
        if b is None:
            return dict(chm13_class="none")
        return dict(chm13_class="anchor_interval", chm13_chrom=b["contig"], chm13_pos=b["g0"] + 1, chm13_end=b["g1"])


def graph_of(header):
    gi = meta(header, "pansoma_graph_index")
    name = GRAPHS.get(Path(gi).name)
    if name is None:
        raise SystemExit(f"unknown graph index {gi}")
    return name, gi


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("--vcf-dir", required=True)
    ap.add_argument("--run", required=True)
    ap.add_argument("--output", required=True)
    ap.add_argument("--max-bp", type=int, default=10000, help="anchor search limit: bases walked on each side")
    ap.add_argument("--max-nodes", type=int, default=5000, help="anchor search limit: nodes expanded on each side")
    ap.add_argument("--max-span", type=int, default=100000, help="largest GRCh38 span between the anchors")
    ap.add_argument("--max-cells", type=int, default=40_000_000, help="largest alignment matrix")
    args = ap.parse_args()
    t0 = time.time()
    vdir, out = Path(args.vcf_dir), Path(args.output)
    out.mkdir(parents=True, exist_ok=True)
    graph_vcf = next(vdir.glob("*.graph.vcf.gz"))
    name = graph_vcf.name[:-len(".graph.vcf.gz")]
    linear_vcf = vdir / f"{name}.linear.vcf.gz"
    unplaced_vcf = vdir / f"{name}.linear.unplaced.vcf.gz"
    gheader, _, _ = read_vcf(graph_vcf)
    gname, graph_index = graph_of(gheader)
    tensors = Path(meta(gheader, "pansoma_tensors"))
    kind_all = meta(gheader, "pansoma_kind")
    labels = json.loads((tensors / kind_all / "labels.manifest.json").read_text())
    set_name = tensors.parent.name
    uheader, ucolumns, urecs = read_vcf(unplaced_vcf)
    urecs = list(urecs)
    stats = dict(run=args.run, set=set_name, kind=kind_all, graph=gname, vcf_dir=str(vdir.resolve()),
                 unplaced_records=len(urecs), unplaced_pass=sum(r[6] == "PASS" for r in urecs))
    genome = Genome(labels["fasta"])
    rows = []
    if gname == "linear" or not urecs:
        graph = None
    else:
        graph = Graph(gname, graph_index, labels["reference_path"])
        proj = Projector(graph, genome, args.max_bp, args.max_nodes, args.max_span, args.max_cells)
    for k, rec in enumerate(urecs):
        info = parse_info(rec[7])
        p = proj.project(rec) if graph is not None else dict(proj_class="iii", proj_detail="graph_not_prepared",
                                                               chm13_class="not_built", node=int(rec[0]),
                                                               offset=int(info["START"]), event_kind=info["KIND"])
        p.update(run=args.run, set=set_name, kind=kind_all, chrom_graph_id=rec[0], ref=rec[3], alt=rec[4],
                 p_somatic=info.get("OFFREF_P_SOMATIC", info.get("P_SOMATIC")), candidate_id=rec[2],
                 block=info.get("BLOCK", ""), filter=rec[6], qual=rec[5], af=rec[9].split(":")[-1],
                 alt_reads=a1_reads(rec), _rec=rec)
        rows.append(p)
        if (k + 1) % 500 == 0:
            print(f"{k + 1:,} / {len(urecs):,} projected, {time.time() - t0:.0f} s", flush=True)
    # placed records and the class-i records -> complete VCF
    lheader, lcolumns, lrecs = read_vcf(linear_vcf)
    placed = {}
    for r in lrecs:
        placed[(r[0], int(r[1]), r[3], r[4])] = r
    n_placed = len(placed)
    dups = 0
    for p in rows:
        p["dup_of_placed"] = ""
        if p["proj_class"] != "i":
            continue
        rec = p["_rec"]
        info = parse_info(rec[7])
        key = (p["grch38_chrom"], int(p["grch38_pos"]), p["grch38_ref"], p["grch38_alt"])
        new_info = dict(NODE=rec[0], **info, LEFTSHIFT=p.get("leftshift") or None,
                        PROJ=p["proj_detail"], PROJ_NEARGAP=bool(p["near_gap"]) or None)
        line = [key[0], str(key[1]), rec[2], key[2], key[3], rec[5], rec[6], format_info(new_info), *rec[8:]]
        p["_line"] = line
        if key in placed:
            dups += 1
            p["dup_of_placed"] = placed[key][2]
            old = placed[key]
            score = lambda r: (r[6] == "PASS", float(parse_info(r[7]).get("P_SOMATIC", 0)))  # noqa: E731
            keep, other = (line, old) if score(line) > score(old) else (old, line)
            kinfo = parse_info(keep[7])
            kinfo["MERGED"] = int(kinfo.get("MERGED", 1)) + 1
            keep = keep[:2] + [keep[2] + ";" + other[2]] + keep[3:7] + [format_info(kinfo)] + keep[8:]
            placed[key] = keep
        else:
            placed[key] = line
    order = {c: i for i, c in enumerate(genome.fasta.references)}
    complete = [placed[k] for k in sorted(placed, key=lambda k: (order[k[0]], k[1], k[2], k[3]))]
    info_end = max(i for i, h in enumerate(lheader) if h.startswith(("##INFO", "##FILTER", "##FORMAT")))
    cheader = lheader[:info_end + 1] + EXTRA_INFO.split("\n") + lheader[info_end + 1:]
    cheader = [h for h in cheader if not h.startswith("##pansoma_graph_vcf=")] + [f"##pansoma_graph_vcf={graph_vcf.resolve()}",
                                                                                   f"##E9a_projection={out.resolve()}"]
    write_vcf(out / "complete.linear.vcf.gz", cheader, lcolumns, complete)
    proj_lines = sorted((p["_line"] for p in rows if "_line" in p), key=lambda r: (order[r[0]], int(r[1]), r[3], r[4]))
    write_vcf(out / "projected.grch38.vcf.gz", cheader, lcolumns, proj_lines)
    remaining = [p["_rec"] for p in rows if p["proj_class"] != "i"]
    write_vcf(out / "complete.linear.unplaced.vcf.gz", uheader, ucolumns, remaining)
    with open(out / "projected.tsv", "w") as f:
        f.write("\t".join(TSV_COLUMNS) + "\n")
        for p in rows:
            f.write("\t".join(str(p.get(c, "")) for c in TSV_COLUMNS) + "\n")
    count = lambda rs, key: dict(Counter(key(p) for p in rs))  # noqa: E731
    passing = [p for p in rows if p["filter"] == "PASS"]
    stats.update(
        placed_records=n_placed, projected_class_i=sum(p["proj_class"] == "i" for p in rows), merged_with_placed=dups,
        complete_records=len(complete),
        classes_all=count(rows, lambda p: p["proj_class"]), classes_pass=count(passing, lambda p: p["proj_class"]),
        details_pass=count(passing, lambda p: f"{p['proj_class']}:{p['proj_detail']}"),
        chm13_pass=count(passing, lambda p: p["chm13_class"]),
        chm13_by_class_pass=count(passing, lambda p: f"{p['proj_class']}:{p['chm13_class']}"),
        near_gap_class_i_pass=sum(p["proj_class"] == "i" and p["near_gap"] == 1 for p in passing),
        seconds=time.time() - t0, limits=dict(max_bp=args.max_bp, max_nodes=args.max_nodes, max_span=args.max_span,
                                              max_cells=args.max_cells))
    # vcfeval reads unplaced_pass from <calls>.stats.json: the PASS records still without a GRCh38 variant
    (out / "complete.linear.stats.json").write_text(json.dumps(dict(
        unplaced=len(remaining), unplaced_pass=sum(r[6] == "PASS" for r in remaining),
        note="E9a complete call set: linear.vcf.gz + class-i projections"), indent=2) + "\n")
    (out / "projection.json").write_text(json.dumps(stats, indent=2) + "\n")
    print(json.dumps(stats, indent=2), flush=True)


if __name__ == "__main__":
    main()
