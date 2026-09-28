"""Germline / somatic truth VCFs -> graph-node candidate keys -> one label per merged tensor.

1. Truth alleles. Every ALT of every record on chr1-22 is split off, trimmed of shared
   prefix/suffix bases, MNVs become one SNV per differing base; complex replacements are
   kept only as "nearby truth" (they cannot equal a candidate).
2. Equivalent positions. An indel in a repeat has several equivalent placements (vg puts
   it at one end of the repeat in read orientation, VCFs are left-aligned); all of them,
   shifted left and right on GRCh38, are enumerated.
3. Node keys. Each placement is converted, through the GRCh38 path (reference_path), into
   the builder's candidate identity `node:start:KIND:REF>ALT` in node-forward coordinates:
   reverse-oriented nodes flip the position and reverse-complement the bases; an insertion
   at a node junction gets a key on both nodes; a deletion over consecutive reference nodes of
   one orientation gets the builder's multi-node key (`node:start:DEL:REF>@n2+n3`); nodes the
   reference visits more than once give no key.
4. Labels. A tensor's representative allele (its candidate_id, A1) is looked up:

       1 somatic    A1 is a somatic truth allele with FILTER PASS/. (inside or outside the BED), or overlaps one by
                    more than MIN_OVERLAP (partial, below)
       2 germline   A1 is a germline truth allele with FILTER PASS/. or only GAP1/GAP2 (dipcall: on one assembled
                    haplotype, the other uncalled), inside or outside the BED, or overlaps one (partial)
       0 non        a tensor whose site has another allele that is a truth allele as in 1/2 which A1 does not
                    overlap enough ("allele", below; inside or outside the BED), and every other tensor on a GRCh38
                    node inside somatic BED ∩ germline BED: no truth allele, a truth allele nearby that A1 does not
                    overlap enough (errors and artifacts next to real variants), a germline allele with another
                    FILTER (dipcall HET1/HET2)
      -1 ignore     a filtered somatic truth allele (A1's or another allele's), outside the BEDs, no position (no
                    unique GRCh38 visit and no reference node within ANCHOR_REACH IDs, the two reference neighbours
                    on different contigs, or more than ANCHOR_GAP GRCh38 bases between them), an off-reference node
                    without a partial somatic match (not a training negative: most branch nodes carry no truth at
                    all), or, whatever the truth, a tensor with AF below the label run's floor of its kind: an SNV
                    below snv_min_af (short-read sets: 0.07), an INDEL (A1 an INS or DEL) below indel_min_af

   The first rule that applies decides (classify, in this order): the AF floors (snv_min_af for SNVs, indel_min_af
   for INDELs); a truth allele of A1 (somatic, then germline PASS/GAP); a truth allele as in 1/2 of another allele
   ("allele" partial, else 0), else a filtered somatic one (-1); no position; outside the BEDs; A1 a germline allele
   with another FILTER (0); a "residual" partial match (somatic, then germline; an off-reference node: somatic only,
   else -1); otherwise 0.

   Partial (labels.ndjson `partial`, `overlap`, `partial_truth`): the same event written differently by the
   graph alignment. "allele": another allele of the site is the truth allele and A1 overlaps it by more than
   MIN_OVERLAP (else 0). "residual": a truth allele at the same place (spans of equivalent placements
   intersecting or adjacent; for an off-reference node, within NEAR_BP of the interval between its reference
   neighbours by node ID) that A1 overlaps by more than MIN_OVERLAP. Overlap: DEL/DEL the deleted bases in
   common, INS/INS the inserted bases in common at the same insertion point, over the longer allele; for
   somatic truths also the haplotype overlap of the A1 reads (haplotype_overlap: a branch plus a residual
   edit, a skipped node plus a mismatch), which is how off-reference nodes are matched.

Outputs next to the merged shards (per chromosome, same order as <chrom>_variant_summary.ndjson):
<chrom>_shard_NNNNN_labels.npy (int8) and <chrom>_labels.ndjson; labels.manifest.json (counts, reasons, the rule
constants, the AF floors snv_min_af / indel_min_af (null when unset), rules_sha256 = SHA-256 of this file, provenance);
<truth dir>/<set>.graph.tsv (every truth allele with its keys) and <set>.recall.tsv / recall.json.
"""
from collections import Counter, defaultdict
from datetime import datetime, timezone
import json
from pathlib import Path
import re

import numpy as np

from ..common import format_name, read_json, sha256_file, write_json
from .chr_index import AUTOSOMES
from .merge_shards import LAYOUT
from .reference_path import ReferencePath, rc

FORMAT = "truth-labels"  # the rules are the module docstring and the constants below; labels.manifest.json
# identifies them by this file's SHA-256 (rules_sha256).
# dipcall FILTER values of a germline allele present on one assembled haplotype while the other is uncalled: the
# allele is in the normal genome (zygosity unknown), so it counts as germline truth.
GAP_FILTERS = frozenset({"GAP1", "GAP2"})
LABELS = {"ignore": -1, "non": 0, "somatic": 1, "germline": 2}
NEAR_BP = 10
MIN_OVERLAP = 0.45  # a tensor overlapping a truth allele by more than this counts as that truth (partial)
ANCHOR_REACH = 200  # node IDs searched on each side for the reference nodes around an off-reference node
ANCHOR_GAP = 1024  # most GRCh38 bases between those reference nodes (the reference a branch replaces; the nodes'
# own lengths, up to 1024 bp each, not counted) for the interval to be taken as the node's position
EVIDENCE_ROWS = 20  # A1 reads (and half as many REF reads) compared with the haplotypes
MIN_EVIDENCE_BASES = 30
HAPLOTYPE_WINDOW = 90  # GRCh38 bases on each side of a truth allele for the haplotype comparison
MAX_SHIFTS = 5000
CANDIDATE = re.compile(rb'"candidate_id": "([^"]+)"')
REASONS = re.compile(rb'"reasons": \[([^\]]*)\]')


# --- truth alleles ------------------------------------------------------------------

def split_alleles(pos, ref, alts):
    """(pos0, ref, alt, kind, allele_index) per ALT after trimming; kind SNP | INS | DEL | COMPLEX."""
    for index, alt in enumerate(alts or (), 1):
        if not alt or alt.startswith("<") or "*" in alt or "." == alt:
            continue
        r, a, p = ref.upper(), alt.upper(), pos - 1
        if "N" in r or "N" in a:
            continue
        while len(r) > 1 and len(a) > 1 and r[-1] == a[-1]:
            r, a = r[:-1], a[:-1]
        while r and a and r[0] == a[0]:
            r, a, p = r[1:], a[1:], p + 1
        if not r and not a:
            continue
        if len(r) == len(a):
            for k, (x, y) in enumerate(zip(r, a)):
                if x != y:
                    yield p + k, x, y, "SNP", index
        elif not r:
            yield p, "", a, "INS", index
        elif not a:
            yield p, r, "", "DEL", index
        else:
            yield p, r, a, "COMPLEX", index


def placements(genome, pos, ref, alt, kind):
    """Every equivalent (pos0, ref, alt) of an allele on `genome` (the contig sequence, upper case)."""
    if kind == "SNP":
        return [(pos, ref, alt)]
    if kind == "DEL":
        m = len(ref)
        while pos > 0 and genome[pos - 1] == ref[-1]:
            pos, ref = pos - 1, genome[pos - 1] + ref[:-1]
        result = [(pos, ref, "")]
        while len(result) < MAX_SHIFTS and pos + m < len(genome) and genome[pos + m] == ref[0]:
            pos, ref = pos + 1, ref[1:] + genome[pos + m]
            result.append((pos, ref, ""))
        return result
    if kind == "INS":
        while pos > 0 and genome[pos - 1] == alt[-1]:
            pos, alt = pos - 1, alt[-1] + alt[:-1]
        result = [(pos, "", alt)]
        while len(result) < MAX_SHIFTS and pos < len(genome) and genome[pos] == alt[0]:
            pos, alt = pos + 1, alt[1:] + alt[0]
            result.append((pos, "", alt))
        return result
    return []


class Locator:
    """GRCh38 contig position -> the unique-visit node covering it (via the reference walk)."""

    def __init__(self, path):
        self.path = path
        self.cache = {}

    def contig(self, chrom):
        """Walk arrays of one contig; only the current contig is kept (inputs are chromosome-sorted)."""
        if chrom not in self.cache:
            if chrom not in self.path.contig_index:
                return None
            else:
                part = self.path.contig_slice(chrom)
                nodes = np.asarray(self.path.path_nodes[part])
                self.cache = {chrom: (nodes, np.asarray(self.path.path_starts[part]),
                                      np.asarray(self.path.path_reverse[part]),
                                      np.asarray(self.path.lengths)[nodes].astype(np.int64),
                                      np.asarray(self.path.visits)[nodes] == 1)}
        return self.cache[chrom]

    def at(self, chrom, x):
        """(node, start, length, reverse) of the node holding base x, or None (off contig or ambiguous node)."""
        data = self.contig(chrom)
        if data is None:
            return None
        nodes, starts, reverse, lengths, unique = data
        k = int(np.searchsorted(starts, x, side="right")) - 1
        if k < 0 or x >= starts[k] + lengths[k] or not unique[k]:
            return None
        return int(nodes[k]), int(starts[k]), int(lengths[k]), bool(reverse[k])

    def keys(self, chrom, pos, ref, alt, kind):
        """Candidate IDs of one linear placement."""
        if kind == "INS":
            hits = [h for h in (self.at(chrom, pos), self.at(chrom, pos - 1) if pos > 0 else None) if h]
            result = set()
            for node, start, length, reverse in hits:
                t = pos - start
                if 0 <= t <= length:
                    fwd = length - t if reverse else t
                    result.add(f"{node}:{fwd}:INS:>{rc(alt) if reverse else alt}")
            return sorted(result)
        if kind == "DEL":
            return self.deletion_keys(chrom, pos, ref)
        hit = self.at(chrom, pos)
        if hit is None:
            return []
        node, start, length, reverse = hit
        fwd = start + length - pos - 1 if reverse else pos - start
        if reverse:
            ref, alt = rc(ref), rc(alt)
        return [f"{node}:{fwd}:{kind}:{ref}>{alt}"]

    def deletion_keys(self, chrom, pos, ref):
        """Key of a deletion of linear [pos, pos + len(ref)), also across consecutive reference nodes.

        Same identity as the builder's Candidate.metadata(): the nodes in forward order (their own
        forward strand), `start` on the forward-first node, REF over all of them, and further
        nodes as an "@n2+n3" suffix. The builder only joins deletions over mappings of one
        orientation, so a span over nodes of mixed orientation on the reference has no key.
        """
        segments, x, end = [], pos, pos + len(ref)
        while x < end:
            hit = self.at(chrom, x)
            if hit is None:
                return []
            node, start, length, reverse = hit
            b = min(end, start + length)
            segments.append((node, start, length, reverse, x, b))
            x = b
        if len({s[3] for s in segments}) != 1:
            return []
        reverse = segments[0][3]
        forward = [(node, length - (b - start) if reverse else a - start, length - (a - start) if reverse else b - start)
                   for node, start, length, _, a, b in segments]
        if reverse:
            forward.reverse()
            ref = rc(ref)
        via = "@" + "+".join(str(node) for node, _, _ in forward[1:]) if len(forward) > 1 else ""
        return [f"{forward[0][0]}:{forward[0][1]}:DEL:{ref}>{via}"]


class Bed:
    """Merged half-open intervals per chromosome."""

    def __init__(self, intervals):
        self.data = {}
        for chrom, items in intervals.items():
            items = sorted(items)
            merged = []
            for s, e in items:
                if merged and s <= merged[-1][1]:
                    merged[-1][1] = max(merged[-1][1], e)
                else:
                    merged.append([s, e])
            self.data[chrom] = (np.array([m[0] for m in merged], dtype=np.int64),
                                np.array([m[1] for m in merged], dtype=np.int64))

    @classmethod
    def read(cls, path):
        intervals = defaultdict(list)
        with open(path) as stream:
            for line in stream:
                if line.strip() and not line.startswith(("#", "track", "browser")):
                    f = line.split("\t")
                    intervals[f[0]].append((int(f[1]), int(f[2])))
        return cls(intervals)

    def intersect(self, other):
        result = defaultdict(list)
        for chrom in self.data.keys() & other.data.keys():
            (a0, a1), (b0, b1) = self.data[chrom], other.data[chrom]
            i = j = 0
            while i < len(a0) and j < len(b0):
                s, e = max(a0[i], b0[j]), min(a1[i], b1[j])
                if s < e:
                    result[chrom].append((int(s), int(e)))
                if a1[i] < b1[j]:
                    i += 1
                else:
                    j += 1
        return Bed(result)

    def contains(self, chrom, s, e):
        if chrom not in self.data:
            return False
        starts, ends = self.data[chrom]
        k = int(np.searchsorted(starts, s, side="right")) - 1
        return bool(k >= 0 and ends[k] >= e)

    def size(self):
        return int(sum((e - s).sum() for s, e in self.data.values()))


class TruthSet:
    """Truth alleles of one VCF with their node keys and linear spans."""

    def __init__(self, name, vcf, bed, fasta, locator, chromosomes=AUTOSOMES):
        import pysam
        self.name, self.vcf, self.bed_path = name, Path(vcf), Path(bed)
        self.bed = Bed.read(bed)
        self.alleles, self.by_key = [], defaultdict(list)
        spans = defaultdict(list)
        self.stats = Counter()
        genome_file = pysam.FastaFile(str(fasta))
        genome, genome_chrom = None, None
        with pysam.VariantFile(str(vcf)) as source:
            sample = list(source.header.samples)[0] if source.header.samples else None
            for record in source:
                if record.chrom not in chromosomes:
                    self.stats["records_other_contigs"] += 1
                    continue
                if record.chrom != genome_chrom:
                    genome, genome_chrom = genome_file.fetch(record.chrom).upper(), record.chrom
                self.stats["records"] += 1
                filters = list(record.filter.keys())
                passed = not filters or filters == ["PASS"]
                gt = record.samples[sample].get("GT") if sample else None
                for pos, ref, alt, kind, index in split_alleles(record.pos, record.ref, record.alts):
                    self.stats[f"alleles_{kind}"] += 1
                    if kind in ("SNP", "DEL") and genome[pos:pos + len(ref)] != ref:
                        self.stats["reference_mismatch"] += 1
                        continue
                    places = placements(genome, pos, ref, alt, kind)
                    keys = sorted({k for p in places for k in locator.keys(record.chrom, *p, kind)})
                    lo = min((p[0] for p in places), default=pos)
                    hi = max((p[0] + len(p[1]) for p in places), default=pos + len(ref))
                    tid = len(self.alleles)
                    self.alleles.append(dict(truth_id=tid, chrom=record.chrom, vcf_pos=record.pos, vcf_ref=record.ref,
                                             vcf_alt=record.alts[index - 1], pos0=pos, ref=ref, alt=alt, kind=kind,
                                             gt="/".join("." if g is None else str(g) for g in gt) if gt else ".",
                                             phased=bool(sample and record.samples[sample].phased),
                                             allele_index=index, filters=filters or ["PASS"], passed=passed,
                                             in_bed=self.bed.contains(record.chrom, lo if kind != "INS" else lo - 1,
                                                                      max(hi, lo + 1)),
                                             placements=len(places), keys=keys, lo=lo, hi=hi,
                                             left_alt=places[0][2] if kind == "INS" and places else alt))
                    # Every equivalent placement counts as "the allele is here" (an INS span runs from its
                    # leftmost to its rightmost boundary: reads may write it anywhere in the repeat).
                    spans[record.chrom].append((lo, max(hi, lo + 1), tid))
                    for key in keys:
                        self.by_key[key].append(tid)
                    self.stats["alleles_with_keys" if keys else "alleles_without_keys"] += 1
        genome_file.close()
        self.spans, self.index = {}, {}
        for chrom, items in spans.items():
            items.sort()
            lo = np.array([i[0] for i in items], dtype=np.int64)
            hi = np.array([i[1] for i in items], dtype=np.int64)
            self.spans[chrom] = (lo, np.maximum.accumulate(hi))
            self.index[chrom] = (lo, hi, np.array([i[2] for i in items], dtype=np.int64), np.maximum.accumulate(hi))

    def overlapping(self, chrom, s, e, pad=0):
        """Truth ids whose span of equivalent placements intersects [s - pad, e + pad)."""
        if chrom not in self.index:
            return []
        lo, hi, tids, maxhi = self.index[chrom]
        found, j = [], int(np.searchsorted(lo, e + pad, side="left")) - 1
        while j >= 0 and maxhi[j] > s - pad:
            if hi[j] > s - pad:
                found.append(int(tids[j]))
            j -= 1
        return found

    def near(self, chrom, s, e, distance=NEAR_BP):
        """Any truth allele span within `distance` bp of [s, e)."""
        if chrom not in self.spans:
            return False
        lo, maxhi = self.spans[chrom]
        k = int(np.searchsorted(lo, e + distance, side="left"))
        return bool(k > 0 and maxhi[k - 1] > s - distance)

    def write_table(self, path):
        columns = ("truth_id", "chrom", "vcf_pos", "vcf_ref", "vcf_alt", "pos0", "ref", "alt", "kind", "gt", "phased",
                   "allele_index", "filters", "passed", "in_bed", "placements", "keys")
        with open(path, "w") as out:
            out.write("\t".join(columns) + "\n")
            for a in self.alleles:
                out.write("\t".join(",".join(a[c]) if isinstance(a[c], list) else str(a[c]) for c in columns) + "\n")

    def summary(self, tid, candidate_id, representative):
        a = self.alleles[tid]
        return dict(truth_id=tid, chrom=a["chrom"], vcf_pos=a["vcf_pos"], vcf_ref=a["vcf_ref"], vcf_alt=a["vcf_alt"],
                    gt=a["gt"], filters=a["filters"], in_bed=a["in_bed"], matched_candidate_id=candidate_id,
                    representative=representative)


# --- partial matches: overlap with a truth allele ------------------------------------

def lcs(a, b):
    """Length of the longest common subsequence of two short strings."""
    row = [0] * (len(b) + 1)
    for x in a:
        diagonal = 0
        for j, y in enumerate(b, 1):
            diagonal, row[j] = row[j], diagonal + 1 if x == y else max(row[j], row[j - 1])
    return row[-1]


def allele_overlap(kind, pos0, ref, alt, truth):
    """Fraction of two alleles in common, positions included (the truth allele over all its equivalent
    placements, lo..hi): DEL/DEL the deleted bases both remove, INS/INS the inserted bases both carry at the
    same insertion point, over the longer allele. Other pairs 0 (a different SNV base; different kinds are
    compared by haplotype_overlap)."""
    if kind != truth["kind"] or kind not in ("DEL", "INS"):
        return 0.0
    if kind == "DEL":
        shared = min(len(ref), len(truth["ref"]), max(0, min(pos0 + len(ref), truth["hi"]) - max(pos0, truth["lo"])))
        return shared / max(len(ref), len(truth["ref"]))
    if not truth["lo"] <= pos0 <= truth["hi"]:
        return 0.0
    left = truth["left_alt"]
    k = (pos0 - truth["lo"]) % len(left)  # the truth insertion rotates as it moves through a repeat
    return lcs(alt, left[k:] + left[:k]) / max(len(alt), len(left))


def fit_distance(read, text):
    """Edits aligning `read` entirely inside `text` (free ends in `text`)."""
    t = np.frombuffer(text.encode(), dtype=np.uint8)
    positions = np.arange(len(t) + 1, dtype=np.int32)
    previous = np.zeros(len(t) + 1, dtype=np.int32)
    for i, base in enumerate(read.encode(), 1):
        current = np.empty_like(previous)
        current[0] = i
        current[1:] = np.minimum(previous[:-1] + (t != base), previous[1:] + 1)
        previous = np.minimum.accumulate(current - positions) + positions
    return int(previous.min())


def haplotype_overlap(alt_rows, ref_rows, reference, haplotype, event, oriented=True):
    """Median over the A1 reads of the truth event they carry, from their local sequences (the tensor window):
    with d_ref, d_truth = edits of a read against the GRCh38 window and against it with the truth allele, both
    minus the REF reads' median error b, shared = (d_ref + event - d_truth) / 2 and overlap = shared /
    max(d_ref, event). A read not closer to the truth haplotype than to GRCh38 before that subtraction
    (d_truth >= d_ref) carries none of the event: overlap 0 (the subtraction alone would give it 0.5 whenever
    b >= both distances). Catches a truth written differently by the graph alignment (a branch plus a residual
    edit, a skipped node plus a mismatch). `oriented=False` (off-reference node): each read is taken in the
    orientation that fits better. None with fewer than 3 A1 reads."""
    def fits(read):
        pairs = [(fit_distance(q, reference), fit_distance(q, haplotype)) for q in ([read] if oriented else [read, rc(read)])]
        return min(pairs, key=min)
    if len(alt_rows) < 3:
        return None
    base = sorted(min(fits(q)) for q in ref_rows)
    b = base[len(base) // 2] if base else 0
    overlaps = []
    for q in alt_rows:
        raw_ref, raw_truth = fits(q)
        if raw_truth >= raw_ref:  # no closer to the truth than to GRCh38: no evidence for it
            overlaps.append(0.0)
            continue
        d_ref, d_truth = max(0, raw_ref - b), max(0, raw_truth - b)
        shared = max(0.0, (d_ref + event - d_truth) / 2)
        overlaps.append(shared / max(d_ref, event))
    return sorted(overlaps)[len(overlaps) // 2]


def anchor(path, node, reach=ANCHOR_REACH, gap=ANCHOR_GAP):
    """(chrom, start, end) from the nearest unique reference nodes below and above an off-reference node by node ID
    (Minigraph-Cactus numbers nodes in topological order: a branch node's ID lies between its flanks'), the nodes
    included; None without one within `reach` IDs, when the two sides are on different contigs, or when more than
    `gap` GRCh38 bases lie between the two nodes (around centromeres the neighbours by ID can be megabases apart, an
    interval that would match any truth)."""
    sides = []
    below = np.flatnonzero(np.asarray(path.visits[max(1, node - reach):node]) == 1)
    if below.size:
        sides.append(max(1, node - reach) + int(below[-1]))
    above = np.flatnonzero(np.asarray(path.visits[node + 1:node + 1 + reach]) == 1)
    if above.size:
        sides.append(node + 1 + int(above[0]))
    contigs = {int(path.chrom[k]) for k in sides}
    if len(contigs) != 1:
        return None
    nodes = sorted((int(path.start0[k]), int(path.start0[k]) + int(path.lengths[k])) for k in sides)
    if nodes[-1][0] - nodes[0][1] > gap:
        return None
    return path.contigs[contigs.pop()], nodes[0][0], max(end for _, end in nodes)


class Evidence:
    """Read sequences (channel 0) of a merged tensor's A1 and REF rows, and GRCh38 windows."""
    CODES = {1: "A", 2: "C", 3: "G", 4: "T", 5: "N"}

    def __init__(self, directory, fasta):
        import pysam
        self.directory, self.genome, self.loaded = Path(directory), pysam.FastaFile(str(fasta)), (None, None)

    def rows(self, record):
        name = record["shard_file"]
        if self.loaded[0] != name:
            path = self.directory / name
            self.loaded = (name, np.load(path, mmap_mode="r") if path.exists() else None)
        if self.loaded[1] is None:
            return [], []
        bases = self.loaded[1][record["index_within_shard"], 0]
        groups = {g["allele"]: (g["start_row"], g["end_row"]) for g in record.get("row_groups", ())}
        return self.take(bases, groups.get("A1"), EVIDENCE_ROWS), self.take(bases, groups.get("REF"), EVIDENCE_ROWS // 2)

    def take(self, bases, span, limit):
        if not span or span[1] <= span[0]:
            return []
        rows = sorted(set(np.linspace(span[0], span[1] - 1, min(limit, span[1] - span[0])).astype(int).tolist()))
        sequences = ["".join(self.CODES[c] for c in np.asarray(bases[i]).tolist() if 1 <= c <= 5) for i in rows]
        return [q for q in sequences if len(q) >= MIN_EVIDENCE_BASES]

    def window(self, chrom, start, end):
        return self.genome.fetch(chrom, max(0, start), end).upper()


def truth_overlap(truth, tid, record, lin, place, evidence, cache):
    """Overlap (0-1) of a tensor's representative allele with truth allele `tid`: the allele overlap on GRCh38,
    and for somatic truths the haplotype overlap of its A1 reads when the alleles alone do not overlap enough."""
    a = truth.alleles[tid]
    value = allele_overlap(record["event_type"], lin["pos0"], lin["ref"], lin["alt"], a) if lin else 0.0
    if value > MIN_OVERLAP or truth.name != "somatic" or evidence is None or place is None:
        return value
    if "rows" not in cache:
        alt_rows, ref_rows = evidence.rows(record)
        if lin and lin["node_reverse"]:
            alt_rows, ref_rows = [rc(q) for q in alt_rows], [rc(q) for q in ref_rows]
        cache["rows"] = (alt_rows, ref_rows)
    alt_rows, ref_rows = cache["rows"]
    start = max(0, a["pos0"] - HAPLOTYPE_WINDOW)
    reference = evidence.window(place[0], start, a["pos0"] + len(a["ref"]) + HAPLOTYPE_WINDOW)
    o = a["pos0"] - start
    haplotype = reference[:o] + a["alt"] + reference[o + len(a["ref"]):]
    event = 1 if a["kind"] == "SNP" else max(len(a["ref"]), len(a["alt"]))
    h = haplotype_overlap(alt_rows, ref_rows, reference, haplotype, event, oriented=lin is not None)
    return max(value, h or 0.0)


# --- labels -------------------------------------------------------------------------

def linear_interval(lin, kind):
    if kind == "INS":
        return lin["pos0"] - 1, lin["pos0"] + 1
    return lin["pos0"], lin["pos0"] + max(1, len(lin["ref"]))


def classify(record, somatic, germline, path, confident, evidence=None, snv_min_af=None, indel_min_af=None):
    """(label value, label name, reason, details) of one merged summary record (rules: module docstring)."""
    representative = record["candidate_id"]
    alleles = record["alleles"]
    hits = {}
    for truth in (somatic, germline):
        hits[truth.name] = [(a["candidate_id"], tid) for a in alleles for tid in truth.by_key.get(a["candidate_id"], ())]
    details = {truth.name: [truth.summary(tid, cid, cid == representative) for cid, tid in hits[truth.name]]
               for truth in (somatic, germline)}
    lin = path.linear(record["node_id"], record["start"], record["ref"], record["alt"], record["event_type"],
                      record.get("path"))
    details.update(grch38=lin, partial=None)
    # The label run's AF floor of the tensor's kind (A1's event type): as if the build had used it, whatever the truth.
    if snv_min_af is not None and record["event_type"] == "SNP" and record.get("af", 1.0) < snv_min_af:
        return LABELS["ignore"], "ignore", "below_snv_min_af", details
    if indel_min_af is not None and record["event_type"] in ("INS", "DEL") and record.get("af", 1.0) < indel_min_af:
        return LABELS["ignore"], "ignore", "below_indel_min_af", details
    rep_somatic = [tid for cid, tid in hits[somatic.name] if cid == representative]
    rep_germline = [tid for cid, tid in hits[germline.name] if cid == representative]
    counts = {somatic.name: lambda t: somatic.alleles[t]["passed"],
              germline.name: lambda t: germline.alleles[t]["passed"] or set(germline.alleles[t]["filters"]) <= GAP_FILTERS}
    # 1 and 2 hold inside or outside the BEDs: a PASS truth allele (germline: or GAP-filtered, see GAP_FILTERS).
    if rep_somatic:
        if any(counts[somatic.name](t) for t in rep_somatic):
            return LABELS["somatic"], "somatic", "representative_allele_is_somatic_truth", details
        return LABELS["ignore"], "ignore", "somatic_truth_filtered", details
    filtered_germline = False
    if rep_germline:
        if any(germline.alleles[t]["passed"] for t in rep_germline):
            return LABELS["germline"], "germline", "representative_allele_is_germline_truth", details
        if any(counts[germline.name](t) for t in rep_germline):
            return LABELS["germline"], "germline", "representative_allele_is_germline_truth_gap_filtered", details
        filtered_germline = True  # e.g. HET1/HET2: the assembled haplotype itself is ambiguous; not a variant call
    cache = {}

    def partial(kind, name, truth, tids, place):
        scored = max((truth_overlap(truth, t, record, lin, place, evidence, cache), t) for t in tids)
        if scored[0] > MIN_OVERLAP:
            details.update(partial=kind, overlap=round(scored[0], 3), partial_truth=truth.summary(scored[1], None, False))
            return LABELS[name], name, f"{kind}_partial_{name}_truth", details
        return None

    # Another allele of the site is a truth allele: the site is labelled by it when A1 overlaps it > MIN_OVERLAP.
    others = {truth.name: [t for cid, t in hits[truth.name] if cid != representative and counts[truth.name](t)]
              for truth in (somatic, germline)}
    if others[somatic.name] or others[germline.name]:
        place = (lin["chrom"],) if lin else None
        for truth, name in ((somatic, "somatic"), (germline, "germline")):
            if others[truth.name]:
                result = partial("allele", name, truth, others[truth.name], place)
                if result:
                    return result
        return LABELS["non"], "non", "truth_matches_non_representative_allele", details
    if hits[somatic.name]:  # another allele is a filtered somatic truth allele
        return LABELS["ignore"], "ignore", "somatic_truth_filtered", details
    # Position: GRCh38, or for an off-reference node the interval between the reference nodes around it.
    if lin is not None:
        chrom, (s, e), pad = lin["chrom"], linear_interval(lin, record["event_type"]), 1
        inside = confident.contains(chrom, s, e)
    else:
        place = anchor(path, int(record["node_id"]))
        if place is None:
            return LABELS["ignore"], "ignore", "not_on_unique_grch38_node", details
        chrom, s, e = place
        pad, details["anchor"] = NEAR_BP, dict(chrom=chrom, start=s, end=e)
        inside = confident.contains(chrom, (s + e) // 2, (s + e) // 2 + 1)
    # -1 only where test-time calling can drop the same tensors without truth: outside the BED, no position.
    if not inside:
        return LABELS["ignore"], "ignore", "outside_confident_region", details
    if filtered_germline:  # only GRCh38 nodes have truth keys, so never an off-reference node
        return LABELS["non"], "non", "germline_truth_filtered", details
    # A truth allele at the same place that the tensor overlaps by more than MIN_OVERLAP: the same event written
    # differently (somatic: allele or read-haplotype overlap; germline: allele overlap on GRCh38).
    tids = [t for t in somatic.overlapping(chrom, s, e, pad) if counts[somatic.name](t)]
    if tids:
        result = partial("residual", "somatic", somatic, tids, (chrom,))
        if result:
            return result
    # An off-reference node without a partial somatic match is not a training negative: most branch
    # nodes carry no truth at all, so its tensors are -1 rather than 0.
    if lin is None:
        return LABELS["ignore"], "ignore", "off_reference_no_truth_match", details
    tids = [t for t in germline.overlapping(chrom, s, e, pad) if counts[germline.name](t)]
    if tids:
        result = partial("residual", "germline", germline, tids, (chrom,))
        if result:
            return result
    # 0 is every other tensor: artifacts and errors next to real variants included, as test-time calling sees them.
    if somatic.near(chrom, s, e) or germline.near(chrom, s, e):
        return LABELS["non"], "non", "near_truth_allele_mismatch", details
    return LABELS["non"], "non", "confident_no_truth_allele", details


def label_directory(directory, somatic, germline, path, confident, provenance, fasta=None, snv_min_af=None,
                    indel_min_af=None):
    """Write <chrom>_shard_*_labels.npy, <chrom>_labels.ndjson and labels.manifest.json for one merged directory."""
    directory = Path(directory)
    manifest = read_json(directory / "manifest.json")
    if format_name(manifest.get("layout")) != LAYOUT:
        raise ValueError(f"{directory} is not a merged per-chromosome directory")
    counts, reasons, partials, written = {}, Counter(), Counter(), []
    matched = {somatic.name: defaultdict(list), germline.name: defaultdict(list)}
    evidence = Evidence(directory, fasta) if fasta else None
    for chrom, info in manifest["chromosomes"].items():
        counts[chrom] = Counter()
        sizes = {s["file"]: s["tensors"] for s in info["shards"]}
        arrays = {f: np.full(n, -128, dtype=np.int8) for f, n in sizes.items()}
        target = directory / f"{chrom}_labels.ndjson"
        with (directory / info["summary"]).open() as source, target.with_name(target.name + ".tmp").open("w") as out:
            for line in source:
                record = json.loads(line)
                value, name, reason, details = classify(record, somatic, germline, path, confident, evidence, snv_min_af,
                                                        indel_min_af)
                if details["partial"]:
                    partials[f"{details['partial']}_{name}"] += 1
                arrays[record["shard_file"]][record["index_within_shard"]] = value
                counts[chrom][name] += 1
                reasons[reason] += 1
                for truth in (somatic, germline):
                    for d in details[truth.name]:
                        matched[truth.name][d["truth_id"]].append((d["representative"], chrom))
                out.write(json.dumps(dict(candidate_id=record["candidate_id"], site_id=record.get("site_id"),
                                          chrom=chrom, shard_file=record["shard_file"],
                                          index_within_shard=record["index_within_shard"], label=value,
                                          label_name=name, reason=reason, **details)) + "\n")
        for file, values in arrays.items():
            if (values == -128).any():
                raise ValueError(f"{directory / file}: summary does not cover every tensor")
            label_file = directory / file.replace("_data.npy", "_labels.npy")
            np.save(label_file.with_name(label_file.name + ".tmp.npy"), values)
            written.append((label_file.with_name(label_file.name + ".tmp.npy"), label_file))
        written.append((target.with_name(target.name + ".tmp"), target))
    for temporary, final in written:
        temporary.replace(final)
    report = dict(format=FORMAT, rules_sha256=sha256_file(__file__), created=datetime.now(timezone.utc).isoformat(),
                  labels=LABELS, near_bp=NEAR_BP, tensors=sum(sum(c.values()) for c in counts.values()),
                  counts={c: dict(v) for c, v in counts.items()},
                  totals=dict(sum(counts.values(), Counter())), reasons=dict(reasons), partial=dict(partials),
                  min_overlap=MIN_OVERLAP, anchor_reach=ANCHOR_REACH, anchor_gap=ANCHOR_GAP, snv_min_af=snv_min_af,
                  indel_min_af=indel_min_af, **provenance)
    write_json(directory / "labels.manifest.json", report)
    return report, matched


def scan_filtered(directory, keys):
    """{candidate_id: reasons} for filtered candidates whose id is a truth key (one streaming pass)."""
    found = {}
    path = Path(directory) / "filtered_candidates.ndjson"
    if not path.exists():
        return found
    with path.open("rb") as stream:
        for line in stream:
            m = CANDIDATE.search(line)
            if m and m.group(1).decode() in keys:
                r = REASONS.search(line)
                found[m.group(1).decode()] = [x.strip().strip('"') for x in r.group(1).decode().split(",")] if r else []
    return found


def recall(truth, matched, filtered, output_dir):
    """Per truth allele: tensor (representative / other allele), filtered (with reasons), or why no candidate."""
    status = Counter()
    by_kind = defaultdict(Counter)
    with open(Path(output_dir) / f"{truth.name}.recall.tsv", "w") as out:
        out.write("truth_id\tchrom\tvcf_pos\tkind\tpassed\tin_bed\tstatus\tdetail\n")
        for a in truth.alleles:
            hits = matched.get(a["truth_id"], [])
            if any(rep for rep, _ in hits):
                state, detail = "tensor_representative", ""
            elif hits:
                state, detail = "tensor_non_representative_allele", ""
            elif a["kind"] == "COMPLEX":
                state, detail = "complex_allele", ""
            elif not a["keys"]:
                state, detail = "no_unique_grch38_node_key", ""
            else:
                reasons = sorted({r for k in a["keys"] for r in filtered.get(k, ())})
                state, detail = ("filtered", ",".join(reasons)) if reasons else ("no_candidate", "")
            status[state] += 1
            by_kind[a["kind"]][state] += 1
            if state == "filtered":
                by_kind[a["kind"]]["filtered:" + detail] += 1
            out.write(f"{a['truth_id']}\t{a['chrom']}\t{a['vcf_pos']}\t{a['kind']}\t{a['passed']}\t{a['in_bed']}"
                      f"\t{state}\t{detail}\n")
    return dict(alleles=len(truth.alleles), status=dict(status), by_kind={k: dict(v) for k, v in by_kind.items()})


def label_run(tensors, kinds, reference_path, fasta, somatic_vcf, somatic_bed, germline_vcf, germline_bed, truth_dir,
              recall_dir=None, snv_min_af=None, indel_min_af=None):
    """Build both truth sets, label every merged autosome directory of `kinds`, write recall reports.

    Truth tables (<set>.graph.tsv) depend only on the graph and the VCFs and go to `truth_dir`;
    recall reports depend on the run and go to `recall_dir` (default: `tensors`).
    """
    truth_dir = Path(truth_dir)
    truth_dir.mkdir(parents=True, exist_ok=True)
    recall_dir = Path(recall_dir or tensors)
    recall_dir.mkdir(parents=True, exist_ok=True)
    path = ReferencePath(reference_path)
    locator = Locator(path)
    somatic = TruthSet("somatic", somatic_vcf, somatic_bed, fasta, locator)
    germline = TruthSet("germline", germline_vcf, germline_bed, fasta, locator)
    confident = somatic.bed.intersect(germline.bed)
    for truth in (somatic, germline):
        truth.write_table(truth_dir / f"{truth.name}.graph.tsv")
    provenance = dict(truth={t.name: dict(vcf=str(t.vcf), vcf_sha256=sha256_file(t.vcf), bed=str(t.bed_path),
                                          bed_sha256=sha256_file(t.bed_path), stats=dict(t.stats),
                                          table=str(truth_dir / f"{t.name}.graph.tsv"))
                             for t in (somatic, germline)},
                      confident_region=dict(definition="somatic BED ∩ germline BED", bp=confident.size()),
                      fasta=str(fasta), reference_path=str(path.directory))
    reports, matched = {}, {somatic.name: defaultdict(list), germline.name: defaultdict(list)}
    for kind in kinds:
        reports[kind], found = label_directory(Path(tensors) / kind, somatic, germline, path, confident, provenance, fasta,
                                               snv_min_af, indel_min_af)
        for name, items in found.items():
            for tid, hits in items.items():
                matched[name][tid] += hits
    keys = {k for truth in (somatic, germline) for a in truth.alleles for k in a["keys"]}
    filtered = {}
    for kind in kinds:
        filtered.update(scan_filtered(Path(tensors) / kind, keys))
    recall_report = {t.name: recall(t, matched[t.name], filtered, recall_dir) for t in (somatic, germline)}
    write_json(recall_dir / "truth_recall.json", recall_report)
    return dict(labels={k: dict(totals=r["totals"], reasons=r["reasons"]) for k, r in reports.items()},
                recall=recall_report, truth_stats={t.name: dict(t.stats) for t in (somatic, germline)})
