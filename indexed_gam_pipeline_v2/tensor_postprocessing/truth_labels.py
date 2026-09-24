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
   at a node junction gets a key on both nodes; a deletion must fit inside one node; nodes
   the reference visits more than once give no key.
4. Labels. A tensor's representative allele (its candidate_id) is looked up:

       1 somatic    representative allele is a somatic truth allele
       2 germline   representative allele is a germline truth allele (FILTER PASS/., inside the germline BED)
       0 non        on a unique GRCh38 node, inside somatic BED ∩ germline BED, no truth allele at the site,
                    and no truth allele (any set, any filter) within NEAR_BP
      -1 ignore     everything else, with a reason (outside the confident region, not on a unique GRCh38 node,
                    near a truth allele with a different allele, filtered germline record, or the truth allele
                    is a non-representative allele of the site)

Outputs next to the merged shards (per chromosome, same order as <chrom>_variant_summary.ndjson):
<chrom>_shard_NNNNN_labels.npy (int8) and <chrom>_labels.ndjson; labels.manifest.json;
<truth dir>/<set>.graph.tsv (every truth allele with its keys) and <set>.recall.tsv / recall.json.
"""
from collections import Counter, defaultdict
from datetime import datetime, timezone
import json
from pathlib import Path
import re

import numpy as np

from indexed_gam_pipeline_v2.common import read_json, sha256_file, write_json
from indexed_gam_pipeline_v2.tensor_postprocessing.chr_index import AUTOSOMES
from indexed_gam_pipeline_v2.tensor_postprocessing.reference_path import ReferencePath, rc

VERSION = "truth-labels-v1"
LABELS = {"ignore": -1, "non": 0, "somatic": 1, "germline": 2}
NEAR_BP = 10
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
        v5 (single-node deletions only) simply never matches the suffixed keys.
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
                                             placements=len(places), keys=keys))
                    spans[record.chrom].append((lo, max(hi, lo + 1) if kind != "INS" else lo + 1, tid))
                    for key in keys:
                        self.by_key[key].append(tid)
                    self.stats["alleles_with_keys" if keys else "alleles_without_keys"] += 1
        genome_file.close()
        self.spans = {}
        for chrom, items in spans.items():
            items.sort()
            lo = np.array([i[0] for i in items], dtype=np.int64)
            hi = np.array([i[1] for i in items], dtype=np.int64)
            self.spans[chrom] = (lo, np.maximum.accumulate(hi))

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


# --- labels -------------------------------------------------------------------------

def linear_interval(lin, kind):
    if kind == "INS":
        return lin["pos0"] - 1, lin["pos0"] + 1
    return lin["pos0"], lin["pos0"] + max(1, len(lin["ref"]))


def classify(record, somatic, germline, path, confident):
    """(label value, label name, reason, details) of one merged summary record."""
    representative = record["candidate_id"]
    alleles = record.get("alleles") or [record]
    hits = {}
    for truth in (somatic, germline):
        hits[truth.name] = [(a["candidate_id"], tid) for a in alleles for tid in truth.by_key.get(a["candidate_id"], ())]
    details = {truth.name: [truth.summary(tid, cid, cid == representative) for cid, tid in hits[truth.name]]
               for truth in (somatic, germline)}
    lin = path.linear(record["node_id"], record["start"], record["ref"], record["alt"], record["event_type"],
                      record.get("path"))
    details["grch38"] = lin
    rep_somatic = [tid for cid, tid in hits[somatic.name] if cid == representative]
    rep_germline = [tid for cid, tid in hits[germline.name] if cid == representative]
    if rep_somatic:
        return LABELS["somatic"], "somatic", "representative_allele_is_somatic_truth", details
    if rep_germline:
        good = [t for t in rep_germline if germline.alleles[t]["passed"] and germline.alleles[t]["in_bed"]]
        if good:
            return LABELS["germline"], "germline", "representative_allele_is_germline_truth", details
        return LABELS["ignore"], "ignore", "germline_truth_filtered_or_outside_bed", details
    if hits[somatic.name] or hits[germline.name]:
        return LABELS["ignore"], "ignore", "truth_matches_non_representative_allele", details
    if lin is None:
        return LABELS["ignore"], "ignore", "not_on_unique_grch38_node", details
    s, e = linear_interval(lin, record["event_type"])
    if not confident.contains(lin["chrom"], s, e):
        return LABELS["ignore"], "ignore", "outside_confident_region", details
    if somatic.near(lin["chrom"], s, e) or germline.near(lin["chrom"], s, e):
        return LABELS["ignore"], "ignore", "near_truth_allele_mismatch", details
    return LABELS["non"], "non", "confident_no_truth_allele", details


def label_directory(directory, somatic, germline, path, confident, provenance):
    """Write <chrom>_shard_*_labels.npy, <chrom>_labels.ndjson and labels.manifest.json for one merged directory."""
    directory = Path(directory)
    manifest = read_json(directory / "manifest.json")
    if manifest.get("layout") != "chromosome-shards-v1":
        raise ValueError(f"{directory} is not a merged per-chromosome directory")
    counts, reasons, written = {}, Counter(), []
    matched = {somatic.name: defaultdict(list), germline.name: defaultdict(list)}
    for chrom, info in manifest["chromosomes"].items():
        counts[chrom] = Counter()
        sizes = {s["file"]: s["tensors"] for s in info["shards"]}
        arrays = {f: np.full(n, -128, dtype=np.int8) for f, n in sizes.items()}
        target = directory / f"{chrom}_labels.ndjson"
        with (directory / info["summary"]).open() as source, target.with_name(target.name + ".tmp").open("w") as out:
            for line in source:
                record = json.loads(line)
                value, name, reason, details = classify(record, somatic, germline, path, confident)
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
    report = dict(version=VERSION, created=datetime.now(timezone.utc).isoformat(), labels=LABELS,
                  near_bp=NEAR_BP, tensors=sum(sum(c.values()) for c in counts.values()),
                  counts={c: dict(v) for c, v in counts.items()},
                  totals=dict(sum(counts.values(), Counter())), reasons=dict(reasons), **provenance)
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
              recall_dir=None):
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
        reports[kind], found = label_directory(Path(tensors) / kind, somatic, germline, path, confident, provenance)
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
