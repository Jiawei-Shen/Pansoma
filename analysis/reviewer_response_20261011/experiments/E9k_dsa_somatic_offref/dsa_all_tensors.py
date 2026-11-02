"""E9k part 5 (user 2026-10-29): how many of ALL the COLO829BL-DSA somatic SNVs (44,810) can give a tensor on the HPRC
graphs (v1.1 d9, v2.1 d46), in the six COLO829T SNV tensor sets. Part 4 (dsa_only_tensors.py) did this for the 5,500
without a GRCh38 position or on a germline allele; this runs the same steps on every SNV, in DATA = .../E9k/all:

  dsa_all_tensors.py windows   events.tsv + windows.fq for all 44,810 SNVs (dsa_only_tensors.windows, every class)
  (job) giraffe to d9 / d46 + E9j's s3_classify.py (tmp/reviewer_cpu_20261029/E9k/giraffe_all.sh)
  dsa_all_tensors.py labels    a tensor of a set "has" the SNV when its representative allele A1 is
                                (key) the candidate key of the SNV's A / B single-base edit (E9j's rule, as part 4), or
                                (GRCh38) the SNV lifted to GRCh38 (chain, dsa_vs_grch38.py): the tensor's GRCh38 position and
                                ALT (labels' `grch38`, ReferencePath.linear) equal the lifted position and ALT;
                                "site" = a tensor at the same GRCh38 position or graph site with another A1. Labels as in
                                tensors/SNV/chr*_labels.ndjson.
  dsa_all_tensors.py pileup    why many SNVs have no tensor: ALT reads at each SNV's GRCh38 position in the COLO829T
                                GRCh38 BAMs the tensor sets come from (Fiber-seq B004 read once per chromosome, Illumina A001 a fetch per site;
                                bases as pysam count_coverage; MAPQ >= 10 as the tensor builder, and any MAPQ). Positions: the GRCh38 lift (lifted SNVs; a
                                random 2,000 of the truth SNVs as a control) or, for DSA-only SNVs of class A on d9, the
                                GRCh38 position of the graph node (part 4b labels.tsv).
Writes DATA/{labels_all.tsv, pileup.tsv} and here dsa_all_tensors.md / .json, dsa_all_pileup.md / .json.
"""
import csv
import glob
import gzip
import json
import os
import random
import re
import sys
from bisect import bisect_left
from collections import Counter, defaultdict
from multiprocessing import Pool
from pathlib import Path

import dsa_only_tensors as T
import pysam

T.DATA = Path("/scratch/jshen/data/reviewer_experiments_20261013/E9k_dsa_somatic_offref/all")
LIN = re.compile(r'"grch38": \{"chrom": "([^"]+)", "pos0": (\d+), "ref": "([^"]*)", "alt": "([^"]*)"')
CLASSES = [("grch38_same", True, "lifts, same base, in the SMaHT GRCh38 truth"),
           ("grch38_same", False, "lifts, same base, not in the truth"),
           ("grch38_diff_ref", None, "lifts, on a germline allele"),
           ("dsa_only", None, "no GRCh38 position in the chain (DSA only)")]


def scan(args):
    """Tensors of one labels file whose A1 is a wanted key or a wanted lifted GRCh38 SNV, and the GRCh38 / graph sites."""
    st, path, keys, sites, lins, lin_sites = args
    by_key, by_lin, site_key, site_lin = {}, {}, set(), set()
    with open(path) as fh:
        for line in fh:
            i = line.find('"candidate_id": "') + 17
            cid = line[i:line.index('"', i)]
            if ":SNP:" not in cid:
                continue
            lab = None
            if cid in keys:
                r = json.loads(line)
                lab = (r["label_name"], r["reason"])
                by_key[cid] = lab
            node_site = cid.rsplit(":", 2)[0] + ":"
            if node_site in sites:
                site_key.add(node_site)
            m = LIN.search(line)
            if m:
                c, p1, alt = m.group(1), int(m.group(2)) + 1, m.group(4)
                if (c, p1) in lin_sites:
                    site_lin.add((c, p1))
                    if (c, p1, alt) in lins:
                        if lab is None:
                            r = json.loads(line)
                            lab = (r["label_name"], r["reason"])
                        by_lin[(c, p1, alt)] = lab
    return st, by_key, by_lin, site_key, site_lin


def labels():
    snv = {}
    with gzip.open(T.HERE / "dsa_snvs.tsv.gz", "rt") as fh:
        for r in csv.DictReader(fh, delimiter="\t"):
            snv[f"{r['contig']}:{r['pos']}"] = r
    lifted = {}  # event -> (chrom, pos1, alt on GRCh38)
    events = {}
    for r in csv.DictReader(open(T.DATA / "events.tsv"), delimiter="\t"):
        d = snv[r["grch38_records"]]
        events[r["event"]] = (d["cls"], d["in_grch38_truth"] == "True")
        if d["grch38"]:
            c, p, strand = d["grch38"].rsplit(":", 2)
            alt = d["alt"] if strand == "+" else d["alt"].translate(T.COMP)
            lifted[r["event"]] = (c, int(p), alt)
    graph_cls, keys = {}, {}
    for g in ("d9", "d46"):
        rows, k = T.event_keys(g)
        for r in rows:
            graph_cls[(g, r["event"])] = r["cls"].split("_")[0] if r["cls"].startswith("D") else r["cls"]
        keys[g] = k
    lins = set(lifted.values())
    lin_sites = {(c, p) for c, p, _ in lins}
    tasks = []
    for g, sets in T.SETS.items():
        kk = {v[0] for v in keys[g].values()}
        ss = {v[1] for v in keys[g].values()}
        for st in sets:
            for f in sorted(glob.glob(str(T.TENSORS / st / "tensors/SNV/chr*_labels.ndjson"))):
                tasks.append((st, f, kk, ss, lins, lin_sites))
    by_key, by_lin = defaultdict(dict), defaultdict(dict)
    site_key, site_lin = defaultdict(set), defaultdict(set)
    with Pool(int(os.environ.get("SLURM_CPUS_PER_TASK", "4"))) as pool:
        for st, bk, bl, sk, sl in pool.imap_unordered(scan, tasks, chunksize=1):
            by_key[st].update(bk)
            by_lin[st].update(bl)
            site_key[st] |= sk
            site_lin[st] |= sl
    out = []
    for e, (cls, in_truth) in events.items():
        row = dict(event=e, dsa_class=cls, in_truth=in_truth, grch38=":".join(map(str, lifted.get(e, ()))))
        for g, sets in T.SETS.items():
            row[f"{g}_graph_class"] = graph_cls.get((g, e), "")
            key, site = keys[g].get(e, ("", ""))
            row[f"{g}_key"] = key
            lin = lifted.get(e)
            for st in sets:
                lab = by_key[st].get(key) if key else None
                how = "key" if lab else ""
                if lin and lin in by_lin[st]:
                    lab = lab or by_lin[st][lin]
                    how = "key+grch38" if how else "grch38"
                if lab:
                    row[st] = f"exact|{how}|{lab[0]}|{lab[1]}"
                elif (key and site in site_key[st]) or (lin and lin[:2] in site_lin[st]):
                    row[st] = "site"
                else:
                    row[st] = "none"
        out.append(row)
    with open(T.DATA / "labels_all.tsv", "w") as w:
        dw = csv.DictWriter(w, fieldnames=list(out[0]), delimiter="\t")
        dw.writeheader()
        dw.writerows(out)
    allsets = [st for v in T.SETS.values() for st in v]
    L = ["# E9k part 5: tensors on the HPRC graphs for all COLO829BL-DSA somatic SNVs (dsa_all_tensors.py)", "",
         f"{len(out):,} DSA somatic SNVs. A tensor has the SNV when its representative allele A1 is the SNV (by the candidate key",
         "of the DSA window's edit on the graph, or by the SNV's GRCh38 lift); \"site\" = a tensor there with another A1.", "",
         "## SNVs with a tensor (exact A1)", "",
         "| DSA class | SNVs | " + " | ".join(allsets) + " | any v1.1 set | any v2.1 set | any of the six |",
         "|---|---:|" + "---:|" * (len(allsets) + 3)]
    summary = {}

    def exact(r, st):
        return r[st].startswith("exact")

    groups = [(lab, [r for r in out if r["dsa_class"] == c and (t is None or r["in_truth"] == t)]) for c, t, lab in CLASSES]
    groups.append(("all", out))
    for lab, sel in groups:
        n = len(sel)
        cells = [sum(exact(r, st) for r in sel) for st in allsets]
        a9 = sum(any(exact(r, st) for st in T.SETS["d9"]) for r in sel)
        a46 = sum(any(exact(r, st) for st in T.SETS["d46"]) for r in sel)
        a = sum(any(exact(r, st) for st in allsets) for r in sel)
        L.append(f"| {lab} | {n:,} | " + " | ".join(f"{c:,} ({100 * c / n:.1f}%)" for c in cells)
                 + f" | {a9:,} ({100 * a9 / n:.1f}%) | {a46:,} ({100 * a46 / n:.1f}%) | {a:,} ({100 * a / n:.1f}%) |")
        summary[lab] = dict(snvs=n, **dict(zip(allsets, cells)), any_d9=a9, any_d46=a46, any=a)
    L += ["", "## Where each SNV lands on the graph (DSA window, E9j classes) and whether a set has its tensor", ""]
    for g, sets in T.SETS.items():
        L += [f"### {g}", "", "| DSA class | graph class | SNVs | " + " | ".join(f"{st}: exact / site" for st in sets) + " |",
              "|---|---|---:|" + "---:|" * len(sets)]
        for lab, sel in groups[:-1]:
            for gc in ("A", "B", "C", "D"):
                s2 = [r for r in sel if r[f"{g}_graph_class"] == gc]
                if not s2:
                    continue
                cells = [f"{sum(exact(r, st) for r in s2):,} / {sum(r[st] == 'site' for r in s2):,}" for st in sets]
                L.append(f"| {lab} | {gc} | {len(s2):,} | " + " | ".join(cells) + " |")
        L.append("")
    L += ["## Labels of the exact tensors (all SNVs)", ""]
    for st in allsets:
        c = Counter("|".join(r[st].split("|")[2:]) for r in out if exact(r, st))
        how = Counter(r[st].split("|")[1] for r in out if exact(r, st))
        L.append(f"- {st}: " + "; ".join(f"{k.replace('|', ' / ')} {v:,}" for k, v in c.most_common())
                 + " (found by " + ", ".join(f"{k} {v:,}" for k, v in how.most_common()) + ")")
        summary[f"labels|{st}"] = dict(c)
    (T.HERE / "dsa_all_tensors.md").write_text("\n".join(L) + "\n")
    (T.HERE / "dsa_all_tensors.json").write_text(json.dumps(summary, indent=1) + "\n")
    print("\n".join(L))


BAMS = {"fiberseq": "/scratch/jshen/data/COLO829T/COLO829T_fiberseq/SMHTCOLO829T-X-X-M45-B004-uwsc-MERGED-pbmm2_1.13.0_GRCh38.aligned.sorted.bam",
        "illumina": "/scratch/qfu/COLO829T_WGS/SMHTCOLO829T-X-X-M45-A001-uwsc-SMAFIZO17Q9I-sentieon_bwamem_202308.01_GRCh38.aligned.sorted.bam"}
BASES = "ACGT"


STREAM = {"fiberseq": True, "illumina": False}  # long reads: one pass over each chromosome (a fetch per site re-reads
PARTS = T.DATA / "pileup_parts"                  # the same blocks, 757 GB in 2 h for job 461142); short reads: a fetch per site


def qpos(read, p0):
    """Index in read.query_sequence of reference position p0, or None (deletion, skip, outside the alignment)."""
    q, r = 0, read.reference_start
    for op, n in read.cigartuples:
        if op in (0, 7, 8):
            if p0 < r + n:
                return q + p0 - r
            q, r = q + n, r + n
        elif op in (1, 4):
            q += n
        elif op in (2, 3):
            if p0 < r + n:
                return None
            r += n
    return None


def count(args):
    """PARTS/{platform}.{chrom}.tsv: event, depth and ALT reads at MAPQ >= 10, depth and ALT reads at any MAPQ. Reads
    without the unmapped / secondary / QC-fail / duplicate flags and bases A/C/G/T (as count_coverage), no base-quality
    threshold. Skipped when the part exists (resume)."""
    platform, chrom, sites = args
    out = PARTS / f"{platform}.{chrom}.tsv"
    if out.exists():
        return out
    at = defaultdict(list)
    for e, p1, alt in sites:
        at[p1 - 1].append((e, alt))
    pos = sorted(at)
    n = {e: [0, 0, 0, 0] for e, _, _ in sites}
    with pysam.AlignmentFile(BAMS[platform]) as bam:
        if chrom in bam.references:
            if STREAM[platform]:
                reads = ((r, None) for r in bam.fetch(chrom, pos[0], pos[-1] + 1))
            else:
                reads = ((r, [p]) for p in pos for r in bam.fetch(chrom, p, p + 1))
            for r, ps in reads:
                if r.flag & 0x704 or r.query_sequence is None:
                    continue
                if ps is None:
                    ps = pos[bisect_left(pos, r.reference_start):bisect_left(pos, r.reference_end)]
                seq, q10 = r.query_sequence, r.mapping_quality >= 10
                for p in ps:
                    i = qpos(r, p)
                    if i is None or seq[i] not in BASES:
                        continue
                    for e, alt in at[p]:
                        c, hit = n[e], seq[i] == alt
                        c[2] += 1
                        c[3] += hit
                        if q10:
                            c[0] += 1
                            c[1] += hit
    tmp = out.with_suffix(".tmp")
    tmp.write_text("".join(f"{e}\t" + "\t".join(map(str, c)) + "\n" for e, c in n.items()))
    tmp.rename(out)
    return out


def pileup():
    rows = {r["event"]: r for r in csv.DictReader(open(T.DATA / "labels_all.tsv"), delimiter="\t")}
    part4 = {}  # DSA-only A events on d9: GRCh38 position of the node (part 4b), by DSA site
    p4ev = {r["event"]: r["grch38_records"] for r in csv.DictReader(open(T.DATA.parent / "dsa_only/events.tsv"), delimiter="\t")}
    for r in csv.DictReader(open(T.DATA.parent / "dsa_only/labels.tsv"), delimiter="\t"):
        if r["graph"] == "d9" and r["graph_class"] == "A" and r["grch38"]:
            c, p, ch = r["grch38"].split(":")
            part4[p4ev[r["event"]]] = (c, int(p), ch.split(">")[1])
    site = {r["event"]: r["grch38_records"] for r in csv.DictReader(open(T.DATA / "events.tsv"), delimiter="\t")}
    random.seed(1)
    truth = sorted(e for e, r in rows.items() if r["dsa_class"] == "grch38_same" and r["in_truth"] == "True")
    keep = set(random.sample(truth, 2000))
    todo = []
    for e, r in rows.items():
        if r["grch38"] and (r["in_truth"] == "False" or e in keep):
            c, p, alt = r["grch38"].split(":")
            todo.append((e, c, int(p), alt, "lift"))
        elif r["dsa_class"] == "dsa_only" and site[e] in part4:
            todo.append((e, *part4[site[e]], "d9 node"))
    todo.sort(key=lambda x: (x[1], x[2]))
    by_chrom = defaultdict(list)
    for e, c, p, alt, how in todo:
        by_chrom[c].append((e, p, alt))
    tasks = sorted(((pf, c, v) for pf in BAMS for c, v in by_chrom.items()), key=lambda t: (not STREAM[t[0]], -len(t[2])))
    PARTS.mkdir(exist_ok=True)
    res = defaultdict(dict)
    with Pool(int(os.environ.get("SLURM_CPUS_PER_TASK", "4"))) as pool:
        for out in pool.imap_unordered(count, tasks):
            print("done", out.name, flush=True)
    for pf in BAMS:
        for c in by_chrom:
            for line in open(PARTS / f"{pf}.{c}.tsv"):
                e, *x = line.split()
                res[e][pf] = tuple(map(int, x))
    cols = ["event", "dsa_class", "in_truth", "pos_from", "grch38"] + [f"{pf}_{k}" for pf in BAMS for k in ("dp_q10", "alt_q10", "dp_all", "alt_all")]
    with open(T.DATA / "pileup.tsv", "w") as w:
        w.write("\t".join(cols) + "\n")
        for e, c, p, alt, how in todo:
            r = rows[e]
            w.write("\t".join([e, r["dsa_class"], r["in_truth"], how, f"{c}:{p}:{alt}"]
                               + [str(x) for pf in BAMS for x in res[e][pf]]) + "\n")
    sets = {"fiberseq": ["COLO829T_fiberseq", "COLO829T_fiberseq_d46"], "illumina": ["COLO829T_Illumina", "COLO829T_Illumina_d46"]}
    auto = {f"chr{i}" for i in range(1, 23)}
    truth_xy = set()  # the SMaHT truth's chrX / chrY SNVs (in_truth above is autosomal, as the chr1-22 tensor sets)
    for v in pysam.VariantFile(str(T.TRUTH)):
        if v.chrom in ("chrX", "chrY") and len(v.ref) == 1 and v.alts and len(v.alts[0]) == 1:
            truth_xy.add(f"{v.chrom}:{v.pos}:{v.alts[0].upper()}")
    g38ref = {}  # GRCh38 base at the lift of each germline-allele SNV (dsa_vs_grch38.py)
    for r in csv.DictReader(gzip.open(T.HERE / "dsa_snvs.tsv.gz", "rt"), delimiter="\t"):
        if r["cls"] == "grch38_diff_ref":
            g38ref[":".join(r["grch38"].split(":")[:2])] = r["grch38_ref"]
    med = {pf: sorted(res[e][pf][0] for e, *_ in todo if rows[e]["in_truth"] == "True")[len(keep) // 2] for pf in BAMS}

    def reads(dq, aq, da, aa):
        if aq >= 3:
            return "ALT >= 3 reads and >= 8% at MAPQ >= 10" if aq >= 0.08 * dq else "ALT >= 3 reads but < 8% at MAPQ >= 10"
        if aq >= 1:
            return "ALT 1-2 reads at MAPQ >= 10"
        if aa >= 1:
            return "ALT only in reads with MAPQ < 10"
        return "no ALT read" if da else "no reads"
    order = ["ALT >= 3 reads and >= 8% at MAPQ >= 10", "ALT >= 3 reads but < 8% at MAPQ >= 10", "ALT 1-2 reads at MAPQ >= 10",
             "ALT only in reads with MAPQ < 10", "no ALT read", "no reads"]
    ch = lambda c: "chr1-22" if c in auto else c if c in ("chrX", "chrY") else "other"
    groups = [("lifts, same base, in the SMaHT truth (chr1-22; random 2,000; control)",
               lambda r, c, p, a: r["dsa_class"] == "grch38_same" and r["in_truth"] == "True"),
              ("lifts, same base, not in the truth: chr1-22",
               lambda r, c, p, a: r["dsa_class"] == "grch38_same" and r["in_truth"] == "False" and c in auto),
              ("lifts, same base, not in the truth: chrX / chrY (outside the chr1-22 tensor sets)",
               lambda r, c, p, a: r["dsa_class"] == "grch38_same" and r["in_truth"] == "False" and c in ("chrX", "chrY")),
              ("lifts, same base, not in the truth: other GRCh38 contigs (outside the tensor sets)",
               lambda r, c, p, a: r["dsa_class"] == "grch38_same" and r["in_truth"] == "False" and ch(c) == "other"),
              ("lifts, on a germline allele: ALT = the GRCh38 base (reads with it are reference reads)",
               lambda r, c, p, a: r["dsa_class"] == "grch38_diff_ref" and g38ref[f"{c}:{p}"] == a),
              ("lifts, on a germline allele: ALT another base",
               lambda r, c, p, a: r["dsa_class"] == "grch38_diff_ref" and g38ref[f"{c}:{p}"] != a),
              ("DSA only, class A on d9, node on chr1-22 (GRCh38 position of the graph node)",
               lambda r, c, p, a: r["dsa_class"] == "dsa_only" and c in auto),
              ("DSA only, class A on d9, node elsewhere", lambda r, c, p, a: r["dsa_class"] == "dsa_only" and c not in auto)]
    lifted = [(e, c, p, alt) for e, c, p, alt, how in todo if rows[e]["dsa_class"] == "grch38_same" and rows[e]["in_truth"] == "False"]
    xy = [x for x in lifted if x[1] in ("chrX", "chrY")]
    L = ["# E9k part 5b: ALT reads at the DSA somatic SNVs in the COLO829T GRCh38 BAMs (dsa_all_tensors.py pileup)", "",
         "Reads at each SNV's GRCh38 position in the BAM of the reads the tensor set was built from (same sample):",
         f"Fiber-seq `{Path(BAMS['fiberseq']).name}`, Illumina `{Path(BAMS['illumina']).name}`. Bases as pysam count_coverage, no base",
         "quality threshold; MAPQ >= 10 = the tensor builder's min_mapq, 8% = its SNV AF threshold. Tensor = a tensor with this",
         "SNV as A1 in the v1.1 or v2.1 set of that platform. Depth = reads at MAPQ >= 10; the median at the truth controls is",
         f"{med['fiberseq']} (Fiber-seq) and {med['illumina']} (Illumina).", "",
         f"\"Not in the truth\" is against the truth's chr1-22 SNVs (the tensor sets are chr1-22): {len(xy):,} of the {len(lifted):,} lifted",
         f"SNVs not in it lift to chrX / chrY, and {sum(f'{c}:{p}:{a}' in truth_xy for e, c, p, a in xy):,} of these are in the SMaHT truth's chrX / chrY",
         "part.", ""]
    summary = dict(median_depth=med, lifted_not_in_truth=len(lifted), lifted_xy=len(xy),
                   lifted_xy_in_truth=sum(f"{c}:{p}:{a}" in truth_xy for e, c, p, a in xy))
    for pf in BAMS:
        L += [f"## {pf}", "", "| SNVs | reads at the SNV | SNVs | depth >= 2x median | with a tensor (v1.1 or v2.1 set) |",
              "|---|---|---:|---:|---:|"]
        for lab, f in groups:
            sel = [e for e, c, p, alt, how in todo if f(rows[e], c, p, alt)]
            if not sel:
                continue
            c, d, t = Counter(), Counter(), Counter()
            for e in sel:
                k = reads(*res[e][pf])
                c[k] += 1
                d[k] += res[e][pf][0] >= 2 * med[pf]
                t[k] += any(rows[e][st].startswith("exact") for st in sets[pf])
            for k in order:
                if c[k]:
                    L.append(f"| {lab} ({len(sel):,}) | {k} | {c[k]:,} ({100 * c[k] / len(sel):.1f}%) | {d[k]:,} | {t[k]:,} ({100 * t[k] / c[k]:.0f}%) |")
                    summary[f"{pf}|{lab}|{k}"] = dict(snvs=c[k], depth_2x=d[k], tensor=t[k])
            n_t = sum(t.values())
            L.append(f"| {lab} ({len(sel):,}) | all | {len(sel):,} | {sum(d.values()):,} | {n_t:,} ({100 * n_t / len(sel):.0f}%) |")
        L.append("")
    (T.HERE / "dsa_all_pileup.json").write_text(json.dumps(summary, indent=1) + "\n")
    print("\n".join(L))


if __name__ == "__main__":
    if sys.argv[1] == "windows":
        T.windows(classes=("grch38_same", "grch38_diff_ref", "dsa_only"))
    elif sys.argv[1] == "pileup":
        pileup()
    else:
        labels()
