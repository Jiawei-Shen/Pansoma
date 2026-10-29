"""E9k part 4 (user 2026-10-29): can Pansoma build a tensor for the DSA somatic SNVs that have no GRCh38 position (dsa_only)
or sit on a germline allele of the normal (grch38_diff_ref)? Placement on the graphs as E9j (COLO829T), from the DSA directly.

  dsa_only_tensors.py windows   events.tsv + windows.fq in DATA: per SNV of dsa_snvs.tsv.gz (dsa_vs_grch38.py) with class
                                dsa_only or grch38_diff_ref, the DSA window of its haplotype (haplotype1-* = hapX,
                                haplotype2-* = hapY), 100-bp flanks: N = the DSA sequence, T = the same with the ALT
  (job) giraffe to d9 / d46 with E9j's settings (tmp/reviewer_cpu_20261029/E9k/giraffe.sh), then E9j's s3_classify.py
        unchanged: A = the somatic edit on a GRCh38 node, B = on an off-reference node, C = the allele is a graph branch,
        D = unmapped / MAPQ < 20 / same alignment
  dsa_only_tensors.py tensors   candidate keys (node:start:SNP:REF>ALT, node-forward) of every A and B event (E9j's key
                                rule, applied to A too) against the candidate ids of the six COLO829T SNV tensor sets
                                (index_cache/*.candidates.txt = the representative allele A1 of each tensor): exact A1, or
                                the same site (node:start) with another A1 (the allele may be a non-representative allele)
  dsa_only_tensors.py labels    GRCh38 position of every A key (ReferencePath.linear, as the labels), whether it is in the
                                SMaHT truth and in the label region (dip.bed), and the label of every exact tensor
                                (tensors/SNV/chr*_labels.ndjson); baseline: SMaHT truth SNVs in dip.bed that are an A1
Writes DATA/{events.tsv, windows.fq, classes_<g>.tsv, keys_tensors.tsv, labels.tsv} and here dsa_only_tensors.md / .json,
dsa_only_labels.md / .json.
"""
import bisect
import csv
import glob
import gzip
import json
import os
import sys
from collections import Counter, defaultdict
from pathlib import Path

import pysam

HERE = Path(__file__).resolve().parent
DATA = Path("/scratch/jshen/data/reviewer_experiments_20261013/E9k_dsa_somatic_offref/dsa_only")
HAPS = {"haplotype1": ("hapX", "/scratch/qfu/COLO829BL_DSA/SMHTCOLO829BL-X-X-M45-XX-uwsc-SMASFTH1B2J6-verkko_2.1_hapX.fa"),
        "haplotype2": ("hapY", "/scratch/qfu/COLO829BL_DSA/SMHTCOLO829BL-X-X-M45-XX-uwsc-SMASFK4V33V8-verkko_2.1_hapY.fa")}
FLANK = 100
CACHE = Path("/scratch/jshen/data/pansoma_net_v2_runs/index_cache")
SETS = {"d9": ["COLO829T_Illumina", "COLO829T_fiberseq", "COLO829T_ONT"],
        "d46": ["COLO829T_Illumina_d46", "COLO829T_fiberseq_d46", "COLO829T_ONT_d46"]}
COMP = str.maketrans("ACGTN", "TGCAN")
E9J = HERE.parent / "E9j_offref_colo829t"
TENSORS = Path("/scratch/jshen/data/pansoma_v2_tensors")
TRUTH = TENSORS / "COLO829T_truth/COLO829T_somatic_snv_indel.vcf.gz"
LABEL_BED = "/scratch/qfu/COLO829BL_DSA/dipcall_hg38/dipcall_hg38.dip.bed"  # label 0 region (COLO829T_truth/README.txt)
AUTOSOMES = {f"chr{i}" for i in range(1, 23)}


def windows():
    DATA.mkdir(parents=True, exist_ok=True)
    fas = {k: pysam.FastaFile(v[1]) for k, v in HAPS.items()}
    n = 0
    with gzip.open(HERE / "dsa_snvs.tsv.gz", "rt") as fh, open(DATA / "events.tsv", "w") as ev, open(DATA / "windows.fq", "w") as fq:
        ev.write("event\tchrom\tgrch38_pos\tref\talt\thap\thow\tcontig\tpos0\tgrch38_records\n")
        for r in csv.DictReader(fh, delimiter="\t"):
            if r["cls"] not in ("dsa_only", "grch38_diff_ref"):
                continue
            pref = r["contig"].split("-")[0]
            hap, fa = HAPS[pref][0], fas[pref]
            p0 = int(r["pos"]) - 1
            lo = max(0, p0 - FLANK)
            seq = fa.fetch(r["contig"], lo, p0 + FLANK + 1).upper()
            if seq[p0 - lo] != r["ref"]:
                raise SystemExit(f"{r['contig']}:{r['pos']}: DSA base {seq[p0 - lo]} != VCF REF {r['ref']}")
            tseq = seq[:p0 - lo] + r["alt"] + seq[p0 - lo + 1:]
            eid = f"k{n:05d}"
            n += 1
            ev.write(f"{eid}\t{r['cls']}\t{r['grch38'] or '.'}\t{r['ref']}\t{r['alt']}\t{hap}\tdsa\t{r['contig']}\t{p0}\t{r['contig']}:{r['pos']}\n")
            for which, s in (("N", seq), ("T", tseq)):
                fq.write(f"@{eid}_{which} {r['contig']}:{p0} flank={p0 - lo}\n{s}\n+\n{'I' * len(s)}\n")
    print(f"{n:,} events -> {DATA}")


def node_seqs(graph, nodes):
    sys.path.insert(0, str(E9J))
    import s3_classify as S
    from indexed_gam_pipeline_v4.graph_index import GraphIndex
    with GraphIndex(S.GRAPHS[graph][1]) as gi:
        return {k: v["sequence"] for k, v in gi.get_nodes(nodes).items()}


def event_keys(g):
    """classes rows of graph g, {event: (candidate key, node:start:)} of the A and B single-base edits."""
    rows = list(csv.DictReader(open(DATA / f"classes_{g}.tsv"), delimiter="\t"))
    ab = [r for r in rows if r["cls"] in ("A", "B")]
    seqs = node_seqs(g, {int(r["edit_node"]) for r in ab})
    keys = {}
    for r in ab:
        node, rev, off = int(r["edit_node"]), r["edit_rev"] == "True", int(r["edit_offset"])
        fl_tl, seq = r["edit"].split(":", 1)
        fl, tl = map(int, fl_tl.split(">"))
        if not (fl == tl == 1):
            continue  # not a single-base substitution on the graph (complex or indel edit)
        s, n = seqs[node], len(seqs[node])
        start = n - off - 1 if rev else off
        alt = seq.translate(COMP) if rev else seq
        keys[r["event"]] = (f"{node}:{start}:SNP:{s[start]}>{alt}", f"{node}:{start}:")
    return rows, keys


def tensors():
    out, summary = [], {}
    L = ["# E9k part 4: tensors for the DSA somatic SNVs without a GRCh38 position (dsa_only_tensors.py)", "",
         "Events: DSA somatic SNVs of class dsa_only (no GRCh38 position in the chain) and grch38_diff_ref (on a germline allele),",
         "placed on each graph with E9j's method (giraffe of the DSA windows, s3_classify.py). Tensor = a tensor of the set whose",
         "representative allele is the event's key (exact), or a tensor at the same site with another representative allele.", ""]
    for g in ("d9", "d46"):
        rows, keys = event_keys(g)
        cls = Counter((r["chrom"], r["cls"].split("_")[0] if r["cls"].startswith("D") else r["cls"]) for r in rows)
        L += [f"## {g}", "", "| event class | events | A (GRCh38 node) | B (off-reference node) | C (graph branch) | D (unmapped / MAPQ < 20 / same) |",
              "|---|---:|---:|---:|---:|---:|"]
        for c in ("dsa_only", "grch38_diff_ref"):
            tot = sum(v for (a, b), v in cls.items() if a == c)
            L.append(f"| {c} | {tot:,} | {cls[(c, 'A')]:,} | {cls[(c, 'B')]:,} | {cls[(c, 'C')]:,} | {cls[(c, 'D')]:,} |")
        L += ["", "| set | A or B SNP keys | tensor, exact A1 | tensor at the site, other A1 | of them B (off-reference) exact |", "|---|---:|---:|---:|---:|"]
        ev_cls = {r["event"]: (r["chrom"], r["cls"]) for r in rows}
        for st in SETS[g]:
            f = glob.glob(str(CACHE / f"{st}.tensors.SNV.*.candidates.txt"))
            cands = set(open(f[0]).read().split()) if f else set()
            sites = {c.rsplit(":", 2)[0] + ":" for c in cands}  # node:start:
            exact = {e for e, (k, s) in keys.items() if k in cands}
            site = {e for e, (k, s) in keys.items() if k not in cands and s in sites}
            b_exact = {e for e in exact if ev_cls[e][1] == "B"}
            summary[f"{g}|{st}"] = dict(keys=len(keys), exact=len(exact), site=len(site), b_exact=len(b_exact))
            L.append(f"| {st} | {len(keys):,} | {len(exact):,} | {len(site):,} | {len(b_exact):,} |")
            for e, (k, s) in keys.items():
                out.append(dict(graph=g, set=st, event=e, event_class=ev_cls[e][0], graph_class=ev_cls[e][1], key=k,
                                tensor="exact" if e in exact else ("site" if e in site else "none")))
        L.append("")
    with open(DATA / "keys_tensors.tsv", "w") as w:
        dw = csv.DictWriter(w, fieldnames=list(out[0]), delimiter="\t")
        dw.writeheader()
        dw.writerows(out)
    (HERE / "dsa_only_tensors.md").write_text("\n".join(L) + "\n")
    (HERE / "dsa_only_tensors.json").write_text(json.dumps(summary, indent=1) + "\n")
    print("\n".join(L))


def scan_labels(args):
    """(set, chrom, {key: (label_name, reason)} of the wanted keys, {(chrom, pos1, alt)} of the A1s that are SMaHT truth)."""
    st, path, wanted = args
    hit, truth_a1 = {}, set()
    with open(path) as fh:
        for line in fh:
            i = line.find('"candidate_id": "') + 17
            cid = line[i:line.index('"', i)]
            if cid in wanted or '"representative_allele_is_somatic_truth"' in line:
                r = json.loads(line)
                if cid in wanted:
                    hit[cid] = (r["label_name"], r["reason"])
                if r["reason"] == "representative_allele_is_somatic_truth":
                    lin = r["grch38"]
                    truth_a1.add((lin["chrom"], lin["pos0"] + 1, lin["alt"]))
    return st, Path(path).name.split("_")[0], hit, truth_a1


def in_bed(bed, chrom, pos0):
    v = bed.get(chrom)
    if not v:
        return False
    i = bisect.bisect_right(v[0], pos0) - 1
    return i >= 0 and pos0 < v[1][i]


def labels():
    from multiprocessing import Pool
    sys.path.insert(0, str(E9J))
    import s3_classify as S
    from indexed_gam_pipeline_v4.tensor_postprocessing.reference_path import ReferencePath
    truth = set()
    with pysam.VariantFile(str(TRUTH)) as vcf:
        for r in vcf:
            if len(r.ref) == 1 and r.alts and len(r.alts[0]) == 1:
                truth.add((r.chrom, r.pos, r.alts[0].upper()))
    raw = defaultdict(list)
    for line in open(LABEL_BED):
        c, a, b = line.split()[:3]
        raw[c].append((int(a), int(b)))
    bed = {c: ([a for a, _ in sorted(v)], [b for _, b in sorted(v)]) for c, v in raw.items()}
    truth_in_bed = {t for t in truth if t[0] in AUTOSOMES and in_bed(bed, t[0], t[1] - 1)}
    ev = {}  # (g, event) -> dict
    keys = {}
    for g in ("d9", "d46"):
        rows, k = event_keys(g)
        rp = ReferencePath(S.GRAPHS[g][0])
        for r in rows:
            if r["event"] not in k:
                continue
            key = k[r["event"]][0]
            node, start, _, change = key.split(":")
            ref, alt = change.split(">")
            d = dict(graph=g, event=r["event"], event_class=r["chrom"], graph_class=r["cls"], key=key,
                     grch38="", where="off-reference node", in_truth=False)
            if r["cls"] == "A":
                lin = rp.linear(int(node), int(start), ref, alt, "SNP")
                c, p1 = lin["chrom"], lin["pos0"] + 1
                d["grch38"] = f"{c}:{p1}:{lin['ref']}>{lin['alt']}"
                d["in_truth"] = (c, p1, lin["alt"]) in truth
                d["where"] = ("chr1-22, label region" if in_bed(bed, c, p1 - 1) else "chr1-22, outside the label region") \
                    if c in AUTOSOMES else ("chrX / chrY" if c in ("chrX", "chrY") else "other GRCh38 contig (Un / random / M)")
            ev[(g, r["event"])] = d
            keys[(g, key)] = r["event"]
    tasks = []
    for g, sets in SETS.items():
        wanted = {k for (gg, k) in keys if gg == g}
        for st in sets:
            for f in sorted(glob.glob(str(TENSORS / st / "tensors/SNV/chr*_labels.ndjson"))):
                tasks.append((st, f, wanted))
    hits, truth_a1 = defaultdict(dict), defaultdict(set)
    with Pool(int(os.environ.get("SLURM_CPUS_PER_TASK", "4"))) as pool:
        for st, chrom, hit, ta in pool.imap_unordered(scan_labels, tasks):
            hits[st].update(hit)
            truth_a1[st] |= ta
    out = []
    for g, sets in SETS.items():
        for (gg, e), d in ev.items():
            if gg != g:
                continue
            row = dict(d)
            for st in sets:
                row[st] = "|".join(hits[st].get(d["key"], ("none", "")))
            out.append(row)
    with open(DATA / "labels.tsv", "w") as w:
        cols = list(dict.fromkeys(c for r in out for c in r))
        dw = csv.DictWriter(w, fieldnames=cols, delimiter="\t", restval="")
        dw.writeheader()
        dw.writerows(out)
    wheres = ["chr1-22, label region", "chr1-22, outside the label region", "chrX / chrY", "other GRCh38 contig (Un / random / M)",
              "off-reference node"]
    L = ["# E9k part 4b: labels of the tensors of the DSA-only / germline-allele somatic SNVs (dsa_only_tensors.py labels)", "",
         "A and B events of dsa_only_tensors.md with a single-base key. GRCh38 position of an A key: ReferencePath.linear, as the",
         f"labels. Label region: dip.bed (somatic label 0 region, COLO829T_truth/README.txt). Truth: SMaHT COLO829T ({len(truth):,} SNVs).", ""]
    summary = {}
    for g, sets in SETS.items():
        rows = [r for r in out if r["graph"] == g]
        L += [f"## {g}", "", "| where the edit lands | event class | events | in the SMaHT truth | " + " | ".join(f"{st}: tensor (somatic / non / germline / ignore)" for st in sets) + " |",
              "|---|---|---:|---:|" + "---:|" * len(sets)]
        for wh in wheres:
            for ec in ("dsa_only", "grch38_diff_ref"):
                sel = [r for r in rows if r["where"] == wh and r["event_class"] == ec]
                if not sel:
                    continue
                cells = []
                for st in sets:
                    lab = Counter(r[st].split("|")[0] for r in sel)
                    n_t = len(sel) - lab["none"]
                    cells.append(f"{n_t:,} ({lab['somatic']} / {lab['non']} / {lab['germline']} / {lab['ignore']})")
                    summary[f"{g}|{wh}|{ec}|{st}"] = dict(events=len(sel), tensor=n_t, **{k: v for k, v in lab.items()})
                L.append(f"| {wh} | {ec} | {len(sel):,} | {sum(r['in_truth'] for r in sel):,} | " + " | ".join(cells) + " |")
        L += ["", "Reasons of the tensors' labels (all events of this graph):", ""]
        for st in sets:
            rs = Counter(r[st] for r in rows if not r[st].startswith("none"))
            L.append(f"- {st}: " + "; ".join(f"{k.replace('|', ' / ')} {v}" for k, v in rs.most_common()))
        L += ["", "Baseline, all COLO829T SMaHT truth SNVs on chr1-22 in the label region "
              f"({len(truth_in_bed):,}): the representative allele of a tensor labelled somatic",
              ""] + [f"- {st}: {len(truth_a1[st] & truth_in_bed):,} ({100 * len(truth_a1[st] & truth_in_bed) / len(truth_in_bed):.1f}%)"
                     for st in sets] + [""]
        for st in sets:
            summary[f"baseline|{st}"] = dict(truth_in_bed=len(truth_in_bed), a1_somatic=len(truth_a1[st] & truth_in_bed))
    (HERE / "dsa_only_labels.md").write_text("\n".join(L) + "\n")
    (HERE / "dsa_only_labels.json").write_text(json.dumps(summary, indent=1) + "\n")
    print("\n".join(L))


if __name__ == "__main__":
    {"windows": windows, "tensors": tensors, "labels": labels}[sys.argv[1]]()
