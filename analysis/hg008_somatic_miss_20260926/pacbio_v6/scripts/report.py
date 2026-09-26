"""Tables and numbers for REPORT.md from reads_part*.jsonl (v6 INDEL misses + controls).

Writes:
  indel_no_candidate.tsv                 every INS/DEL truth allele without a v6 candidate (one row each)
  controls.tsv                           INDEL truth alleles that got a v6 tensor, analysed the same way
  residual_edits_labelled_germline.tsv   residual edits (>= 3 ALT reads) of misses whose tensor is labelled germline
  transitions.tsv                        v5 status / v5 class -> v6 status / v6 class for every INDEL truth allele
  summary.json                           counts used in the report
"""
import csv, json, re, sys
from collections import Counter, defaultdict
from pathlib import Path

import numpy as np

csv.field_size_limit(sys.maxsize)
HERE = Path(__file__).resolve().parent
RUN = Path("/scratch/jshen/data/pansoma_v2_tensors/Liss_lab_PacBio_Revio_20240125")
TRUTH_DIR = Path("/scratch/jshen/data/pansoma_v2_tensors/truth")
rows = list({r["truth_id"]: r for f in sorted(HERE.glob("reads_part*.jsonl"))
             for r in (json.loads(l) for l in open(f) if l.strip())}.values())
truth = {r["truth_id"]: r for r in csv.DictReader(open(TRUTH_DIR / "somatic.graph.tsv"), delimiter="\t")}
somatic_pos = defaultdict(list)
for t in truth.values():
    somatic_pos[t["chrom"]].append(int(t["pos0"]))
somatic_pos = {c: np.array(sorted(v)) for c, v in somatic_pos.items()}
germ, germ_alleles = defaultdict(list), {}
with open(TRUTH_DIR / "germline.graph.tsv") as f:
    next(f)
    for line in f:
        x = line.split("\t", 14)
        germ[x[1]].append(int(x[5]))
        germ_alleles.setdefault((x[1], int(x[5])), []).append(
            f"{x[6] or '-'}>{x[7] or '-'}({x[9]},{'PASS' if x[13] == 'True' else x[12]})")
germ = {c: np.array(sorted(v)) for c, v in germ.items()}


def tensor_labels():
    """{candidate_id of every allele of every labelled v6 tensor: (label name, reason, representative?)}."""
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


LABEL_OF = tensor_labels()

# Primary reason -> class (the v5 report's mapping; the last three reasons did not occur among v5 misses).
CLASS = {
    "no_read_spans_site": "I5_low_mapq_or_no_reads",
    "no_ALT_read": "I4_no_or_few_alt_reads",
    "fewer_than_3_ALT_reads": "I4_no_or_few_alt_reads",
    "site_bypassed:branch:no_edit": "I1_allele_fully_in_graph",
    "site_bypassed:skip_edge:no_edit": "I1_allele_fully_in_graph",
    "site_bypassed:branch:edit_on_path_taken": "I2_partly_in_graph_residual_edit",
    "site_bypassed:skip_edge:edit_on_path_taken": "I2_partly_in_graph_residual_edit",
    "site_on_grch38:other_edit": "I3_grch38_path_different_edit",
    "indel_longer_than_50": "I6_longer_than_50bp",
    "site_on_grch38:no_local_edit": "I7_unclear",
    "complex_replacement_edit": "I8_complex_replacement_edit",
    "key_observed_on_non_target_node": "I9_truth_edit_on_unbuilt_node",
    "key_observed_on_target_node": "I10_truth_edit_on_built_node",
    "key_observed_low_base_quality": "I10_truth_edit_on_built_node",
}


def count_near(table, chrom, pos0, d=50, exclude_self=False):
    a = table.get(chrom, np.array([]))
    n = int(np.searchsorted(a, pos0 + d, "right") - np.searchsorted(a, pos0 - d, "left"))
    return n - 1 if exclude_self else n


def primary(r):
    """Majority representation among the error-free ALT reads (else among all ALT-like reads)."""
    alt = r["alt_exact"] + r["alt_like"]
    if r["spanning"] == 0:
        return "no_read_spans_site", False
    if alt == 0:
        return "no_ALT_read", False
    if alt < 3:
        return "fewer_than_3_ALT_reads", False
    exact = {k: v for k, v in r["alt_representation"].items() if k.endswith(":exact_read")}
    rep = Counter()
    for k, v in (exact or r["alt_representation"]).items():
        rep[k.rsplit(":", 1)[0]] += v
    return rep.most_common(1)[0][0], bool(exact)


def residual(r, min_reads=3):
    """Edits the ALT reads make at the site (>= min_reads reads), with the label their v6 tensor got."""
    out = []
    for c, d in r["alt_candidates"].items():
        if d["reads"] >= min_reads:
            lab = LABEL_OF.get(c)
            out.append((c, d["reads"], (lab[0] + ("" if lab[2] else "(other site allele)") + ":" + lab[1]) if lab else "no tensor"))
    return out


def bucket(n):
    return "1" if n == 1 else "2-5" if n <= 5 else "6-20" if n <= 20 else "21-50" if n <= 50 else ">50"


for r in rows:
    t = truth[r["truth_id"]]
    pos0 = int(t["pos0"])
    r["reason"], r["from_exact_reads"] = primary(r)
    r["final_class"] = CLASS[r["reason"]]
    r["residual"] = residual(r)
    r["somatic_within_50bp"] = count_near(somatic_pos, r["chrom"], pos0, exclude_self=True)
    r["germline_within_50bp"] = count_near(germ, r["chrom"], pos0)
    r["germline_same_position"] = ";".join(germ_alleles.get((r["chrom"], pos0), []))
    r["germline_same_allele"] = any(a.startswith(f"{t['ref'] or '-'}>{t['alt'] or '-'}(")
                                    for a in germ_alleles.get((r["chrom"], pos0), []))
    r["passed"], r["in_bed"] = t["passed"], t["in_bed"]
    alt = r["alt_exact"] + r["alt_like"]
    r["alt_fraction"] = round(alt / r["spanning"], 3) if r["spanning"] else ""

COLUMNS = ["final_class", "chrom", "vcf_pos", "vcf_ref", "vcf_alt", "kind", "length", "gt", "passed", "in_bed", "reason",
           "from_exact_reads", "records", "low_mapq", "spanning", "ref_exact", "ref_like", "alt_exact", "alt_like", "unclear",
           "alt_fraction", "placements", "somatic_within_50bp", "germline_within_50bp", "germline_same_position",
           "germline_same_allele"]


def write(path, items):
    with open(HERE / path, "w", newline="") as out:
        w = csv.writer(out, delimiter="\t", lineterminator="\n")
        w.writerow(COLUMNS + ["residual_edits(reads)=tensor_label", "alt_read_representation", "context"])
        for r in items:
            w.writerow([r[c] for c in COLUMNS] + ["; ".join(f"{c}({n})={lab}" for c, n, lab in r["residual"]),
                                                  json.dumps(r["alt_representation"]), r["context"]])


misses = [r for r in rows if r["status"] == "no_candidate"]
controls = [r for r in rows if r["status"] != "no_candidate"]
order = lambda r: (r["kind"], r["final_class"], r["chrom"], r["vcf_pos"])  # noqa: E731
write("indel_no_candidate.tsv", sorted(misses, key=order))
write("controls.tsv", sorted(controls, key=order))
with open(HERE / "residual_edits_labelled_germline.tsv", "w") as out:
    out.write("chrom\tvcf_pos\tvcf_ref\tvcf_alt\tkind\tfinal_class\tresidual_candidate\talt_reads\ttensor_label\n")
    for r in sorted(misses, key=order):
        for c, n, lab in r["residual"]:
            if lab.startswith("germline"):
                out.write(f"{r['chrom']}\t{r['vcf_pos']}\t{r['vcf_ref']}\t{r['vcf_alt']}\t{r['kind']}\t{r['final_class']}"
                          f"\t{c}\t{n}\t{lab}\n")

# ---- v5 -> v6 transitions -----------------------------------------------------------------
v5_status = {r["truth_id"]: r["status"] for r in csv.DictReader(open(HERE / "v5_reference/v5_somatic.recall.tsv"), delimiter="\t")}
v6_status = {r["truth_id"]: r["status"] for r in csv.DictReader(open(RUN / "v6_tensors/somatic.recall.tsv"), delimiter="\t")}
v5_class = {(r["chrom"], r["vcf_pos"], r["vcf_ref"], r["vcf_alt"]): r["final_class"]
            for r in csv.DictReader(open(HERE / "v5_reference/v5_indel_no_candidate.tsv"), delimiter="\t")}
v6_class = {r["truth_id"]: r["final_class"] for r in misses}
with open(HERE / "transitions.tsv", "w") as out:
    out.write("truth_id\tchrom\tvcf_pos\tkind\tpassed\tin_bed\tv5_status\tv5_class\tv6_status\tv6_class\n")
    for tid, t in truth.items():
        if t["kind"] not in ("INS", "DEL"):
            continue
        k5 = v5_class.get((t["chrom"], t["vcf_pos"], t["vcf_ref"], t["vcf_alt"]), "")
        out.write(f"{tid}\t{t['chrom']}\t{t['vcf_pos']}\t{t['kind']}\t{t['passed']}\t{t['in_bed']}\t{v5_status[tid]}\t{k5}"
                  f"\t{v6_status[tid]}\t{v6_class.get(tid, '')}\n")

# ---- numbers for the report ----------------------------------------------------------------
summary = dict(loci=len(rows), misses=len(misses), controls=len(controls), by_kind={})
for label, group in (("miss", misses), ("control", controls)):
    for kind in ("INS", "DEL"):
        sel = [r for r in group if r["kind"] == kind]
        detail = {}
        for key in ("final_class", "reason"):
            for why, n in Counter(r[key] for r in sel).items():
                sub = [r for r in sel if r[key] == why]
                fr = sorted(r["alt_fraction"] for r in sub if r["alt_fraction"] != "")
                detail[f"{key}={why}"] = dict(
                    n=n, pct=round(100 * n / len(sel), 1), median_alt_fraction=fr[len(fr) // 2] if fr else None,
                    from_exact_reads=sum(r["from_exact_reads"] for r in sub),
                    pass_in_bed=sum(r["passed"] == "True" and r["in_bed"] == "True" for r in sub),
                    by_length=dict(Counter(bucket(r["length"]) for r in sub)))
        summary["by_kind"][f"{label}:{kind}"] = dict(total=len(sel), groups=dict(sorted(detail.items(), key=lambda x: -x[1]["n"])))
res = defaultdict(Counter)
for r in misses:
    labels = {lab.split(":")[0] for _, _, lab in r["residual"]} or {"(no residual edit in >=3 ALT reads)"}
    for lab in labels:
        res[f"{r['kind']}:{r['final_class'][:2]}"][lab] += 1
summary["residual_edit_tensor_labels"] = {k: dict(v) for k, v in sorted(res.items())}
json.dump(summary, open(HERE / "summary.json", "w"), indent=1)
print(json.dumps({k: {g: d["n"] for g, d in v["groups"].items() if g.startswith("final_class")}
                  for k, v in summary["by_kind"].items()}, indent=1))
print(json.dumps(summary["residual_edit_tensor_labels"], indent=1))
