"""Tables and numbers for REPORT.md from reads_part*.jsonl + snv_detail.jsonl + snv_nodes.json (HG008-T Illumina WGS, v3 tensors).

Writes:
  indel_no_candidate.tsv                 every INS/DEL truth allele without an Illumina candidate (one row each)
  snv_no_candidate.tsv                   every SNV truth allele without an Illumina candidate (+ SNV detail columns)
  controls.tsv                           truth alleles that got an Illumina tensor, analysed the same way
  residual_edits_labelled_germline.tsv   residual edits (>= 3 ALT reads) of misses whose tensor is labelled germline
  transitions.tsv                        PacBio v6, ONT and Illumina status / class for every truth allele
  summary.json                           counts used in the report
"""
import csv, json, re, sys
from collections import Counter, defaultdict
from pathlib import Path

import numpy as np

csv.field_size_limit(sys.maxsize)
HERE = Path(__file__).resolve().parent
RUN = Path("/scratch/jshen/data/pansoma_v2_tensors/Liss_lab_BCM_Illumina-WGS_20240313")
PACBIO = Path("/scratch/jshen/data/pansoma_v2_tensors/Liss_lab_PacBio_Revio_20240125")
PACBIO_ANALYSIS = HERE.parent / "somatic_miss_analysis_v6"
ONT = Path("/scratch/jshen/data/pansoma_v2_tensors/Liss_lab_Northeastern-ONT-UL-20241216")
ONT_ANALYSIS = HERE.parent / "somatic_miss_analysis_ont"
TRUTH_DIR = RUN / "truth"
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
    """{candidate_id of every allele of every labelled Illumina tensor: (label name, reason, representative?)}."""
    result = {}
    ids = re.compile(r'"candidate_id": "([^"]+)"')
    for kind in ("SNV", "INDEL"):
        for f in sorted((RUN / "v3_tensors" / kind).glob("chr*_labels.ndjson")):
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
    """Edits the ALT reads make at the site (>= min_reads reads), with the label their Illumina tensor got."""
    out = []
    for c, d in r["alt_candidates"].items():
        if d["reads"] >= min_reads:
            lab = LABEL_OF.get(c)
            out.append((c, d["reads"], (lab[0] + ("" if lab[2] else "(other site allele)") + ":" + lab[1]) if lab else "no tensor"))
    return out


def bucket(n):
    return "1" if n == 1 else "2-5" if n <= 5 else "6-20" if n <= 20 else "21-50" if n <= 50 else ">50"


snv_detail = {json.loads(l)["truth_id"]: json.loads(l) for l in open(HERE / "snv_detail.jsonl")}
snv_nodes = json.load(open(HERE / "snv_nodes.json"))


def snv_class(r):
    """S1-S4c from the SNV detail (same-length reads, base at the SNV, node used) and the read-level reason."""
    d = snv_detail[r["truth_id"]]
    if d["records"] == 0 or d["low_mapq"] / d["records"] >= 0.8 or d["spanning"] == 0:
        return "S1_low_mapq_or_no_reads"
    if d["same_length"] / d["spanning"] >= 0.5:
        alt = d["base"].get("ALT", 0)
        if alt < 3:
            return "S2_alt_absent_in_reads"
        if d["alt_via"].get("other_node_or_skip", 0) > d["alt_via"].get("grch38_node", 0):
            return "S3_alt_is_existing_graph_allele"
    if r["reason"] in ("no_ALT_read", "fewer_than_3_ALT_reads"):
        return "S4c_repeat_no_clear_alt"
    return "S4a_repeat_snv_edit_on_branch" if any(":SNP:" in c for c, _, _ in r["residual"]) else "S4b_repeat_absorbed_by_graph_path"


for r in rows:
    t = truth[r["truth_id"]]
    pos0 = int(t["pos0"])
    r["reason"], r["from_exact_reads"] = primary(r)
    r["residual"] = residual(r)
    r["final_class"] = snv_class(r) if r["kind"] == "SNP" else CLASS[r["reason"]]
    if r["kind"] == "SNP":
        d, n = snv_detail[r["truth_id"]], snv_nodes.get(r["truth_id"], {})
        r["reads_same_length_as_ref"] = f"{d['same_length']}/{d['spanning']}"
        r["base_at_snv_in_same_length_reads"] = json.dumps(d["base"])
        r["alt_reads_node_at_snv"] = json.dumps(d["alt_via"])
        r.update(node_length=n.get("node_length"), node_is_target=n.get("target"),
                 discovery_perfect=n.get("perfect"), discovery_not_perfect=n.get("not_perfect"))
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


SNV_COLUMNS = ["reads_same_length_as_ref", "base_at_snv_in_same_length_reads", "alt_reads_node_at_snv", "node_length",
               "node_is_target", "discovery_perfect", "discovery_not_perfect"]


def write(path, items, extra=()):
    with open(HERE / path, "w", newline="") as out:
        w = csv.writer(out, delimiter="\t", lineterminator="\n")
        w.writerow(COLUMNS + list(extra) + ["residual_edits(reads)=tensor_label", "alt_read_representation", "context"])
        for r in items:
            w.writerow([r[c] for c in COLUMNS] + [r.get(c, "") for c in extra]
                       + ["; ".join(f"{c}({n})={lab}" for c, n, lab in r["residual"]),
                          json.dumps(r["alt_representation"]), r["context"]])


misses = [r for r in rows if r["status"] == "no_candidate"]
controls = [r for r in rows if r["status"] != "no_candidate"]
order = lambda r: (r["kind"], r["final_class"], r["chrom"], r["vcf_pos"])  # noqa: E731
write("indel_no_candidate.tsv", sorted((r for r in misses if r["kind"] != "SNP"), key=order))
write("snv_no_candidate.tsv", sorted((r for r in misses if r["kind"] == "SNP"), key=order), SNV_COLUMNS)
write("controls.tsv", sorted(controls, key=order))
with open(HERE / "residual_edits_labelled_germline.tsv", "w") as out:
    out.write("chrom\tvcf_pos\tvcf_ref\tvcf_alt\tkind\tfinal_class\tresidual_candidate\talt_reads\ttensor_label\n")
    for r in sorted(misses, key=order):
        for c, n, lab in r["residual"]:
            if lab.startswith("germline"):
                out.write(f"{r['chrom']}\t{r['vcf_pos']}\t{r['vcf_ref']}\t{r['vcf_alt']}\t{r['kind']}\t{r['final_class']}"
                          f"\t{c}\t{n}\t{lab}\n")

# ---- PacBio v6 / ONT -> Illumina transitions ---------------------------------------------
def statuses(path):
    return {r["truth_id"]: r["status"] for r in csv.DictReader(open(path), delimiter="\t")}


def classes(folder):
    found = {}
    for f in ("indel_no_candidate.tsv", "snv_no_candidate.tsv"):
        for r in csv.DictReader(open(folder / f), delimiter="\t"):
            found[(r["chrom"], r["vcf_pos"], r["vcf_ref"], r["vcf_alt"])] = r["final_class"]
    return found


pb_status, ont_status = statuses(PACBIO / "v6_tensors/somatic.recall.tsv"), statuses(ONT / "v3_tensors/somatic.recall.tsv")
il_status = statuses(RUN / "v3_tensors/somatic.recall.tsv")
pb_class, ont_class = classes(PACBIO_ANALYSIS), classes(ONT_ANALYSIS)
il_class = {r["truth_id"]: r["final_class"] for r in misses}
with open(HERE / "transitions.tsv", "w") as out:
    out.write("truth_id\tchrom\tvcf_pos\tkind\tpassed\tin_bed\tpacbio_status\tpacbio_class\tont_status\tont_class"
              "\tillumina_status\tillumina_class\n")
    for tid, t in truth.items():
        key = (t["chrom"], t["vcf_pos"], t["vcf_ref"], t["vcf_alt"])
        out.write(f"{tid}\t{t['chrom']}\t{t['vcf_pos']}\t{t['kind']}\t{t['passed']}\t{t['in_bed']}\t{pb_status[tid]}"
                  f"\t{pb_class.get(key, '')}\t{ont_status[tid]}\t{ont_class.get(key, '')}\t{il_status[tid]}"
                  f"\t{il_class.get(tid, '')}\n")

# ---- numbers for the report ----------------------------------------------------------------
summary = dict(loci=len(rows), misses=len(misses), controls=len(controls), by_kind={})
for label, group in (("miss", misses), ("control", controls)):
    for kind in ("SNP", "INS", "DEL"):
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
