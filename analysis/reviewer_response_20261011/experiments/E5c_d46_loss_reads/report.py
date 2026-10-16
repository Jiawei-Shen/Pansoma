"""E5c tables: how the ALT reads of the truth alleles lost on d46 (against GRCh38-only) are aligned on d46.

Reads DATA/<set>/reads*.jsonl (reads.py) for the six d46 sets, and joins
  * the repeat class of each truth INDEL (E7 s2 sample alleles: HP>=7 / STR / other),
  * the PoN tag and the d46 spelling test (E5b per_truth.tsv.gz: d46_spelled = a d46 path spells GRCh38 + the truth
    allele between the window anchors, independent of the reads).
Per locus, the ALT-like reads are counted by representation (reads.py):
  entire  = site_bypassed:{branch|skip_edge}:no_edit        (the read walks graph alt nodes / a skip edge, no edit)
  partial = site_bypassed:{branch|skip_edge}:edit_on_path_taken (alt nodes / skip edge plus a residual edit)
  grch38_other = site_on_grch38:other_edit                    (stays on GRCh38 nodes, a different edit)
  truth_edit = key_observed_*                                 (the truth edit itself is made on d46)
  other = site_on_grch38:no_local_edit, complex, longer than 50 bp
and the locus gets the class of the majority representation of its error-free ALT reads (else of all ALT-like reads),
as in analysis/hg008_somatic_miss_20260926 (I1-I10).
Outputs: tables.md, summary.json (this folder), DATA/per_locus.tsv.gz.
"""
import csv, gzip, json, sys
from collections import Counter, defaultdict
from pathlib import Path

csv.field_size_limit(sys.maxsize)
HERE = Path(__file__).resolve().parent
DATA = Path("/scratch/jshen/data/reviewer_experiments_20261013/E5c_d46_loss_reads")
E5B = "/scratch/jshen/data/reviewer_experiments_20261013/E5b_absorbed_d46/per_truth.tsv.gz"
E7 = "/scratch/jshen/data/reviewer_experiments_20261013/E7_graph_composition/sample_alleles/{}.somatic.tsv.gz"
SETS = ["HG008T_Illumina_d46", "HG008T_PacBio_d46", "HG008T_ONT_d46",
        "COLO829T_Illumina_d46", "COLO829T_fiberseq_d46", "COLO829T_ONT_d46"]
CLASS = {
    "no_read_spans_site": "I5_low_mapq_or_no_reads",
    "no_ALT_read": "I4_no_or_few_alt_reads",
    "fewer_than_3_ALT_reads": "I4_no_or_few_alt_reads",
    "site_bypassed:branch:no_edit": "I1_entire_on_graph_alt_path",
    "site_bypassed:skip_edge:no_edit": "I1_entire_on_graph_alt_path",
    "site_bypassed:branch:edit_on_path_taken": "I2_partial_alt_path_plus_edit",
    "site_bypassed:skip_edge:edit_on_path_taken": "I2_partial_alt_path_plus_edit",
    "site_on_grch38:other_edit": "I3_grch38_path_different_edit",
    "indel_longer_than_50": "I6_longer_than_50bp",
    "site_on_grch38:no_local_edit": "I7_unclear",
    "complex_replacement_edit": "I8_complex_replacement_edit",
    "key_observed_on_non_target_node": "I9_truth_edit_on_unbuilt_node",
    "key_observed_on_target_node": "I10_truth_edit_made",
    "key_observed_low_base_quality": "I10_truth_edit_made",
}
CLASSES = ["I1_entire_on_graph_alt_path", "I2_partial_alt_path_plus_edit", "I3_grch38_path_different_edit",
           "I10_truth_edit_made", "I4_no_or_few_alt_reads", "I5_low_mapq_or_no_reads", "I6_longer_than_50bp",
           "I7_unclear", "I8_complex_replacement_edit", "I9_truth_edit_on_unbuilt_node"]
SHORT = {c: c.split("_")[0] for c in CLASSES}


def group(kind):
    if kind.startswith("site_bypassed") and kind.endswith(":no_edit"):
        return "entire"
    if kind.startswith("site_bypassed"):
        return "partial"
    if kind == "site_on_grch38:other_edit":
        return "grch38_other"
    if kind.startswith("key_observed"):
        return "truth_edit"
    return "other"


def primary(r):
    alt = r["alt_exact"] + r["alt_like"]
    if r["spanning"] == 0:
        return "no_read_spans_site"
    if alt == 0:
        return "no_ALT_read"
    if alt < 3:
        return "fewer_than_3_ALT_reads"
    exact = {k: v for k, v in r["alt_representation"].items() if k.endswith(":exact_read")}
    rep = Counter()
    for k, v in (exact or r["alt_representation"]).items():
        rep[k.rsplit(":", 1)[0]] += v
    return rep.most_common(1)[0][0]


def load_joins():
    e5b = {}
    with gzip.open(E5B, "rt") as fh:
        hdr = fh.readline().rstrip("\n").split("\t")
        for line in fh:
            f = dict(zip(hdr, line.rstrip("\n").split("\t")))
            e5b[(f["sample"], f["truth_id"])] = (f["pon"] != "", f["d46_spelled"] == "1")
    rep = {}
    for s5, s7 in (("HG008T", "HG008"), ("COLO829T", "COLO829")):
        with gzip.open(E7.format(s7), "rt") as fh:
            fh.readline()
            for line in fh:
                f = line.rstrip("\n").split("\t")
                rep[(s5, f[0])] = {"HP2-6": "other"}.get(f[7], f[7]) if f[6] in ("INS", "DEL") else ""
    return e5b, rep


def pct(a, b):
    return f"{100 * a / b:.1f}" if b else "-"


def main():
    e5b, rep = load_joins()
    per_locus, summary, md = [], {}, []
    for s in SETS:
        tumor = s.split("_")[0]
        files = sorted((DATA / s).glob("reads_part*.jsonl")) + sorted((DATA / s).glob("reads.jsonl"))
        rows = list({r["truth_id"]: r for f in files for r in (json.loads(l) for l in open(f) if l.strip())}.values())
        if not rows:
            continue
        for r in rows:
            g = Counter()
            for k, v in r["alt_representation"].items():
                g[group(k.rsplit(":", 1)[0])] += v
            pon, spelled = e5b.get((tumor, r["truth_id"]), (False, False))
            per_locus.append(dict(set=s, truth_id=r["truth_id"], locus_set=r["set"], kind=r["kind"],
                                  status_linear=r["status_linear"], status_d46=r["status"],
                                  repeat=rep.get((tumor, r["truth_id"]), ""), pon_tagged=int(pon), d46_spelled=int(spelled),
                                  spanning=r["spanning"], alt_reads=r["alt_exact"] + r["alt_like"],
                                  entire=g["entire"], partial=g["partial"], grch38_other=g["grch38_other"],
                                  truth_edit=g["truth_edit"], other=g["other"],
                                  primary=primary(r), cls=CLASS.get(primary(r), primary(r))))
    with gzip.open(DATA / "per_locus.tsv.gz", "wt") as fo:
        w = csv.DictWriter(fo, fieldnames=list(per_locus[0]), delimiter="\t")
        w.writeheader()
        w.writerows(per_locus)

    md.append("# E5c tables: truth alleles lost on HPRC v2.1 d46 against GRCh38-only, by how their ALT reads align on d46\n")
    md.append("Lost = a tensor allele equals the truth on GRCh38-only (status tensor_representative or "
              "tensor_non_representative_allele) and none does on d46 (no_candidate or filtered). chr1-22, all truth "
              "records. Class = majority representation of the locus's error-free ALT reads (reads.py / report.py).\n")
    for kind_label, kinds in (("INDEL", ("INS", "DEL")), ("SNV", ("SNP",))):
        md.append(f"\n## {kind_label}: class of each lost locus (% of lost loci)\n")
        md.append("| set | d46 status | n lost | " + " | ".join(SHORT[c] for c in CLASSES) + " | I1 + I2 |")
        md.append("|---|---|---|" + "---|" * len(CLASSES) + "---|")
        for s in SETS:
            for st in ("all", "no_candidate", "filtered"):
                L = [x for x in per_locus if x["set"] == s and x["locus_set"] == "lost" and x["kind"] in kinds
                     and (st == "all" or x["status_d46"] == st)]
                if not L:
                    continue
                c = Counter(x["cls"] for x in L)
                md.append(f"| {s} | {st} | {len(L):,} | " + " | ".join(pct(c[k], len(L)) for k in CLASSES)
                          + f" | {pct(c[CLASSES[0]] + c[CLASSES[1]], len(L))} |")
                summary.setdefault(kind_label, {}).setdefault(s, {})[st] = dict(n=len(L), **{SHORT[k]: c[k] for k in CLASSES})
        md.append(f"\n## {kind_label}: lost loci whose ALT reads go through graph alt nodes (entire or partial)\n")
        md.append("| set | n lost with >= 3 ALT reads | any ALT read on alt nodes | majority of ALT reads on alt nodes | "
                  "ALT reads pooled: entire / partial / GRCh38 other edit / truth edit / other (%) |")
        md.append("|---|---|---|---|---|")
        for s in SETS:
            L = [x for x in per_locus if x["set"] == s and x["locus_set"] == "lost" and x["kind"] in kinds and x["alt_reads"] >= 3]
            if not L:
                continue
            anyalt = sum(1 for x in L if x["entire"] + x["partial"] > 0)
            major = sum(1 for x in L if 2 * (x["entire"] + x["partial"]) > x["alt_reads"])
            tot = Counter()
            for x in L:
                for k in ("entire", "partial", "grch38_other", "truth_edit", "other"):
                    tot[k] += x[k]
            T = sum(tot.values())
            md.append(f"| {s} | {len(L):,} | {pct(anyalt, len(L))} | {pct(major, len(L))} | "
                      + " / ".join(pct(tot[k], T) for k in ("entire", "partial", "grch38_other", "truth_edit", "other")) + " |")
    md.append("\n## INDEL: class by repeat context (lost loci; % of the class row)\n")
    md.append("| set | repeat | n lost | I1 | I2 | I3 | I10 | I4 + I5 | rest |")
    md.append("|---|---|---|---|---|---|---|---|---|")
    for s in SETS:
        for rc in ("HP>=7", "STR", "other"):
            L = [x for x in per_locus if x["set"] == s and x["locus_set"] == "lost" and x["kind"] in ("INS", "DEL") and x["repeat"] == rc]
            if not L:
                continue
            c = Counter(SHORT[x["cls"]] if x["cls"] in SHORT else x["cls"] for x in L)
            rest = len(L) - c["I1"] - c["I2"] - c["I3"] - c["I10"] - c["I4"] - c["I5"]
            md.append(f"| {s} | {rc} | {len(L):,} | {pct(c['I1'], len(L))} | {pct(c['I2'], len(L))} | {pct(c['I3'], len(L))} | "
                      f"{pct(c['I10'], len(L))} | {pct(c['I4'] + c['I5'], len(L))} | {pct(rest, len(L))} |")
    md.append("\n## INDEL: lost loci not PoN-tagged (what the graph costs after population filtering), by class\n")
    md.append("| set | n lost, not tagged | I1 | I2 | I3 | I10 | I4 + I5 | rest |")
    md.append("|---|---|---|---|---|---|---|---|")
    for s in SETS:
        L = [x for x in per_locus if x["set"] == s and x["locus_set"] == "lost" and x["kind"] in ("INS", "DEL") and not x["pon_tagged"]]
        if not L:
            continue
        c = Counter(SHORT[x["cls"]] if x["cls"] in SHORT else x["cls"] for x in L)
        rest = len(L) - c["I1"] - c["I2"] - c["I3"] - c["I10"] - c["I4"] - c["I5"]
        md.append(f"| {s} | {len(L):,} | {pct(c['I1'], len(L))} | {pct(c['I2'], len(L))} | {pct(c['I3'], len(L))} | "
                  f"{pct(c['I10'], len(L))} | {pct(c['I4'] + c['I5'], len(L))} | {pct(rest, len(L))} |")
    md.append("\n## Controls (tensor on both sets): class of the locus (%), same method\n")
    md.append("| set | kind | n | I10 (truth edit made) | I1 | I2 | I3 | other |")
    md.append("|---|---|---|---|---|---|---|---|")
    for s in SETS:
        for kl, kinds in (("INDEL", ("INS", "DEL")), ("SNV", ("SNP",))):
            L = [x for x in per_locus if x["set"] == s and x["locus_set"] == "control" and x["kind"] in kinds]
            if not L:
                continue
            c = Counter(SHORT[x["cls"]] if x["cls"] in SHORT else x["cls"] for x in L)
            other = len(L) - c["I10"] - c["I1"] - c["I2"] - c["I3"]
            md.append(f"| {s} | {kl} | {len(L)} | {pct(c['I10'], len(L))} | {pct(c['I1'], len(L))} | {pct(c['I2'], len(L))} | "
                      f"{pct(c['I3'], len(L))} | {pct(other, len(L))} |")
    (HERE / "tables.md").write_text("\n".join(md) + "\n")
    (HERE / "summary.json").write_text(json.dumps(summary, indent=1))
    print("wrote", HERE / "tables.md", len(per_locus), "loci")


if __name__ == "__main__":
    main()
