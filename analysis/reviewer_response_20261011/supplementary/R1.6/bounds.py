"""chr1 recall bounds for Supplementary Fig. R1.6 (HG008T, COLO829T, HCC1395T x three platforms, HPRC v2.1 d46 sets).

Per set and kind, over the PASS chr1 truth alleles whose REF span overlaps the truth BED (the BED of the sets' labels
manifest, which is also rtg's -e BED; overlap is rtg's rule: HCC1395T has 9 chr1 truth INDELs whose VCF position lies just
outside the BED but whose deleted bases are inside, so somatic.recall.tsv's in_bed (position only) gives 124, rtg 133):
  n         truth alleles
  tagged    truth alleles our population filter tags (E5b per_truth.tsv.gz 'pon')
  exact     truth alleles with an exact candidate: a tensor labelled 1 by representative_allele_is_somatic_truth
            carries the truth id (the E21 index caches)
  exact_untagged
-> bound_any = (n - tagged) / n (any caller after the filter), bound_exact = exact_untagged / n (Pansoma).
Same inputs and rules as E27 denominators/label_status.py, whose load() is imported; HG008T / COLO829T must reproduce
label_status.json (asserted). Checks against the rtg reports (Pansoma six-tumor INDEL seed 0, chr1):
  n == tp_baseline + fn of eval_bed (the evaluated truth; HG008T INDEL 588 against rtg's 587, the one-record difference
  already in Supplementary Table R1.6), and the PoN tags on the whole chr1 == report truth_pon_tags.
Writes bounds.json next to this script.
Run: /wanglab/jshen/anaconda3/envs/sambcfenv/bin/python bounds.py   (login node, a few minutes: loads 18 index caches)
"""
import csv
import gzip
import json
import sys
from collections import defaultdict
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
DEN = HERE.parents[1] / "experiments" / "E27_assessment_20261030" / "denominators"
sys.path.insert(0, str(DEN))
import label_status as LS  # noqa: E402

RUNS = Path("/scratch/jshen/data/pansoma_net_v2_runs")
SETS = [(t, p) for t in ("HG008T", "COLO829T", "HCC1395T") for p in ("Illumina", "fiberseq" if t == "COLO829T" else "PacBio", "ONT")]


def main():
    pon, per_truth = {}, defaultdict(list)
    for r in csv.DictReader(gzip.open(LS.PT, "rt"), delimiter="\t"):
        pon[(r["sample"], int(r["truth_id"]))] = bool(r["pon"])
        if r["chrom"] == "chr1" and r["passed"] == "True":
            per_truth[r["sample"]].append(r)
    old = json.loads((DEN / "label_status.json").read_text())
    out = {}
    for sample, plat in SETS:
        set_name = f"{sample}_{plat}_d46"
        for kind in ("SNV", "INDEL"):
            idx = LS.load(set_name, kind)
            a, rname = idx.arrays, np.array(idx.meta["reasons"])[idx.arrays["reason"]]
            exact_ids = set()
            for p in np.flatnonzero((a["label"] == 1) & (rname == LS.EXACT)).tolist():
                exact_ids.update(a["truth_ids"][a["truth_ptr"][p]:a["truth_ptr"][p + 1]].tolist())
            starts, ends = LS.somatic_bed(idx.dir)[1]["chr1"]
            chr1 = [r for r in per_truth[sample] if r["kind"] in LS.KINDS[kind]]

            def overlaps(r):
                s0 = int(r["vcf_pos"]) - 1
                k = np.searchsorted(starts, s0 + len(r["vcf_ref"]), side="left") - 1
                return k >= 0 and ends[k] > s0
            truth = [r for r in chr1 if overlaps(r)]
            tags = [pon.get((sample, int(r["truth_id"]))) for r in truth]
            assert None not in tags, (set_name, kind, "truth without a PoN tag")
            n, tagged = len(truth), sum(tags)
            exact_untagged = sum(int(r["truth_id"]) in exact_ids and not t for r, t in zip(truth, tags))
            res = dict(n=n, tagged=tagged, exact=sum(int(r["truth_id"]) in exact_ids for r in truth),
                       exact_untagged=exact_untagged, bound_any=(n - tagged) / n, bound_exact=exact_untagged / n,
                       tagged_whole_chr1=sum(bool(pon.get((sample, int(r["truth_id"])))) for r in chr1), n_whole_chr1=len(chr1))
            key = f"{set_name}|{kind}"
            if key in old:
                c = old[key]["chr1"]
                assert (n, tagged, exact_untagged) == (c["all"]["n"], c["all"]["_tagged"], c["untagged"]["exact"]), (key, res, c)
            ds = f"{sample}_{plat}"
            rep = json.loads((RUNS / f"release_s0_{'PacBio' if plat == 'fiberseq' else plat}_INDEL" / f"vcf_ppr8_chr1_{ds}_d46"
                              / "eval_bed" / "report.json").read_text()) if kind == "INDEL" else None
            if kind == "SNV":
                rep = json.loads((RUNS / f"multicancer6_{'PacBio' if plat == 'fiberseq' else plat}_SNV" /
                                  f"vcf_chr1_{ds}_d46_normfix" / "eval_bed" / "report.json").read_text())
            at = rep["pon"]["at_threshold"]
            res["rtg_truth"] = at["tp_baseline"] + at["fn"]
            res["report_truth_pon_tags"] = rep["truth_pon_tags"]["tagged"]
            res["report_truth_records"] = rep["truth_records"]
            out[key] = res
            known = (sample, kind, n - res["rtg_truth"]) == ("HG008T", "INDEL", 1)  # the table's one-record difference
            ok = (res["rtg_truth"] == n or known) and res["report_truth_pon_tags"] == res["tagged_whole_chr1"]
            flag = ("  (known one-record difference)" if known else "") if ok else "  <-- CHECK"
            print(f"{key}: n {n} (rtg {res['rtg_truth']}), tagged {tagged} (whole chr1 {res['tagged_whole_chr1']} of "
                  f"{len(chr1)}; report {res['report_truth_pon_tags']} of {res['report_truth_records']}), exact {res['exact']}, "
                  f"exact untagged {exact_untagged} -> {res['bound_exact']:.3f} / {res['bound_any']:.3f}{flag}", flush=True)
    (HERE / "bounds.json").write_text(json.dumps(out, indent=1) + "\n")


if __name__ == "__main__":
    main()
