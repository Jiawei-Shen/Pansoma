#!/usr/bin/env python3
"""ver4.6 (2026-11-26) of Response Tables 3.8 (= 2.6) and 3.9: the false-positive decomposition of
../../ver4_20261123/fpdecomp/fp_decomp.py with the HPRC v2.1 side on the six-tumor models WITHOUT the population site
features as a model input (seed 0), so that neither reference's model takes them:

  GRCh38-only  R/m4m6_graphs_20261110/multicancer6_<P>_<KIND>_linear_nositef/vcf_[ppr8_|pp_]chr1_<TUMOR>   (unchanged)
  HPRC v2.1    SNV   R/m6_nositef_graphs_20261118/multicancer6_<P>_SNV_d46_nositef/vcf_chr1_<SET>
                     (trained and called on RIS, pulled to the lab cluster 2026-11-26)
               INDEL R/m6_nositef_graphs_20261118/multicancer6_<P>_INDEL_d46_nositef/vcf_[ppr8_|pp_]chr1_<TUMOR>
                     (lab cluster; ONT INDEL: pp)

Everything else is fp_decomp.py's code (imported, not copied): compute(), cells_json(), the tables and the ranges.
Three functions are replaced:
- ver4_cells: the HPRC v2.1 call-set directories above;
- set_of: the RIS SNV call sets name their files <SET>.tensors.SNV.graph.vcf.gz (the INDEL ones <SET>.INDEL...); the
  ".tensors" suffix is dropped so that the v2.1 test set is found;
- check_reports: the RIS reports list the four population VCFs under /storage2; the two arms' pon_files are compared
  by file name (same four names in the same order) instead of by path. The other checks are unchanged.

Which columns depend on the HPRC v2.1 MODEL: only the last one (HPRC v2.1 FP = records of the v2.1 arm's
raw_calls/fp.vcf.gz). "Candidates retained on HPRC v2.1" is `k in c["tk"]`, and c["tk"] = tensor_keys(c["d46_set"]),
the chr1 tensors of the v2.1 TEST SET (chr1_labels.ndjson), not the model's calls; columns 1-5 otherwise use the
GRCh38-only arm, the germline VCF and the v2.1 haplotype carriers.

usage: /wanglab/jshen/anaconda3/bin/python fp_decomp_nositef.py   (writes next to itself; then check_nositef.py)
"""
import json
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
OLD = HERE.parents[1] / "ver4_20261123" / "fpdecomp"
sys.path.insert(0, str(OLD))
import fp_decomp as F  # noqa: E402

R = F.R
NS = R / "m6_nositef_graphs_20261118"


def cells():
    out = []
    for kind in F.KINDS:
        for tumor in F.TUMORS:
            for p in F.PLATFORMS:
                r = "" if kind == "SNV" else ("pp_" if p == "ONT" else "ppr8_")
                lin = R / "m4m6_graphs_20261110" / f"multicancer6_{p}_{kind}_linear_nositef" / f"vcf_{r}chr1_{tumor}"
                run = NS / f"multicancer6_{p}_{kind}_d46_nositef"
                if kind == "SNV":
                    d = sorted(run.glob(f"vcf_chr1_{tumor}_{F.ds(tumor, p)}_d46*"))
                    assert len(d) == 1, (run, tumor, d)
                    d46 = d[0]
                else:
                    d46 = run / f"vcf_{r}chr1_{tumor}"
                out.append(dict(tumor=tumor, platform=p, kind=kind, lin_dir=lin, d46_dir=d46))
    return out


_set_of = F.set_of


def set_of(vcf_dir, kind):
    s = _set_of(vcf_dir, kind)
    return s[: -len(".tensors")] if s.endswith(".tensors") else s


def check_reports(c):
    """fp_decomp.check_reports with pon_files compared by file name."""
    out = {}
    reps = {arm: json.loads((c[f"{arm}_dir"] / "eval_bed" / "report.json").read_text()) for arm in ("lin", "d46")}
    for key in ("truth", "bed", "truth_records"):
        assert reps["lin"][key] == reps["d46"][key], (c["lin_dir"], key, reps["lin"][key], reps["d46"][key])
    names = {arm: [Path(x).name for x in reps[arm]["pon_files"]] for arm in reps}
    assert names["lin"] == names["d46"] and len(names["lin"]) == 4, (c["lin_dir"], names)
    for stage, which in F.STAGES.items():
        truth = set()
        for arm in ("lin", "d46"):
            a = reps[arm][stage]["at_threshold"]
            base = reps[arm][stage]["baseline"]
            assert a["tp_baseline"] + a["fn"] == base, (c[f"{arm}_dir"], stage, a, base)
            d = c[f"{arm}_dir"] / "eval_bed" / "vcfeval" / which
            assert (d / "done").exists(), d
            n_tp, n_fn = F.n_records(d / "tp-baseline.vcf.gz"), F.n_records(d / "fn.vcf.gz")
            assert (n_tp, n_fn) == (a["tp_baseline"], a["fn"]), (d, n_tp, n_fn, a)
            truth.add(base)
            out[f"{arm}_{stage}"] = dict(tp=a["tp_baseline"], fn=a["fn"], fp=a["fp"], truth=base,
                                         precision=a["precision"], recall=a["recall"], f1=a["f1"])
        assert len(truth) == 1, (c["lin_dir"], stage, truth)
    for arm in ("lin", "d46"):
        out[f"{arm}_call_rule"] = reps[arm].get("steps", {}).get("graph_vcf", {}).get("call")
        out[f"{arm}_threshold"] = reps[arm].get("threshold")
        out[f"{arm}_checkpoint"] = reps[arm].get("steps", {}).get("graph_vcf", {}).get("checkpoint")
    out["truth_vcf"], out["truth_bed"] = reps["lin"]["truth"], reps["lin"]["bed"]
    return out


F.set_of = set_of
F.check_reports = check_reports


def main():
    cs = cells()
    rows = F.compute(cs)
    with open(HERE / "fp_decomp.tsv", "w") as f:
        f.write("\t".join(rows[0]) + "\n")
        for r in rows:
            f.write("\t".join(str(v) for v in r.values()) + "\n")
    J = F.cells_json(cs, rows)
    (HERE / "fp_decomp.json").write_text(json.dumps(J, indent=1) + "\n")
    md, r = F.write_tables(J, HERE)
    # the headers written by fp_decomp.write_tables name the ver4 models: rewrite the two header lines
    t = (HERE / "tables.md").read_text().split("\n")
    assert t[0].startswith("# False-positive decomposition") and t[2].startswith("chr1, benchmark regions")
    t[0] = "# False-positive decomposition, six-tumor models without population site features on both references (ver4.6, 2026-11-26)"
    t[2] = t[2].replace(
        "HPRC v2.1 = six-tumor current model (SNV snv_newpipe_20261109/multicancer6_<P>_SNV_newpipe, INDEL release_s0_<P>_INDEL).",
        "HPRC v2.1 = six-tumor model without population site features as a model input "
        "(m6_nositef_graphs_20261118/multicancer6_<P>_<KIND>_d46_nositef; INDEL PASS ppr8, ONT pp).")
    assert "d46_nositef" in t[2]
    (HERE / "tables.md").write_text("\n".join(t))
    rr = (HERE / "ranges.md").read_text().replace("recomputed on the ver4 tables", "recomputed on the ver4.6 tables (no site features on either reference)")
    (HERE / "ranges.md").write_text(rr)
    print("\n".join(t))
    print(rr)


if __name__ == "__main__":
    main()
