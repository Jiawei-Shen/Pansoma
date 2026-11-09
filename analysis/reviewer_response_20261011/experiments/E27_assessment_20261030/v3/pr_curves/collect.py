"""R3.9 point 4: chr1 precision-recall curves of every Pansoma model family and both comparators, with operating points.

Panels: tumor (HG008T, COLO829T, HCC1395T) x platform (Illumina, PacBio HiFi / Fiber-seq, ONT) x kind (SNV, INDEL).
Curves: rtg vcfeval --all-records over every scored record (benchmark BED, chr1), before ('raw') and after ('pon') our
population filter, as written by pansoma_net_v2.vcfeval (report.json raw/pon curve_roc; each point carries p_somatic).

Model families (validation-selected checkpoints, as in the letter):
- single-tumor, three references (HG008T and COLO829T only; trained on the tumor's chr2-22):
  SNV  = Table R3.1 / Supplementary Table R3.9a models (`--norm train`, 3 seeds; ablation/main_runs.json 'train', the
         run directory recorded there: RIS copy or CPU fp32 re-prediction), PASS = validation recall 0.9 before the filter.
  INDEL = final recipe per tumor (RIS ppr8, analysis/indel_ppr8_9919d4a_20261026), PASS = validation recall 0.8 after
         the filter (ppr8; the RIS text report vcf_chr1.eval_bed.txt, not vcf_chr1_postpon = post-PoN best-F1 threshold).
         Curve: the own-chr1 re-prediction of the checkpoint in ris_ppr8_20261026 (HPRC v2.1, all seeds:
         <run>/vcf_chr1_<SET>/), else the RIS report copied from RIS (<run>/vcf_chr1/eval_bed/report.json, GRCh38-only /
         HPRC v1.1), else only the PASS and best points of the RIS text report. Every report is checked against the RIS
         text (PASS and best F1) in check.md.
- four-line (HPRC v2.1; HCC1395T, HCC1954T, H1437T, H2009T; HG008T and COLO829T unseen):
  SNV multicancer4 (one run, normfix), INDEL unseen_s{0,1,2} (ppr8).
- six-tumor (HPRC v2.1; the four lines + HG008T + COLO829T chr2-22):
  SNV multicancer6 (one run, normfix), INDEL release_s{0,1,2} (ppr8).
- comparators: DeepSomatic-TO, ClairS-TO (E8 score): curve = every record except RefCall / NoCall ranked by QUAL
  ('all' before, 'all_pon' after our filter); operating point = tool PASS ('pass' / 'pass_pon').
Output: $DATA/curves.json (every point), and check.md here (INDEL single-tumor curves against the RIS text reports).
Run: python collect.py   (login node, ~1 min)
"""
import json
import re
from pathlib import Path

HERE = Path(__file__).resolve().parent
DATA = Path("/scratch/jshen/data/reviewer_experiments_20261013/E27_pr_curves")
RUNS = Path("/scratch/jshen/data/pansoma_net_v2_runs")
E8 = Path("/scratch/jshen/data/reviewer_experiments_20261013/E8_comparators/score")
MAIN = HERE.parent.parent / "ablation" / "main_runs.json"
RISTXT = Path("/scratch/jshen/Github/Pansoma/analysis/indel_ppr8_9919d4a_20261026/runs")
RISCKPT = RUNS / "ris_ppr8_20261026"

TUMORS = ("HG008T", "COLO829T", "HCC1395T")
PLATFORMS = ("Illumina", "PacBio", "ONT")
REFS = (("linear", "GRCh38-only"), ("d9", "HPRC v1.1"), ("d46", "HPRC v2.1"))


def dataset(tumor, p):
    return "COLO829T_fiberseq" if (tumor, p) == ("COLO829T", "PacBio") else f"{tumor}_{p}"


def short(tumor, p):
    """Run-name prefix of the single-tumor models: HG008 / COLO829 and fiberseq."""
    return f"{tumor[:-1]}_{'fiberseq' if (tumor, p) == ('COLO829T', 'PacBio') else p}"


def pts(roc):
    return [[r.get("p_somatic"), r["precision"], r["recall"], r["f1"]] for r in roc]


def op(a):
    return None if a is None else dict(P=a["precision"], R=a["recall"], F1=a["f1"], tp=a.get("tp_baseline"), fp=a.get("fp"),
                                       fn=a.get("fn"))


def pansoma_report(f, seed, family, ref, threshold_override=None):
    r = json.loads(Path(f).read_text())
    out = dict(family=family, ref=ref, seed=seed, report=str(f), threshold=r.get("threshold"), unplaced_pass=r.get("unplaced_pass"))
    for stage in ("raw", "pon"):
        out[stage] = dict(curve=pts(r[stage]["curve_roc"]), op=op(r[stage]["at_threshold"]), best=op(r[stage]["best"]),
                          best_p=r[stage]["best"]["p_somatic"] if r[stage]["best"] else None, ceiling=r[stage]["ceiling"])
    return out


TXT_ROW = re.compile(r"^(raw|pon)\s+[0-9,]+\s+[0-9,]+\s+[0-9,]+\s+[0-9,]+\s+([0-9.]+)\s+([0-9.]+)\s+([0-9.]+)\s+\|\s+([0-9.]+)\s+([0-9.]+)\s+"
                     r"([0-9.]+)\s+([0-9.]+)")


def ris_txt(run):
    """RIS ppr8 text report (PASS = ppr8 threshold): {stage: op, best, best_p}."""
    f = RISTXT / run / "vcf_chr1.eval_bed.txt"
    if not f.exists():
        return None
    out = {}
    for line in f.read_text().splitlines():
        m = TXT_ROW.match(line)
        if m:
            s, P, R, F, bF, bp, bP, bR = m.groups()
            out[s] = dict(op=dict(P=float(P), R=float(R), F1=float(F)), best=dict(P=float(bP), R=float(bR), F1=float(bF)),
                          best_p=float(bp))
    return dict(report=str(f), **out)


def main():
    DATA.mkdir(parents=True, exist_ok=True)
    main_runs = json.loads(MAIN.read_text())["train"]
    panels, missing, checks = {}, [], []
    for tumor in TUMORS:
        for p in PLATFORMS:
            ds = dataset(tumor, p)
            for kind in ("SNV", "INDEL"):
                series = []
                # single-tumor, three references
                if tumor != "HCC1395T":
                    for ref, _ in REFS:
                        for seed in (0, 1, 2):
                            if kind == "SNV":
                                key = f"{tumor[:-1]}|{'fiberseq' if (tumor, p) == ('COLO829T', 'PacBio') else p}|{ref}"
                                row = main_runs[key][seed]
                                f = RUNS / row["src"] / row["run"] / "vcf_chr1" / "eval_bed" / "report.json"
                                if f.exists():
                                    series.append(pansoma_report(f, seed, "single", ref))
                                else:
                                    missing.append(str(f))
                            else:
                                run = f"{short(tumor, p)}_INDEL_{ref}_ppr8" + ("" if seed == 0 else f"_s{seed}")
                                t = ris_txt(run)
                                if t is None:
                                    missing.append(f"{RISTXT}/{run}"); continue
                                local = RISCKPT / run / f"vcf_chr1_{ds}_d46" / "eval_bed" / "report.json"
                                copied = RISCKPT / run / "vcf_chr1" / "eval_bed" / "report.json"
                                f = local if (ref == "d46" and local.exists()) else copied
                                if f.exists():
                                    s = pansoma_report(f, seed, "single", ref)
                                    s["source"] = "local re-prediction" if f == local else "RIS report"
                                    for stage in ("raw", "pon"):
                                        checks.append(dict(run=run, stage=stage, source=s["source"],
                                                           op=s[stage]["op"]["F1"], ris_op=t[stage]["op"]["F1"],
                                                           best=s[stage]["best"]["F1"], ris_best=t[stage]["best"]["F1"],
                                                           best_p=s[stage]["best_p"], ris_best_p=t[stage]["best_p"]))
                                        # PASS point = the RIS score (the letter's numbers); a local re-prediction can
                                        # differ by a call at the threshold (GPU rounding), kept as op_local
                                        s[stage]["op_local"], s[stage]["op"] = s[stage]["op"], t[stage]["op"]
                                else:
                                    s = dict(family="single", ref=ref, seed=seed, report=t["report"], op_source=t["report"],
                                             **{stage: dict(curve=None, op=t[stage]["op"], best=t[stage]["best"],
                                                            best_p=t[stage]["best_p"]) for stage in ("raw", "pon")})
                                series.append(s)
                # multi-tumor models (HPRC v2.1)
                for fam, snv_run, indel_proto in (("four-line", "multicancer4", "unseen"), ("six-tumor", "multicancer6", "release")):
                    if kind == "SNV":
                        f = RUNS / f"{snv_run}_{p}_SNV" / f"vcf_chr1_{ds}_d46_normfix" / "eval_bed" / "report.json"
                        if f.exists():
                            series.append(pansoma_report(f, 0, fam, "d46"))
                        else:
                            missing.append(str(f))
                    else:
                        for seed in (0, 1, 2):
                            f = RUNS / f"{indel_proto}_s{seed}_{p}_INDEL" / f"vcf_ppr8_chr1_{ds}_d46" / "eval_bed" / "report.json"
                            if f.exists():
                                series.append(pansoma_report(f, seed, fam, "d46"))
                            else:
                                missing.append(str(f))
                # comparators
                for caller in ("deepsomatic", "clairs_to"):
                    f = E8 / caller / ds / kind / "report.json"
                    if not f.exists():
                        missing.append(str(f)); continue
                    e = json.loads(f.read_text())["eval_bed"]
                    series.append(dict(family=caller, ref="GRCh38", seed=0, report=str(f),
                                       raw=dict(curve=pts(e["all"]["curve_roc"]), op=op(e["pass"]["at_threshold"]),
                                                best=op(e["all"]["best"]), ceiling=e["all"]["ceiling"]),
                                       pon=dict(curve=pts(e["all_pon"]["curve_roc"]), op=op(e["pass_pon"]["at_threshold"]),
                                                best=op(e["all_pon"]["best"]), ceiling=e["all_pon"]["ceiling"])))
                panels[f"{tumor}|{p}|{kind}"] = series
    (DATA / "curves.json").write_text(json.dumps(dict(panels=panels, missing=missing)))
    L = ["# Single-tumor INDEL curves against the RIS text reports (vcf_chr1.eval_bed.txt): PASS F1 and best F1 on the chr1 curve",
         "PASS points in curves.json and the figures are the RIS scores; the local PASS F1 is shown for the check.",
         "", "| run | stage | curve source | PASS F1 | RIS PASS F1 | best F1 (p) | RIS best F1 (p) |", "|---|---|---|---|---|---|---|"]
    L += [f"| {c['run']} | {c['stage']} | {c['source']} | {c['op']:.3f} | {c['ris_op']:.3f} | {c['best']:.3f} ({c['best_p']:.4f}) | "
          f"{c['ris_best']:.3f} ({c['ris_best_p']:.4f}) |" for c in checks]
    if checks:
        d_op = max(abs(c["op"] - c["ris_op"]) for c in checks); d_best = max(abs(c["best"] - c["ris_best"]) for c in checks)
        L += ["", f"{len(checks) // 2} runs; largest difference: PASS F1 {d_op:.3f}, best F1 {d_best:.3f}"]
    nocurve = [f"{k} {x['family']} {x['ref']} s{x['seed']}" for k, v in panels.items() for x in v if x["pon"]["curve"] is None]
    L += ["", f"Series without a curve (PASS and best points only): {len(nocurve)}", ""] + [f"- {m}" for m in nocurve]
    L += ["", f"Missing inputs: {len(missing)}", ""] + [f"- {m}" for m in missing]
    (HERE / "check.md").write_text("\n".join(L) + "\n")
    n = sum(len(v) for v in panels.values())
    print(f"{len(panels)} panels, {n} series, {len(missing)} missing; check.md written")


if __name__ == "__main__":
    main()
