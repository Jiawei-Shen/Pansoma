"""Response Figure 1.2 (ver4, 2026-11-23): graph context outside the candidate node masked at inference, on the
six-tumor current SNV models.

Models: snv_newpipe_20261109/multicancer6_<P>_SNV_newpipe (P = Illumina, PacBio, ONT; seed 0, best.pth), HPRC v2.1.
Test sets (chr1, held out of training): HG008T_<P>_d46, COLO829T_{Illumina,fiberseq}_d46_snvaf05, COLO829T_ONT_d46,
HCC1395T_{Illumina,PacBio}_d46_snvaf05, HCC1395T_ONT_d46.
Full tensor  = the existing <model>/vcf_chr1_<set>/eval_bed/report.json (not recomputed).
Own node only = ver4_masking_20261123/multicancer6_<P>_SNV_newpipe__mask/vcf_chr1_<set>/eval_bed/report.json: the same
  checkpoint, sets, site features, BQ statistics, call threshold (rule 2, read from the checkpoint) and VCF / scoring
  script (tmp/snv_sitefeat_20261103/vcf_code.sh); predict --mask-off-node sets every cell outside the candidate's own
  node visit to 0 in all 8 channels (the site columns are kept); tmp/ver4_20261123/masking.
Numbers: chr1 SNV, benchmark regions (eval_bed), PASS calls; "after" = after the population filter (report pon.at_threshold),
"before" = before it (raw.at_threshold). TP = truth alleles found (tp_baseline), FN = truth missed, FP = false calls;
P = tp_call / (tp_call + FP), R = TP / (TP + FN). P / R / F1 are the report's own values (4 decimals, checked here
against the counts) printed with 3 decimals, as in the other tables of the current models (latest_results_20261123,
supplementary/R1.6 points files); the differences full - own node are taken on the report's 4-decimal values and then
rounded to 3 decimals (results.json also keeps the count-based unrounded values, *_exact).
Check: the patched code without the mask flag (__nomask, HG008T PacBio) must reproduce the existing report.
Writes results.md / results.json next to this script. usage: python collect.py
"""
import json
from pathlib import Path

R = Path("/scratch/jshen/data/pansoma_net_v2_runs")
FULL = R / "snv_newpipe_20261109"
MASK = R / "ver4_masking_20261123"
HERE = Path(__file__).resolve().parent
# tumor, platform label, model platform, test set
ROWS = [("HG008T", "Illumina", "Illumina", "HG008T_Illumina_d46"),
        ("HG008T", "PacBio HiFi", "PacBio", "HG008T_PacBio_d46"),
        ("HG008T", "ONT", "ONT", "HG008T_ONT_d46"),
        ("COLO829T", "Illumina", "Illumina", "COLO829T_Illumina_d46_snvaf05"),
        ("COLO829T", "Fiber-seq", "PacBio", "COLO829T_fiberseq_d46_snvaf05"),
        ("COLO829T", "ONT", "ONT", "COLO829T_ONT_d46"),
        ("HCC1395T", "Illumina", "Illumina", "HCC1395T_Illumina_d46_snvaf05"),
        ("HCC1395T", "PacBio HiFi", "PacBio", "HCC1395T_PacBio_d46_snvaf05"),
        ("HCC1395T", "ONT", "ONT", "HCC1395T_ONT_d46")]


def prf(t):
    p = t["tp_call"] / (t["tp_call"] + t["fp"]) if t["tp_call"] + t["fp"] else 0.0
    r = t["tp_baseline"] / (t["tp_baseline"] + t["fn"])
    f = 2 * p * r / (p + r) if p + r else 0.0
    for name, v in (("precision", p), ("recall", r), ("f1", f)):  # the report's own (4-decimal) values
        assert abs(v - t[name]) < 6e-5, (name, v, t[name])
    return (t["precision"], t["recall"], t["f1"]), (p, r, f)


def arm(d):
    f = d / "eval_bed/report.json"
    if not f.exists():
        return None
    j = json.load(open(f))
    a, b = j["pon"]["at_threshold"], j["raw"]["at_threshold"]
    (p, r, f1), exact = prf(a)
    (pb, rb, fb), exact_before = prf(b)
    return dict(P=p, R=r, F1=f1, P_exact=exact[0], R_exact=exact[1], F1_exact=exact[2], F1_before_exact=exact_before[2], TP=a["tp_baseline"], TP_call=a["tp_call"], FP=a["fp"], FN=a["fn"],
                truth=a["tp_baseline"] + a["fn"], P_before=pb, R_before=rb, F1_before=fb, TP_before=b["tp_baseline"],
                FP_before=b["fp"], FN_before=b["fn"], truth_before=b["tp_baseline"] + b["fn"], threshold=j["threshold"],
                truth_vcf=j["truth"], bed=j["bed"], truth_records=j["truth_records"],
                checkpoint=j["steps"]["graph_vcf"]["checkpoint"], call_rule=j["steps"]["graph_vcf"]["call"],
                predictions=j["steps"]["graph_vcf"]["predictions"], report=str(f))


out, problems = [], []
for tumor, plat, mp, s in ROWS:
    v = "vcf_chr1_" + s.replace("_snvaf05", "")
    full = arm(FULL / f"multicancer6_{mp}_SNV_newpipe" / v)
    own = arm(MASK / f"multicancer6_{mp}_SNV_newpipe__mask" / v)
    row = dict(tumor=tumor, platform=plat, model=f"snv_newpipe_20261109/multicancer6_{mp}_SNV_newpipe", set=s, full=full, own_node=own)
    if own:
        for k in ("truth", "truth_before", "threshold", "truth_vcf", "bed", "truth_records", "checkpoint", "call_rule", "predictions"):
            if full[k] != own[k]:
                problems.append(f"{s}: {k} differs between the arms: {full[k]} vs {own[k]}")
        row["delta"] = {k: full[k] - own[k] for k in ("P", "R", "F1", "F1_before")}
        row["delta"] = {k: round(d, 4) for k, d in row["delta"].items()}  # 4-decimal inputs
        row["delta_exact"] = {k: full[k + "_exact"] - own[k + "_exact"] for k in ("P", "R", "F1", "F1_before")}
    out.append(row)

# reproduction check: the patched code without the flag
chk = arm(MASK / "multicancer6_PacBio_SNV_newpipe__nomask/vcf_chr1_HG008T_PacBio_d46")
ref = arm(FULL / "multicancer6_PacBio_SNV_newpipe/vcf_chr1_HG008T_PacBio_d46")
check = None
if chk:
    keys = ("TP", "TP_call", "FP", "FN", "TP_before", "FP_before", "FN_before", "threshold", "predictions")
    check = dict(set="HG008T_PacBio_d46", existing={k: ref[k] for k in keys}, patched_no_mask={k: chk[k] for k in keys},
                 identical=all(ref[k] == chk[k] for k in keys))
    if not check["identical"]:
        problems.append("the patched code without the mask does NOT reproduce the existing HG008T PacBio report")

f3 = lambda x: f"{x:.3f}"  # noqa: E731
d3 = lambda x: f"{x:+.3f}".replace("-", "−")  # noqa: E731
L = ["# Response Figure 1.2 (ver4): graph context outside the candidate node masked at inference", "",
     "Six-tumor current SNV models (`snv_newpipe_20261109/multicancer6_<P>_SNV_newpipe`, seed 0, HPRC v2.1), chr1 SNV, "
     "benchmark regions, PASS calls (rule 2: the checkpoint's threshold at validation callable recall 0.95 after the "
     "population filter). Full tensor = the existing call sets; own node only = the same checkpoint and threshold with "
     "every tensor cell outside the candidate's own node visit set to 0 at inference (the site columns are kept). "
     "Δ = full tensor − own node only (difference of the reports' 4-decimal values, rounded to 3 decimals).", "",
     "## chr1, after the population filter", "",
     "| Tumor (chr1) | Platform | Truth (chr1) | Full tensor P / R / F1 (chr1) | Own node only P / R / F1 (chr1) | ΔP | ΔR | ΔF1 | Full TP / FP / FN | Own node TP / FP / FN |",
     "|---|---|---|---|---|---|---|---|---|---|"]
done = [r for r in out if r["own_node"]]
for r in out:
    a, b = r["full"], r["own_node"]
    if b:
        d = r["delta"]
        L.append(f"| {r['tumor']} | {r['platform']} | {a['truth']:,} | {f3(a['P'])} / {f3(a['R'])} / {f3(a['F1'])} | "
                 f"{f3(b['P'])} / {f3(b['R'])} / {f3(b['F1'])} | {d3(d['P'])} | {d3(d['R'])} | {d3(d['F1'])} | "
                 f"{a['TP']:,} / {a['FP']:,} / {a['FN']:,} | {b['TP']:,} / {b['FP']:,} / {b['FN']:,} |")
    else:
        L.append(f"| {r['tumor']} | {r['platform']} | {a['truth']:,} | {f3(a['P'])} / {f3(a['R'])} / {f3(a['F1'])} | pending | | | | "
                 f"{a['TP']:,} / {a['FP']:,} / {a['FN']:,} | pending |")
L += ["", "## chr1, before the population filter", "",
      "| Tumor (chr1) | Platform | Full tensor P / R / F1 (chr1, before) | Own node only P / R / F1 (chr1, before) | ΔF1 (before) | Full TP / FP / FN (before) | Own node TP / FP / FN (before) |",
      "|---|---|---|---|---|---|---|"]
for r in out:
    a, b = r["full"], r["own_node"]
    if b:
        L.append(f"| {r['tumor']} | {r['platform']} | {f3(a['P_before'])} / {f3(a['R_before'])} / {f3(a['F1_before'])} | "
                 f"{f3(b['P_before'])} / {f3(b['R_before'])} / {f3(b['F1_before'])} | {d3(r['delta']['F1_before'])} | "
                 f"{a['TP_before']:,} / {a['FP_before']:,} / {a['FN_before']:,} | {b['TP_before']:,} / {b['FP_before']:,} / {b['FN_before']:,} |")
    else:
        L.append(f"| {r['tumor']} | {r['platform']} | {f3(a['P_before'])} / {f3(a['R_before'])} / {f3(a['F1_before'])} | pending | | "
                 f"{a['TP_before']:,} / {a['FP_before']:,} / {a['FN_before']:,} | pending |")
summary = {}
if done:
    L += ["", f"## Direction of the change ({len(done)} comparisons, after the population filter)", "",
          "| Metric | Lower with own node only (Δ > 0) | Higher with own node only (Δ < 0) | Equal | Range of Δ | Mean Δ |", "|---|---|---|---|---|---|"]
    for k, name in (("R", "Recall"), ("P", "Precision"), ("F1", "F1"), ("F1_before", "F1 before the filter")):
        ds = [r["delta"][k] for r in done]
        counts = [(r["full"]["TP"], r["own_node"]["TP"]) for r in done] if k == "R" else None
        lower = sum(d > 1e-12 for d in ds); higher = sum(d < -1e-12 for d in ds)
        summary[k] = dict(lower_with_mask=lower, higher_with_mask=higher, equal=len(ds) - lower - higher, min=min(ds), max=max(ds),
                          mean=sum(ds) / len(ds))
        L.append(f"| {name} | {lower} | {higher} | {len(ds) - lower - higher} | {d3(min(ds))} to {d3(max(ds))} | {d3(sum(ds) / len(ds))} |")
    tp = [(r["full"]["TP"], r["own_node"]["TP"]) for r in done]; fp = [(r["full"]["FP"], r["own_node"]["FP"]) for r in done]
    summary["tp_lost"] = [a - b for a, b in tp]; summary["fp_change"] = [b - a for a, b in fp]
    L += ["", f"Truth alleles lost by the mask per comparison (full TP − own-node TP): {', '.join(str(a - b) for a, b in tp)} "
          f"(sum {sum(a - b for a, b in tp):,} of {sum(a for a, _ in tp):,} full-tensor TP). "
          f"False calls, own node − full: {', '.join(f'{b - a:+d}' for a, b in fp)} (sum {sum(b - a for a, b in fp):+,} on {sum(a for a, _ in fp):,})."]
fr = HERE / "mask_fraction.json"
if fr.exists():
    rows = {(r["set"], r["group"]): r for r in json.load(open(fr))}
    L += ["", "## Share of the tensor cells the mask zeroes (chr1 samples, mask_fraction.py)", "",
          "| Set (chr1) | Tensors sampled (all / somatic label) | Tensors with a masked cell (all / somatic) | Non-empty cells zeroed (all / somatic) |", "|---|---|---|---|"]
    for _, _, _, s in ROWS:
        a, b = rows.get((s, "all")), rows.get((s, "somatic"))
        if a and b:
            L.append(f"| {s} | {a['sampled']:,} / {b['sampled']:,} | {a['tensors_with_masked']:.1%} / {b['tensors_with_masked']:.1%} | "
                     f"{a['masked_share']:.1%} / {b['masked_share']:.1%} |")
            for r in out:
                if r["set"] == s:
                    r["mask_fraction"] = dict(all=a, somatic=b)
L += ["", "## Checks", ""]
if check:
    e, n = check["existing"], check["patched_no_mask"]
    L.append(f"- Patched code without the mask flag, HG008T PacBio HiFi chr1 (model PacBio): after the filter TP / FP / FN "
             f"{n['TP']} / {n['FP']} / {n['FN']} vs the existing report {e['TP']} / {e['FP']} / {e['FN']}; before the filter "
             f"{n['TP_before']} / {n['FP_before']} / {n['FN_before']} vs {e['TP_before']} / {e['FP_before']} / {e['FN_before']}; "
             f"threshold {n['threshold']} vs {e['threshold']}: {'identical' if check['identical'] else 'DIFFERENT'}.")
else:
    L.append("- Patched code without the mask flag (HG008T PacBio HiFi): pending.")
L.append(f"- Truth counts, truth VCF, BED, checkpoint, call rule, threshold and number of predictions equal between the two arms in "
         f"{len(done) - len({p.split(':')[0] for p in problems if 'arms' in p})} of {len(done)} finished comparisons"
         + ("." if len(done) == len(ROWS) else f"; {len(ROWS) - len(done)} pending."))
for p in problems:
    L.append(f"- PROBLEM: {p}")
(HERE / "results.md").write_text("\n".join(L) + "\n")
(HERE / "results.json").write_text(json.dumps(dict(rows=out, check=check, summary=summary, problems=problems), indent=1) + "\n")
print("\n".join(L))
