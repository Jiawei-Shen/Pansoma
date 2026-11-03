"""Tables of the ppr8 INDEL runs of submit.sh, from whatever has finished. Standard library only (login node):

    PANSOMA_RUNS=/storage2/fs1/epigenome/Active/jiawei.s/pansoma_net_v2_runs python3 collect.py

chr1 after the PoN, in the benchmark BED (COLO829T: the full SMaHT union). Every cell: PASS F1 at the checkpoint's call
rule (post-PoN validation recall 0.8 of the callable positives, `vcf_chr1*`), the post-PoN validation best-F1 threshold
(`vcf_chr1*_postpon`) and the best F1 of the chr1 curve (an upper bound, chosen on the test chromosome). The earlier
runs beside them: the `--norm dataset` INDEL runs of 2026-10-19/20 (`*_INDEL_<arm>_normbq`, PASS at the validation
recall-0.9 threshold, 71cb51a INDEL PoN rule) and the 7fc7396 purity joint models (`vcf_chr1_<SET>_pon71cb51a`).
Writes tables.md and summary.json next to this script, and copies each run's train.log and rtg reports to runs/<run>/.
"""
import json
import os
import re
import shutil
import statistics
from pathlib import Path

RUNS = Path(os.environ.get("PANSOMA_RUNS", "/storage2/fs1/epigenome/Active/jiawei.s/pansoma_net_v2_runs"))
HERE = Path(__file__).resolve().parent
SETS = [("HG008", "Illumina"), ("HG008", "PacBio"), ("HG008", "ONT"),
        ("COLO829", "Illumina"), ("COLO829", "fiberseq"), ("COLO829", "ONT")]
ARMS = ["d46", "d9", "linear"]
ARM_NAME = {"d46": "HPRC v2.1 d46", "d9": "HPRC v1.1 d9", "linear": "GRCh38-only"}
SEEDS = ["", "_s1", "_s2"]
PURITIES = ["", "_P75", "_P50", "_P25"]
CROSS = {"Illumina": "Illumina", "PacBio": "fiberseq", "ONT": "ONT"}  # HG008T platform -> COLO829T platform


def set_of(s, p, arm):
    return f"{s}T_{p}" + {"d9": "", "d46": "_d46", "linear": "_linear"}[arm]


def report(run, vcf, bed="eval_bed"):
    f = RUNS / run / vcf / bed / "report.json"
    return json.loads(f.read_text()) if f.exists() else None


def score(run, vcf="vcf_chr1"):
    """{f1, p, r, tp, fp, threshold, best, postpon} of a run's call set (and its _postpon twin), or None."""
    r = report(run, vcf)
    if r is None:
        return None
    t, b = r["pon"]["at_threshold"], r["pon"]["best"]
    out = dict(f1=t["f1"], p=t["precision"], r=t["recall"], tp=t["tp_baseline"], fp=t["fp"], threshold=r["threshold"],
               best=b["f1"], best_threshold=b["p_somatic"], truth=r["pon"]["baseline"], ceiling=r["pon"]["ceiling"],
               truth_tagged=r["truth_pon_tags"]["tagged"])
    nb = report(run, vcf, "eval_nobed")
    out["f1_nobed"] = nb["pon"]["at_threshold"]["f1"] if nb else None
    pp = report(run, f"{vcf}_postpon")
    out["postpon"] = pp["pon"]["at_threshold"]["f1"] if pp else None
    out["postpon_threshold"] = pp["threshold"] if pp else None
    out["postpon_p"] = pp["pon"]["at_threshold"]["precision"] if pp else None
    out["postpon_r"] = pp["pon"]["at_threshold"]["recall"] if pp else None
    return out


def earlier(run, vcf="vcf_chr1"):
    """The 71cb51a-rule call set of an earlier run (<vcf>_pon71cb51a when it exists)."""
    for v in (f"{vcf}_pon71cb51a", vcf):
        r = report(run, v)
        if r is not None:
            t = r["pon"]["at_threshold"]
            return dict(f1=t["f1"], p=t["precision"], r=t["recall"], best=r["pon"]["best"]["f1"], vcf=v)
    return None


def best_epoch(run):
    m = re.search(r"best epoch (\d+)", (RUNS / run / "train.log").read_text()) if (RUNS / run / "train.log").exists() \
        else None
    return int(m.group(1)) if m else None


def f3(x):
    return "-" if x is None else f"{x:.3f}"


def d3(a, b):
    return "-" if a is None or b is None else f"{b - a:+.3f}"


def cell(s):
    """PASS F1 (P / R) | post-PoN best-F1 threshold F1 | best F1"""
    if s is None:
        return "-", "-", "-"
    return f"{s['f1']:.3f} ({s['p']:.3f} / {s['r']:.3f})", f3(s["postpon"]), f3(s["best"])


def seeds_summary(vals):
    vals = [v for v in vals if v is not None]
    if not vals:
        return None, "-"
    m = statistics.mean(vals)
    return m, (f"{m:.3f} [{min(vals):.3f}-{max(vals):.3f}]" if len(vals) > 1 else f"{m:.3f}") + \
        ("" if len(vals) == 3 else f" (n={len(vals)})")


def copy_logs(run, vcfs):
    d = HERE / "runs" / run
    if (RUNS / run / "train.log").exists():
        d.mkdir(parents=True, exist_ok=True)
        shutil.copy2(RUNS / run / "train.log", d / "train.log")
    for v in vcfs:
        for bed in ("eval_bed", "eval_nobed"):
            f = RUNS / run / v / bed / "report.txt"
            if f.exists():
                d.mkdir(parents=True, exist_ok=True)
                shutil.copy2(f, d / f"{v}.{bed}.txt")


def main():
    out = dict(single={}, nopc={}, joint={}, cross={})
    L = ["# ppr8 INDEL runs (9919d4a INDEL pipeline): results\n",
         "(`collect.py`; chr1 after the PoN, in the benchmark BED. PASS: the checkpoint's call rule, post-PoN validation "
         "recall 0.8 of the callable positives; post-PoN best: the post-PoN validation best-F1 threshold; best: the chr1 "
         "curve's best point (upper bound, chosen on chr1). Earlier: the `--norm dataset` runs `*_INDEL_<arm>_normbq` "
         "(PASS at validation recall 0.9; 71cb51a PoN rule). '-' = not finished.)\n"]

    # 1. seed 0 per set and arm against the earlier run
    L += ["\n## Seed 0 per set and reference, against the earlier run\n",
          "| set | reference | best epoch | PASS F1 (P / R) | post-PoN best | best | earlier PASS (P / R) | earlier best "
          "| PASS - earlier PASS | run |", "|---|---|---|---|---|---|---|---|---|---|"]
    for s, p in SETS:
        for arm in ARMS:
            run = f"{s}_{p}_INDEL_{arm}_ppr8"
            sc, old = score(run), earlier(f"{s}_{p}_INDEL_{arm}_normbq")
            a, b, c = cell(sc)
            o = "-" if old is None else f"{old['f1']:.3f} ({old['p']:.3f} / {old['r']:.3f})"
            L.append(f"| {s}T {p} | {arm} | {best_epoch(run) or '-'} | {a} | {b} | {c} | {o} | "
                     f"{'-' if old is None else f3(old['best'])} | {d3(old and old['f1'], sc and sc['f1'])} | `{run}` |")
            out["single"].setdefault(f"{s}_{p}", {}).setdefault(arm, {})["s0"] = sc
            out["single"][f"{s}_{p}"][arm]["earlier"] = old
            copy_logs(run, ["vcf_chr1", "vcf_chr1_postpon"])

    # 2. three seeds per set and arm, and the arm differences of the three-seed means
    L += ["\n## Three seeds per set and reference (mean [min-max])\n",
          "| set | reference | PASS F1 | post-PoN best | best | false calls at PASS | seed PASS F1 |",
          "|---|---|---|---|---|---|---|"]
    means = {}
    for s, p in SETS:
        for arm in ARMS:
            scs = [score(f"{s}_{p}_INDEL_{arm}_ppr8{x}") for x in SEEDS]
            for x, sc in zip(SEEDS, scs):
                out["single"][f"{s}_{p}"][arm]["s" + (x[2:] if x else "0")] = sc
                if x:
                    copy_logs(f"{s}_{p}_INDEL_{arm}_ppr8{x}", ["vcf_chr1", "vcf_chr1_postpon"])
            mp, tp = seeds_summary([sc and sc["f1"] for sc in scs])
            mq, tq = seeds_summary([sc and sc["postpon"] for sc in scs])
            mb, tb = seeds_summary([sc and sc["best"] for sc in scs])
            fps = [sc["fp"] for sc in scs if sc]
            tfp = f"{statistics.mean(fps):.0f} [{min(fps)}-{max(fps)}]" if fps else "-"
            means[(s, p, arm)] = (mp, mq, mb)
            L.append(f"| {s}T {p} | {arm} | {tp} | {tq} | {tb} | {tfp} | "
                     f"{' / '.join(f3(sc and sc['f1']) for sc in scs)} |")
    L += ["\nDifferences of the three-seed means, PASS F1 (post-PoN best; best):\n",
          "| set | d46 - d9 | d46 - GRCh38-only | d9 - GRCh38-only |", "|---|---|---|---|"]
    for s, p in SETS:
        row = []
        for a1, a2 in (("d9", "d46"), ("linear", "d46"), ("linear", "d9")):
            x, y = means[(s, p, a1)], means[(s, p, a2)]
            row.append(f"{d3(x[0], y[0])} ({d3(x[1], y[1])}; {d3(x[2], y[2])})")
        L.append(f"| {s}T {p} | " + " | ".join(row) + " |")
    out["seed_means"] = {f"{s}_{p}_{a}": dict(zip(("pass", "postpon", "best"), v)) for (s, p, a), v in means.items()}

    # 3. no ch6 on d46
    L += ["\n## No ch6 on d46 (seed 0) against the three with-ch6 seeds\n",
          "| set | no ch6 PASS (P / R) | post-PoN best | best | with ch6 PASS seeds | no ch6 - with ch6 seed 0 PASS | "
          "below every with-ch6 seed (PASS / best) | earlier no-ch6 PASS / best |",
          "|---|---|---|---|---|---|---|---|"]
    for s, p in SETS:
        run = f"{s}_{p}_INDEL_nopc_d46_ppr8"
        sc = score(run)
        withs = [score(f"{s}_{p}_INDEL_d46_ppr8{x}") for x in SEEDS]
        a, b, c = cell(sc)
        below = "-"
        if sc and all(withs):
            below = ("yes" if sc["f1"] < min(w["f1"] for w in withs) else "no") + "/" + \
                    ("yes" if sc["best"] < min(w["best"] for w in withs) else "no")
        old = earlier(f"{s}_{p}_INDEL_nopc_d46_normbq")
        L.append(f"| {s}T {p} | {a} | {b} | {c} | {' / '.join(f3(w and w['f1']) for w in withs)} | "
                 f"{d3(withs[0] and withs[0]['f1'], sc and sc['f1'])} | {below} | "
                 f"{'-' if old is None else f3(old['f1']) + ' / ' + f3(old['best'])} |")
        out["nopc"][f"{s}_{p}"] = dict(nopc=sc, earlier=old)
        copy_logs(run, ["vcf_chr1", "vcf_chr1_postpon"])

    # 4. COLO829T purity joint models (d9)
    L += ["\n## COLO829T purity: one model per platform on P100 + P75 + P50 + P25 (d9)\n",
          "| platform | purity | PASS F1 (P / R) | post-PoN best | best | 7fc7396 joint PASS / best | single-set d9 seeds "
          "PASS (P100) |", "|---|---|---|---|---|---|---|"]
    for p in ("Illumina", "fiberseq", "ONT"):
        run = f"COLO829_{p}_INDEL_joint_ppr8"
        vcfs = []
        for d in PURITIES:
            sset = f"COLO829T_{p}{d}"
            sc = score(run, f"vcf_chr1_{sset}")
            old = earlier(f"COLO829_{p}_INDEL_joint_7fc7396", f"vcf_chr1_{sset}")
            single = "-"
            if not d:
                single = " / ".join(f3((score(f"COLO829_{p}_INDEL_d9_ppr8{x}") or {}).get("f1")) for x in SEEDS)
            a, b, c = cell(sc)
            L.append(f"| {p} | {('P' + d[2:]) if d else 'P100'} | {a} | {b} | {c} | "
                     f"{'-' if old is None else f3(old['f1']) + ' / ' + f3(old['best'])} | {single} |")
            out["joint"][sset] = dict(ppr8=sc, earlier=old)
            vcfs += [f"vcf_chr1_{sset}", f"vcf_chr1_{sset}_postpon"]
        copy_logs(run, vcfs)
    L.append("\nOne PASS threshold per model (pooled validation of the four purities). The Illumina and ONT dilutions "
             "carry the read-source BQ cue (`../benchmark_7fc7396_20261012/colo829t_purity_joint_7fc7396_20261012`): "
             "their scores are likely optimistic.")

    # 5. cross-tumor d46
    L += ["\n## Cross-tumor (d46): a tumor's models on chr1 of the other tumor, against the target's own models\n",
          "| model (trained on) | tested on | PASS F1 seeds 0 / 1 / 2 | post-PoN best | best | target's own models PASS "
          "(mean) | cross - own (PASS means) |", "|---|---|---|---|---|---|---|"]
    for hp, cp in CROSS.items():
        for src, dst, own in ((f"HG008_{hp}", f"COLO829T_{cp}_d46", f"COLO829_{cp}"),
                              (f"COLO829_{cp}", f"HG008T_{hp}_d46", f"HG008_{hp}")):
            scs = [score(f"{src}_INDEL_d46_ppr8{x}", f"vcf_chr1_{dst}") for x in SEEDS]
            for x in SEEDS:
                copy_logs(f"{src}_INDEL_d46_ppr8{x}", [f"vcf_chr1_{dst}", f"vcf_chr1_{dst}_postpon"])
            mp, _ = seeds_summary([sc and sc["f1"] for sc in scs])
            _, tq = seeds_summary([sc and sc["postpon"] for sc in scs])
            _, tb = seeds_summary([sc and sc["best"] for sc in scs])
            mo, to = seeds_summary([(score(f"{own}_INDEL_d46_ppr8{x}") or {}).get("f1") for x in SEEDS])
            L.append(f"| {src.replace('_', 'T ', 1)} | {dst} | "
                     f"{' / '.join(f3(sc and sc['f1']) for sc in scs)} | {tq} | {tb} | {to} | {d3(mo, mp)} |")
            out["cross"][f"{src}->{dst}"] = scs

    # chr1 details
    L += ["\n## chr1 details (seed 0)\n",
          "| run | PASS threshold | post-PoN best threshold | truth (BED) | truth PoN-tagged | recall ceiling after the PoN | "
          "TP / FP at PASS | PASS F1 without the BED |", "|---|---|---|---|---|---|---|---|"]
    for s, p in SETS:
        for arm in ARMS:
            run = f"{s}_{p}_INDEL_{arm}_ppr8"
            sc = score(run)
            if sc:
                pt = "-" if sc["postpon_threshold"] is None else f"{sc['postpon_threshold']:.4f}"
                L.append(f"| `{run}` | {sc['threshold']:.4f} | {pt} | {sc['truth']} | {sc['truth_tagged']} | "
                         f"{sc['ceiling']:.3f} | {sc['tp']} / {sc['fp']} | {f3(sc['f1_nobed'])} |")
    (HERE / "tables.md").write_text("\n".join(L) + "\n")
    (HERE / "summary.json").write_text(json.dumps(out, indent=1) + "\n")
    print("\n".join(L))


if __name__ == "__main__":
    main()
