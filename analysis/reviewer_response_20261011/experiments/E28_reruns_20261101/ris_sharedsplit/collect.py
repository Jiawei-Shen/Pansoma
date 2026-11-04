"""E28 (3) on RIS: R3.1 shared validation split (chr20-22), 18 SNV models (HG008T / COLO829T x three platforms x
GRCh38-only / HPRC v1.1 / HPRC v2.1, seed 0), trained by submit.sh. Copies each run's chr1 vcfeval report and the
"done; best epoch" line of its train.log here, and tabulates PASS after the PoN in the benchmark BED.

    python3 collect.py   (stdlib only; PANSOMA_RUNS defaults to the RIS runs directory)

Writes runs/<run>/report.json (= <run>/vcf_chr1/eval_bed/report.json), best_epoch.tsv, tables.md, summary.json.
"""
import json
import os
import shutil
from pathlib import Path

HERE = Path(__file__).resolve().parent
RUNS = Path(os.environ.get("PANSOMA_RUNS", "/storage2/fs1/epigenome/Active/jiawei.s/pansoma_net_v2_runs"))
ARMS = (("linear", "GRCh38-only", "_linear"), ("d9", "HPRC v1.1 (d9)", ""), ("d46", "HPRC v2.1 (d46)", "_d46"))


def runs():
    for t in ("HG008", "COLO829"):
        for p in ("Illumina", "PacBio", "ONT"):
            m = "fiberseq" if t == "COLO829" and p == "PacBio" else p
            for arm, label, suffix in ARMS:
                yield f"{t}_{m}_SNV_{arm}_valchr20-22", f"{t}T_{m}{suffix}", f"{t}T {m}", label


def main():
    rows, out = [], {}
    L = ["# E28 (3) on RIS: R3.1 SNV, validation chr20-22 shared by the three references (seed 0)\n",
         "(`collect.py`; chr1, after the PoN, in the benchmark BED. PASS = the checkpoint's call rule (validation "
         "recall 0.9); best = the chr1 curve's best point. Recipe of Table R3.1: 50.4 M, 20 epochs, pow:0.75, lr 1e-4, "
         "4 x H100 x 1024, `--norm train`, `--val-chroms chr20 chr21 chr22`, training chr2-19.)\n",
         "| set | reference | best epoch | PASS P / R / F1 | best F1 | run |", "|---|---|---|---|---|---|"]
    for run, tset, name, label in runs():
        d = RUNS / run
        rep = d / "vcf_chr1" / "eval_bed" / "report.json"
        done = next((l.strip() for l in (d / "train.log").read_text().splitlines() if l.startswith("done; best epoch")),
                    None) if (d / "train.log").exists() else None
        rows.append(f"{run}\t{tset}\t{done or '-'}")
        if rep.exists():
            (HERE / "runs" / run).mkdir(parents=True, exist_ok=True)
            shutil.copy2(rep, HERE / "runs" / run / "report.json")
            r = json.loads(rep.read_text())["pon"]
            a, b = r["at_threshold"], r["best"] or {}
            epoch = done.split("best epoch ")[1].split(" ")[0] if done else "-"
            out[run] = dict(set=tset, best_epoch=epoch, done=done, p=a["precision"], r=a["recall"], f1=a["f1"],
                            best_f1=b.get("f1"))
            L.append(f"| {name} | {label} | {epoch} | {a['precision']:.3f} / {a['recall']:.3f} / **{a['f1']:.3f}** | "
                     f"{b.get('f1', float('nan')):.3f} | `{run}` |")
        else:
            L.append(f"| {name} | {label} | - | - | - | `{run}` (not finished) |")
    (HERE / "best_epoch.tsv").write_text("run\tset\tdone line\n" + "\n".join(rows) + "\n")
    (HERE / "tables.md").write_text("\n".join(L) + "\n")
    (HERE / "summary.json").write_text(json.dumps(out, indent=1) + "\n")
    print("\n".join(L))


if __name__ == "__main__":
    main()
