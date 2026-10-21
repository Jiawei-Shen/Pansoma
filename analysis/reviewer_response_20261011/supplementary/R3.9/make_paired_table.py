"""R3.9 Supplementary Table b: paired 1-Mb block-bootstrap 95% intervals (2,000 resamples, seed 20261015) of the
differences in chr1 precision, recall and F1 after population filtering (benchmark regions), from
experiments/E16_bootstrap/paired.tsv: HPRC graph minus the GRCh38-only graph (SNVs), and Pansoma on the HPRC v2.1 graph
minus DeepSomatic-TO 1.9.0 and ClairS-TO (0.4.2 on HG008T, 0.4.1 on COLO829T; the user's existing tumor-only runs,
re-scored by E8). Pansoma call sets: the first-seed in-sample models (COLO829T Illumina v1.1: epoch 6, E16's ARM_RUN). Writes paired_table.tsv and paired_table.md.
usage: python3 make_paired_table.py
"""
import csv
from pathlib import Path

HERE = Path(__file__).resolve().parent
SRC = HERE.parents[1] / "experiments" / "E16_bootstrap" / "paired.tsv"
GRAPH = {"d9": "HPRC v1.1 (d9)", "d46": "HPRC v2.1 (d46)"}


def arm(uid):
    run = uid.split("/")[0]
    return "d46" if "_d46_" in run else "linear" if "_linear_" in run else "d9"


def fmt(r, m):
    return f"{float(r[m]):+.3f} [{float(r[m + '_lo']):+.3f}, {float(r[m + '_hi']):+.3f}]"


def main():
    rows = [r for r in csv.DictReader(open(SRC), delimiter="\t") if r["mode"] == "bed" and r["block_mb"] == "1"]
    out = []
    for r in rows:
        if r["type"] == "arms" and " SNV " in r["A_label"] and r["B_label"].endswith("linear"):
            out.append(("HPRC - GRCh38-only", r["A_label"].split(" ")[0], "SNV", GRAPH[arm(r["A"])], "GRCh38-only graph", r))
        elif (r["type"] == "Pansoma vs comparator" and r["A"].endswith("/vcf_chr1") and arm(r["A"]) == "d46"
              and "_joint_" not in r["A"]):
            lab = r["B_label"].split(" on ")[0]
            s, k = r["B_label"].split(" on ")[1].rsplit(" ", 1)
            out.append(("Pansoma - comparator", s, k, GRAPH["d46"], lab, r))
    out.sort(key=lambda x: (x[0] != "HPRC - GRCh38-only", x[2] != "SNV", x[1], x[3], x[4]))
    with open(HERE / "paired_table.tsv", "w") as f:
        w = csv.writer(f, delimiter="\t", lineterminator="\n")
        w.writerow(["comparison", "set", "kind", "A", "B", "dP", "dP_lo", "dP_hi", "dR", "dR_lo", "dR_hi", "dF1", "dF1_lo", "dF1_hi"])
        for c, s, k, a, b, r in out:
            w.writerow([c, s, k, a, b] + [r[x] for x in ("dP", "dP_lo", "dP_hi", "dR", "dR_lo", "dR_hi", "dF1", "dF1_lo", "dF1_hi")])
    md = ["| comparison | dataset | kind | A | B | Δprecision [95% CI] | Δrecall [95% CI] | ΔF1 [95% CI] |",
          "|---|---|---|---|---|---|---|---|"]
    for c, s, k, a, b, r in out:
        md.append(f"| {c} | {s.replace('_', ' ')} | {k} | {a} | {b} | {fmt(r, 'dP')} | {fmt(r, 'dR')} | {fmt(r, 'dF1')} |")
    (HERE / "paired_table.md").write_text("\n".join(md) + "\n")
    print(len(out), "rows")


if __name__ == "__main__":
    main()
