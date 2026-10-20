"""R3.9 Supplementary Table: chr1 SNV precision, recall and F1 after population filtering (benchmark regions, PASS at the
recall-0.9 validation threshold) of the in-sample models on the HPRC v1.1 (d9) and v2.1 (d46) graphs, three training
seeds each, with mean +- SD (sample SD over the three seeds). Source: the RIS summary of 0a07f91,
analysis/benchmark_7fc7396_20261012/reviewer_gpu_7fc7396_20261016/summary.json, `normbq_seeds` (models with base
quality standardised per dataset over all its chromosomes, f344bca). COLO829T Illumina on d9: epochs 6 / 5 / 6 of
seeds 0 / 1 / 2 (the user, 2026-10-20; RIS b662f00), with the validation-selected epoch 8 in an extra row.
Writes seed_table.tsv and seed_table.md here.
usage: python3 make_seed_table.py
"""
import json
import statistics as st
from pathlib import Path

HERE = Path(__file__).resolve().parent
SRC = HERE.parents[2] / "benchmark_7fc7396_20261012" / "reviewer_gpu_7fc7396_20261016" / "summary.json"
ARM = {"d9": "HPRC v1.1 (d9)", "d46": "HPRC v2.1 (d46)",
       "d9 validation-selected (epoch 8)": "HPRC v1.1 (d9), validation-selected epoch 8"}


def main():
    rows = [r for r in json.loads(SRC.read_text())["normbq_seeds"] if r["arm"] in ARM]
    order = list(dict.fromkeys(r["set"] for r in rows))
    rows.sort(key=lambda r: (order.index(r["set"]), list(ARM).index(r["arm"])))
    tsv = ["tumor\tplatform\tgraph\tmetric\tseed0\tseed1\tseed2\tmean\tsd"]
    md = ["| tumor | platform | graph | precision: seeds 0 / 1 / 2 (mean ± SD) | recall | F1 |", "|---|---|---|---|---|---|"]
    for r in rows:
        tumor, plat = r["set"].split(" ", 1)
        sc = r["scores"]
        assert all(sc), r["set"]
        cells = []
        for key, name in (("p", "precision"), ("r", "recall"), ("pass", "F1")):
            v = [s[key] for s in sc]
            m, sd = st.mean(v), st.stdev(v)
            tsv.append("\t".join([tumor, plat, ARM[r["arm"]], name] + [f"{x:.4f}" for x in v] + [f"{m:.4f}", f"{sd:.4f}"]))
            cells.append(f"{' / '.join(f'{x:.3f}' for x in v)} ({m:.3f} ± {sd:.3f})")
        lab = ARM[r['arm']] + (" at epochs 6 / 5 / 6" if r["arm"] == "d9" and r["set"] == "COLO829T Illumina" else "")
        md.append(f"| {tumor} | {plat} | {lab} | " + " | ".join(cells) + " |")
    (HERE / "seed_table.tsv").write_text("\n".join(tsv) + "\n")
    (HERE / "seed_table.md").write_text("\n".join(md) + "\n")
    print("\n".join(md))


if __name__ == "__main__":
    main()
