"""R3.9 Supplementary Table c: the seed-matched three-reference comparison behind "35 of 36" (R3.1): chr1 SNV precision,
recall and F1 after population filtering (benchmark regions, PASS at the recall-0.9 validation threshold) of the
in-sample models on the GRCh38-only graph, HPRC v1.1 (d9) and HPRC v2.1 (d46), three seeds each (base and mapping
quality standardised with the training statistics, the setting that has GRCh38-only seeds). Source: the RIS summary
(analysis/benchmark_7fc7396_20261012/reviewer_gpu_7fc7396_20261016/summary.json, `seeds`; COLO829T Illumina v1.1 at
its validation-selected epoch). Writes three_ref_table.tsv and three_ref_table.md here.
usage: python3 make_three_ref_table.py
"""
import json
import statistics as st
from pathlib import Path

HERE = Path(__file__).resolve().parent
SRC = HERE.parents[2] / "benchmark_7fc7396_20261012" / "reviewer_gpu_7fc7396_20261016" / "summary.json"
ARM = {"linear": "GRCh38-only", "d9": "HPRC v1.1 (d9)", "d46": "HPRC v2.1 (d46)"}


def main():
    rows = json.loads(SRC.read_text())["seeds"]
    by = {(r["set"], r["arm"]): r["scores"] for r in rows}
    tsv = ["tumor\tplatform\treference\tmetric\tseed0\tseed1\tseed2\tmean\tsd"]
    md = ["| tumor | platform | reference | precision: seeds 0 / 1 / 2 (mean ± SD) | recall | F1 | HPRC > GRCh38-only, same seed |",
          "|---|---|---|---|---|---|---|"]
    wins = total = 0
    for (sset, arm), sc in by.items():
        tumor, plat = sset.split(" ", 1)
        cells = []
        for key, name in (("p", "precision"), ("r", "recall"), ("pass", "F1")):
            v = [s[key] for s in sc]
            tsv.append("\t".join([tumor, plat, ARM[arm], name] + [f"{x:.4f}" for x in v] + [f"{st.mean(v):.4f}", f"{st.stdev(v):.4f}"]))
            cells.append(f"{' / '.join(f'{x:.3f}' for x in v)} ({st.mean(v):.3f} ± {st.stdev(v):.3f})")
        cmp = "-"
        if arm != "linear":
            lin = by[(sset, "linear")]
            w = sum(a["pass"] > b["pass"] for a, b in zip(sc, lin))
            wins += w; total += 3
            cmp = f"{w} / 3"
        md.append(f"| {tumor} | {plat} | {ARM[arm]} | " + " | ".join(cells) + f" | {cmp} |")
    (HERE / "three_ref_table.tsv").write_text("\n".join(tsv) + "\n")
    (HERE / "three_ref_table.md").write_text("\n".join(md) + "\n")
    print(f"HPRC above GRCh38-only, same seed: {wins} / {total}")


if __name__ == "__main__":
    main()
