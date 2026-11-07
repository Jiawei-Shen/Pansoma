"""R3.9 Supplementary Table: per-seed chr1 INDEL precision, recall and F1 of the final INDEL models (pre-registered
batch, PASS at the post-filter validation recall 0.8 = rule ppr8; HPRC v2.1 d46; benchmark regions; after the shared
gnomAD / CoLoRSdb INDEL filter) on HG008T and COLO829T, six-tumor (release) and four-line (unseen) models, with mean +- SD
(sample SD over the three seeds) and the paired 1-Mb chr1 block-bootstrap difference (2,000 resamples) of the seed-mean
F1 minus each comparator's F1, as computed by tmp/indel_gap_20261020/final_20261023/analyze.py.
Source: tmp/indel_gap_20261020/final_20261023/results_ppr8.json (cells; per-seed values from the rtg report.json it lists,
pon.at_threshold). Writes indel_seed_table.tsv and indel_seed_table.md here.
usage: python3 make_indel_seed_table.py
"""
import json
import statistics as st
from pathlib import Path

HERE = Path(__file__).resolve().parent
SRC = HERE.parents[3] / "tmp" / "indel_gap_20261020" / "final_20261023" / "results_ppr8.json"
PLAT = {"HG008T": [("Illumina", "Illumina"), ("PacBio", "PacBio HiFi"), ("ONT", "ONT")],
        "COLO829T": [("Illumina", "Illumina"), ("PacBio", "Fiber-seq"), ("ONT", "ONT")]}
PROTO = [("release", "six-tumor"), ("unseen", "four-line")]


def main():
    cells = {(c["protocol"], c["tumor"], c["platform"]): c for c in json.loads(SRC.read_text())["cells"] if c["chrom"] == "chr1"}
    tsv = ["tumor\tplatform\ttraining_set\tP_seeds\tR_seeds\tF1_seeds\tP_mean\tP_sd\tR_mean\tR_sd\tF1_mean\tF1_sd\t"
           "d_deepsomatic\td_deepsomatic_ci\td_clairs_to\td_clairs_to_ci"]
    md = ["| tumor | platform | training set | precision: seeds 0 / 1 / 2 (mean ± SD) | recall | F1 | "
          "ΔF1 vs DeepSomatic-TO [95% CI] | ΔF1 vs ClairS-TO [95% CI] |", "|---|---|---|---|---|---|---|---|"]
    for tumor in ("HG008T", "COLO829T"):
        for p, pname in PLAT[tumor]:
            for proto, tname in PROTO:
                c = cells[(proto, tumor, p)]
                assert sorted(c["seeds"]) == ["0", "1", "2"], c["seeds"].keys()
                a = [json.loads(Path(c["seeds"][s]["report"]).read_text())["pon"]["at_threshold"] for s in ("0", "1", "2")]
                cols = []
                for k in ("precision", "recall", "f1"):
                    v = [x[k] for x in a]
                    cols.append((v, st.mean(v), st.stdev(v)))
                d = [(c[f"d_{k}"], c[f"d_{k}_ci"]) for k in ("deepsomatic", "clairs_to")]
                tsv.append("\t".join([tumor, pname, tname] + ["/".join(f"{x:.4f}" for x in v) for v, _, _ in cols] +
                                     [f"{y:.4f}" for _, m, s in cols for y in (m, s)] +
                                     [z for dv, ci in d for z in (f"{dv:+.4f}", f"{ci[0]:+.4f},{ci[1]:+.4f}")]))
                md.append(f"| {tumor} | {pname} | {tname} | " +
                          " | ".join(" / ".join(f"{x:.3f}" for x in v) + f" ({m:.3f} ± {s:.3f})" for v, m, s in cols) + " | " +
                          " | ".join(f"{dv:+.3f} [{ci[0]:+.3f}, {ci[1]:+.3f}]" for dv, ci in d) + " |")
    (HERE / "indel_seed_table.tsv").write_text("\n".join(tsv) + "\n")
    (HERE / "indel_seed_table.md").write_text("\n".join(md) + "\n")
    print("\n".join(md))


if __name__ == "__main__":
    main()
