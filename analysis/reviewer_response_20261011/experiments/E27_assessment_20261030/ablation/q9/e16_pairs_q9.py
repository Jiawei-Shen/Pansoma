"""q9 (assessment 2026-10-30): E16 paired 1-Mb / 5-Mb block bootstrap (2,000 resamples, seed 20261015, eval_bed,
PASS after the PoN) of COLO829T Illumina SNV arm pairs, with the HPRC v1.1 (d9) arm at the validation-selected epoch 8
and at the chr1-chosen epochs, using E16's own functions (experiments/E16_bootstrap/e16_bootstrap.py; read only).
Same truth / BED path for every run, so the block resamples are E16's. Writes only next to this script.
usage: python e16_pairs_q9.py
"""
import csv
import json
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
E16 = Path("/scratch/jshen/Github/Pansoma/analysis/reviewer_response_20261011/experiments/E16_bootstrap")
sys.path.insert(0, str(E16))
import callsets as C  # noqa: E402
import e16_bootstrap as E  # noqa: E402

RC = Path("/scratch/jshen/data/pansoma_net_v2_runs/reviewer_cpu_20261013")
RIS = Path("/scratch/jshen/data/pansoma_net_v2_runs/ris_7fc7396_20261012")
RUNS = {
    # seed 0
    "linear_s0": RC / "COLO829_Illumina_SNV_linear_7fc7396",
    "d46_s0": RC / "COLO829_Illumina_SNV_d46_7fc7396",
    "d9_ep8_ddp4_s0": RC / "COLO829_Illumina_SNV_7fc7396_ddp4_b1024",
    "d9_ep8_normbq_s0": RIS / "COLO829_Illumina_SNV_d9_normbq",
    "d9_ep6_normbq_s0": RC / "COLO829_Illumina_SNV_d9_normbq_sweep_e06",
    # seeds 1 / 2 (norm train linear / d46 of the same seed)
    "linear_s1": RIS / "COLO829_Illumina_SNV_linear_7fc7396_s1",
    "linear_s2": RIS / "COLO829_Illumina_SNV_linear_7fc7396_s2",
    "d46_s1": RIS / "COLO829_Illumina_SNV_d46_7fc7396_s1",
    "d46_s2": RIS / "COLO829_Illumina_SNV_d46_7fc7396_s2",
    "d9_ep8_ddp4_s1": RIS / "COLO829_Illumina_SNV_7fc7396_ddp4_b1024_s1",
    "d9_ep8_ddp4_s2": RIS / "COLO829_Illumina_SNV_7fc7396_ddp4_b1024_s2",
    "d9_ep8_normbq_s1": RIS / "COLO829_Illumina_SNV_d9_normbq_s1",
    "d9_ep8_normbq_s2": RIS / "COLO829_Illumina_SNV_d9_normbq_s2",
    "d9_ep5_normbq_s1": RIS / "COLO829_Illumina_SNV_d9_normbq_s1_sweep_e05",
    "d9_ep6_normbq_s2": RIS / "COLO829_Illumina_SNV_d9_normbq_s2_sweep_e06",
}
PAIRS = [("d46_s0", "linear_s0")]
for s in ("s0", "s1", "s2"):
    for d9 in [k for k in RUNS if k.startswith("d9_") and k.endswith(s)]:
        PAIRS += [(d9, f"linear_{s}"), (f"d46_{s}", d9)]


class U:  # the attributes e16_bootstrap.positions() needs
    def __init__(self, name, run_dir):
        self.name, self.run_dir = name, run_dir
        ev = run_dir / "vcf_chr1" / "eval_bed"
        self.rtg = {"bed": ev / "vcfeval" / "pon_calls"}
        assert (self.rtg["bed"] / "done").exists(), run_dir
        self.report = json.loads((ev / "report.json").read_text())
        self.bed = C.resolved(self.report["bed"])
        self.truth = C.resolved(self.report["truth"])


def main():
    units = {k: U(k, v) for k, v in RUNS.items()}
    beds = {(u.truth, u.bed) for u in units.values()}
    assert len(beds) == 1, beds
    fasta = C.manifest("COLO829T_Illumina", "SNV")["fasta"]
    bs = E.BlockSets(fasta, 2000, 20261015)
    boot, point, checks = {}, {}, []
    for k, u in units.items():
        pos = E.positions(u, "bed")
        for width in E.WIDTHS:
            blocks = bs.get("bed", u.bed, width)
            m = E.block_counts(pos, blocks)
            boot[(k, width)] = (blocks.key, blocks.mult @ m)
            point[(k, width)] = m.sum(axis=0)
        at = u.report["pon"]["at_threshold"]
        tot = [int(x) for x in point[(k, E.WIDTHS[0])]]
        checks.append(dict(name=k, run=u.run_dir.name, rtg=tot, report=[at["tp_baseline"], at["fn"], at["tp_call"], at["fp"]],
                           ok=tot == [at["tp_baseline"], at["fn"], at["tp_call"], at["fp"]], f1=at["f1"]))
    rows = []
    for a, b in PAIRS:
        for width in E.WIDTHS:
            assert boot[(a, width)][0] == boot[(b, width)][0]
            row = dict(A=a, B=b, A_run=units[a].run_dir.name, B_run=units[b].run_dir.name, block_mb=width // 1_000_000)
            ba, bb = E.metrics(boot[(a, width)][1]), E.metrics(boot[(b, width)][1])
            pa, pb = E.metrics(point[(a, width)]), E.metrics(point[(b, width)])
            for i, name in enumerate(("P", "R", "F1")):
                d = ba[i] - bb[i]
                lo, hi, se, _ = E.ci(d)
                ok = d[~(d != d)]
                row.update({f"{name}_A": float(pa[i]), f"{name}_B": float(pb[i]), f"d{name}": float(pa[i] - pb[i]),
                            f"d{name}_lo": lo, f"d{name}_hi": hi, f"d{name}_p_gt0": float((ok > 0).mean())})
            rows.append(row)
    with open(HERE / "e16_pairs_q9.tsv", "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=list(rows[0]), delimiter="\t", lineterminator="\n")
        w.writeheader()
        w.writerows(rows)
    md = ["# q9: COLO829T Illumina SNV paired block bootstrap (E16 functions; chr1, BED, PASS after PoN; 2,000 resamples, seed 20261015)", "",
          "Report check (rtg totals TPb / FN / TPc / FP vs report.json): " +
          "; ".join(f"{c['name']} {'ok' if c['ok'] else 'MISMATCH'} F1 {c['f1']:.4f}" for c in checks), "",
          "| A | B | block Mb | F1 A | F1 B | dP [95% CI] | dR [95% CI] | dF1 [95% CI] | P(dF1>0) |", "|---|---|---|---|---|---|---|---|---|"]
    for r in rows:
        f = lambda m: f"{r[m]:+.3f} [{r[m + '_lo']:+.3f}, {r[m + '_hi']:+.3f}]"  # noqa: E731
        md.append(f"| {r['A']} (`{r['A_run']}`) | {r['B']} (`{r['B_run']}`) | {r['block_mb']} | {r['F1_A']:.3f} | {r['F1_B']:.3f} | "
                  f"{f('dP')} | {f('dR')} | {f('dF1')} | {r['dF1_p_gt0']:.3f} |")
    (HERE / "e16_pairs_q9.md").write_text("\n".join(md) + "\n")
    print("\n".join(md))


if __name__ == "__main__":
    main()
