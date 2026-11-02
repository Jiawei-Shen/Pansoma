"""Assessment 2026-10-30: paired 1-Mb block bootstrap (E16 functions, 2,000 resamples, seed 20261015; chr1, benchmark
BED, PASS after the PoN) of the seed-0 `--norm train` call sets of the main three-reference table (validation-selected
checkpoints everywhere): HPRC v1.1 matched (ddp4) - GRCh38-only, HPRC v2.1 - GRCh38-only, HPRC v2.1 - v1.1, six datasets.
Reads E16 code and run outputs only; writes next to this script. usage: python pairs.py"""
import csv
import json
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
E16 = Path("/scratch/jshen/Github/Pansoma/analysis/reviewer_response_20261011/experiments/E16_bootstrap")
sys.path.insert(0, str(E16))
import callsets as C  # noqa: E402
import e16_bootstrap as E  # noqa: E402

BASES = [Path("/scratch/jshen/data/pansoma_net_v2_runs/reviewer_cpu_20261013"),
         Path("/scratch/jshen/data/pansoma_net_v2_runs/ris_7fc7396_20261012")]
SETS = {"HG008_Illumina": "HG008T_Illumina", "HG008_PacBio": "HG008T_PacBio", "HG008_ONT": "HG008T_ONT",
        "COLO829_Illumina": "COLO829T_Illumina", "COLO829_fiberseq": "COLO829T_fiberseq", "COLO829_ONT": "COLO829T_ONT"}
ARMS = {"linear": "_SNV_linear_7fc7396", "d9": "_SNV_7fc7396_ddp4_b1024", "d46": "_SNV_d46_7fc7396"}
PAIRS = [("d9", "linear"), ("d46", "linear"), ("d46", "d9")]


def run_dir(name):
    for b in BASES:
        if (b / name / "vcf_chr1/eval_bed/vcfeval/pon_calls/done").exists():
            return b / name
    raise FileNotFoundError(name)


class U:
    def __init__(self, d):
        self.name, self.run_dir = d.name, d
        ev = d / "vcf_chr1" / "eval_bed"
        self.rtg = {"bed": ev / "vcfeval" / "pon_calls"}
        self.report = json.loads((ev / "report.json").read_text())
        self.bed = C.resolved(self.report["bed"])
        self.truth = C.resolved(self.report["truth"])


def main():
    width = E.WIDTHS[0]
    rows = []
    for prefix, tset in SETS.items():
        units = {a: U(run_dir(prefix + s)) for a, s in ARMS.items()}
        assert len({(u.truth, u.bed) for u in units.values()}) == 1, prefix
        bs = E.BlockSets(C.manifest(tset, "SNV")["fasta"], 2000, 20261015)
        boot, point = {}, {}
        for a, u in units.items():
            blocks = bs.get("bed", u.bed, width)
            m = E.block_counts(E.positions(u, "bed"), blocks)
            boot[a], point[a] = (blocks.key, blocks.mult @ m), m.sum(axis=0)
            at = u.report["pon"]["at_threshold"]
            assert [int(x) for x in point[a]] == [at["tp_baseline"], at["fn"], at["tp_call"], at["fp"]], u.name
        for a, b in PAIRS:
            assert boot[a][0] == boot[b][0]
            ba, bb, pa, pb = E.metrics(boot[a][1]), E.metrics(boot[b][1]), E.metrics(point[a]), E.metrics(point[b])
            row = dict(dataset=tset, A=a, B=b, A_run=units[a].name, B_run=units[b].name)
            for i, k in enumerate(("P", "R", "F1")):
                lo, hi, _, _ = E.ci(ba[i] - bb[i])
                row.update({f"{k}_A": round(float(pa[i]), 4), f"{k}_B": round(float(pb[i]), 4),
                            f"d{k}": round(float(pa[i] - pb[i]), 4), f"d{k}_lo": round(lo, 4), f"d{k}_hi": round(hi, 4)})
            rows.append(row)
            print(row, flush=True)
    with open(HERE / "pairs.tsv", "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=list(rows[0]), delimiter="\t", lineterminator="\n")
        w.writeheader()
        w.writerows(rows)


if __name__ == "__main__":
    main()
