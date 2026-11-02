"""q9 (assessment 2026-10-30): E16 paired 1-Mb / 5-Mb block bootstrap (2,000 resamples, seed 20261015, eval_bed,
PASS after the PoN) of COLO829T Illumina SNV arm pairs, with the HPRC v1.1 (d9) arm at the validation-selected epoch 8
and at the chr1-chosen epochs, using E16's own functions (experiments/E16_bootstrap/e16_bootstrap.py; read only).
Same truth / BED path for every run, so the block resamples are E16's. Writes only next to this script.
usage: python e16_pairs_q4.py
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
SETS = {"HG008T_Illumina": "HG008_Illumina", "HG008T_PacBio": "HG008_PacBio", "HG008T_ONT": "HG008_ONT",
        "COLO829T_Illumina": "COLO829_Illumina", "COLO829T_fiberseq": "COLO829_fiberseq", "COLO829T_ONT": "COLO829_ONT"}


def run_dir(name):
    for base in (RIS, RC):
        d = base / name
        if (d / "vcf_chr1" / "eval_bed" / "vcfeval" / "pon_calls" / "done").exists():
            return d
    raise SystemExit(f"no rtg output for {name}")


def runs_of(p):
    out = {}
    for s, sfx in (("s0", ""), ("s1", "_s1"), ("s2", "_s2")):
        out[f"linear_{s}"] = run_dir(f"{p}_SNV_linear_7fc7396{sfx}")
        out[f"d9_{s}"] = run_dir(f"{p}_SNV_7fc7396_ddp4_b1024{sfx}")
        out[f"d46_{s}"] = run_dir(f"{p}_SNV_d46_7fc7396{sfx}")
    return out


PAIRS = [(f"{a}_{s}", f"{b}_{s}") for s in ("s0", "s1", "s2") for a, b in (("d9", "linear"), ("d46", "linear"), ("d46", "d9"))]


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
    allrows, allchecks = [], []
    for setname, prefix in SETS.items():
        rows, checks = one(setname, prefix)
        allrows += rows
        allchecks += checks
    write(allrows, allchecks)


def one(setname, prefix):
    RUNS = runs_of(prefix)
    units = {k: U(k, v) for k, v in RUNS.items()}
    beds = {(u.truth, u.bed) for u in units.values()}
    assert len(beds) == 1, beds
    fasta = C.manifest(setname, "SNV")["fasta"]
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
            row = dict(set=setname, A=a, B=b, A_run=units[a].run_dir.name, B_run=units[b].run_dir.name, block_mb=width // 1_000_000)
            ba, bb = E.metrics(boot[(a, width)][1]), E.metrics(boot[(b, width)][1])
            pa, pb = E.metrics(point[(a, width)]), E.metrics(point[(b, width)])
            for i, name in enumerate(("P", "R", "F1")):
                d = ba[i] - bb[i]
                lo, hi, se, _ = E.ci(d)
                ok = d[~(d != d)]
                row.update({f"{name}_A": float(pa[i]), f"{name}_B": float(pb[i]), f"d{name}": float(pa[i] - pb[i]),
                            f"d{name}_lo": lo, f"d{name}_hi": hi, f"d{name}_p_gt0": float((ok > 0).mean())})
            rows.append(row)
    for c in checks:
        c["set"] = setname
    return rows, checks


def write(rows, checks):
    with open(HERE / "e16_pairs_q4.tsv", "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=list(rows[0]), delimiter="\t", lineterminator="\n")
        w.writeheader()
        w.writerows(rows)
    md = ["# q9: SNV paired block bootstrap, all six sets x 3 seeds, --norm train, validation-selected checkpoints (E16 functions; chr1, BED, PASS after PoN; 2,000 resamples, seed 20261015)", "",
          "Report check (rtg totals TPb / FN / TPc / FP vs report.json): " +
          "; ".join(f"{c['set']} {c['name']} {'ok' if c['ok'] else 'MISMATCH'} F1 {c['f1']:.4f}" for c in checks), "",
          "| set | A | B | block Mb | F1 A | F1 B | dP [95% CI] | dR [95% CI] | dF1 [95% CI] | P(dF1>0) |", "|---|---|---|---|---|---|---|---|---|---|"]
    for r in rows:
        f = lambda m: f"{r[m]:+.3f} [{r[m + '_lo']:+.3f}, {r[m + '_hi']:+.3f}]"  # noqa: E731
        md.append(f"| {r['set']} | {r['A']} (`{r['A_run']}`) | {r['B']} (`{r['B_run']}`) | {r['block_mb']} | {r['F1_A']:.3f} | {r['F1_B']:.3f} | "
                  f"{f('dP')} | {f('dR')} | {f('dF1')} | {r['dF1_p_gt0']:.3f} |")
    (HERE / "e16_pairs_q4.md").write_text("\n".join(md) + "\n")
    print("\n".join(md))


if __name__ == "__main__":
    main()
