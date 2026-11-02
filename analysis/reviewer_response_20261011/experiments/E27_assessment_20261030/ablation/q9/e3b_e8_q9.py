"""q9 (assessment 2026-10-30), COLO829T Illumina SNV on HPRC v1.1 (d9), epoch 8 vs chr1-chosen epochs:
1. E3b complete-output F1 = 2TP / (2TP + FP + FN + U) (E3b collect.py complete_f1, imported), U = graph-coordinate PASS;
   rescue ON (<run>/vcf_chr1) and OFF (<run>__norescue/vcf_chr1, seed 0 only: no __norescue exists for seeds 1 / 2).
2. E8 latest_models.py VAF split (imported vaf_split): recall of truth with Illumina VAF < 0.08 / >= 0.08 from
   rtg tp-baseline / fn, for the epoch-8 runs E8 does not list.
Reads report.json / rtg VCFs only; writes e3b_e8_q9.md next to this script.
"""
import importlib.util
from pathlib import Path

HERE = Path(__file__).resolve().parent
EXP = Path("/scratch/jshen/Github/Pansoma/analysis/reviewer_response_20261011/experiments")


def load(name, path):
    spec = importlib.util.spec_from_file_location(name, path)
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


E3B = load("e3b_collect", EXP / "E3b_rescue_off_and_indel_threshold" / "collect.py")
E8 = load("e8_latest_models", EXP / "E8_comparators" / "latest_models.py")

RC = Path("/scratch/jshen/data/pansoma_net_v2_runs/reviewer_cpu_20261013")
RIS = Path("/scratch/jshen/data/pansoma_net_v2_runs/ris_7fc7396_20261012")
E3B_RUNS = [("epoch 8, norm train, s0", RC, "COLO829_Illumina_SNV_7fc7396_ddp4_b1024"),
            ("epoch 6, normbq, s0", RC, "COLO829_Illumina_SNV_d9_normbq_sweep_e06"),
            ("d46 s0 (reference)", RC, "COLO829_Illumina_SNV_d46_7fc7396"),
            ("epoch 8, norm train, s1", RIS, "COLO829_Illumina_SNV_7fc7396_ddp4_b1024_s1"),
            ("epoch 8, norm train, s2", RIS, "COLO829_Illumina_SNV_7fc7396_ddp4_b1024_s2"),
            ("epoch 8, normbq, s0", RIS, "COLO829_Illumina_SNV_d9_normbq"),
            ("epoch 8, normbq, s1", RIS, "COLO829_Illumina_SNV_d9_normbq_s1"),
            ("epoch 8, normbq, s2", RIS, "COLO829_Illumina_SNV_d9_normbq_s2"),
            ("epoch 5, normbq, s1", RIS, "COLO829_Illumina_SNV_d9_normbq_s1_sweep_e05"),
            ("epoch 6, normbq, s2", RIS, "COLO829_Illumina_SNV_d9_normbq_s2_sweep_e06")]
VAF_RUNS = [r for r in E3B_RUNS if "d46" not in r[0]]


def main():
    L = ["# q9: COLO829T Illumina SNV d9, E3b complete-output F1 and E8 VAF split (chr1, eval_bed, PASS after PoN)", "",
         "## E3b complete-output F1 (every graph-coordinate PASS call counted as FP)", "",
         "| label | run | placed F1 | ON: graph-coord PASS / complete F1 | OFF: graph-coord PASS / complete F1 |", "|---|---|---|---|---|"]
    for lab, root, run in E3B_RUNS:
        on = E3B.rep(root / run / "vcf_chr1", "eval_bed")
        offd = root / f"{run}__norescue" / "vcf_chr1"
        off = E3B.rep(offd, "eval_bed") if offd.exists() else None
        a = E3B.complete_f1(on)
        b = E3B.complete_f1(off) if off else None
        L.append(f"| {lab} | `{run}` | {on['pon']['at_threshold']['f1']:.4f} | {a[3]:,} / {a[4]:.4f} | "
                 + (f"{b[3]:,} / {b[4]:.4f} |" if b else "- (no __norescue run) |"))
    L += ["", "## E8 VAF split (truth INFO VAF_Ill; recall below / at or above 0.08)", "",
          "| label | run | recall VAF < 0.08 (n) | recall VAF >= 0.08 (n) |", "|---|---|---|---|"]
    for lab, root, run in VAF_RUNS:
        lo, nlo, hi, nhi = E8.vaf_split(root / run / "vcf_chr1" / "eval_bed" / "vcfeval" / "pon_calls")
        L.append(f"| {lab} | `{run}` | {lo:.3f} ({nlo}) | {hi:.3f} ({nhi}) |")
    (HERE / "e3b_e8_q9.md").write_text("\n".join(L) + "\n")
    print("\n".join(L))


if __name__ == "__main__":
    main()
