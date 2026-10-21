"""SMaHT strata (chr1, benchmark BED, PASS after the PoN) of the COLO829T Illumina HPRC v1.1 (d9) runs at epochs
6 / 5 / 6, the three seeds R3.1 / R3.9 use for COLO829T Illumina v1.1 (--norm dataset BQ-only, RIS 7fc7396 sweeps).

region_strata.py has COLO829T Illumina d9 only at the validation-selected epoch 8 (COLO829_Illumina_SNV_7fc7396_ddp4_b1024),
which calls unlabelled tensors in chr1:121-126 / 143-144 Mb (R3.9). Same stratum assignment and tally as region_strata.py.
Writes colo829t_illumina_d9_epochs.tsv next to this script.
usage: /wanglab/jshen/anaconda3/bin/python colo829t_illumina_d9_epochs.py
"""
import csv
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parent / "call_level_common"))
sys.path.insert(0, str(HERE))
import calls_lib as L  # noqa: E402
import region_strata as RS  # noqa: E402

RIS = Path("/scratch/jshen/data/pansoma_net_v2_runs/ris_7fc7396_20261012")
RUNS = ("COLO829_Illumina_SNV_d9_normbq_sweep_e06", "COLO829_Illumina_SNV_d9_normbq_s1_sweep_e05",
        "COLO829_Illumina_SNV_d9_normbq_s2_sweep_e06")


def main():
    st = RS.Strata()
    seq = L.chrom_seq(L.manifest("COLO829T_Illumina", "SNV")["fasta"])
    out = []
    for run in RUNS:
        rows = []
        for part in RS.PARTS:
            for rec in L.read_records(RIS / run / "vcf_chr1" / "eval_bed" / "vcfeval" / "pon_calls" / f"{part}.vcf.gz"):
                key = next((k for k in (L.normalize(seq, rec.pos, rec.ref, a) for a in rec.alts) if k), None)
                if key is not None:
                    rows.append(dict(eval="eval_bed", part=part, smaht=st.smaht_of(L.interval(*key)[0])))
        d = RS.tally(rows, "eval_bed", "smaht", RS.SMAHT)
        for s in RS.SMAHT:
            x = d[s]
            out.append([run, s, x["truth"], x["calls"], x["fp"], f"{x['recall']:.3f}", f"{x['precision']:.3f}", f"{x['f1']:.3f}"])
            print("\t".join(map(str, out[-1])))
    with open(HERE / "colo829t_illumina_d9_epochs.tsv", "w", newline="") as f:
        w = csv.writer(f, delimiter="\t", lineterminator="\n")
        w.writerow(["run", "stratum", "truth", "calls", "fp", "recall", "precision", "f1"])
        w.writerows(out)


if __name__ == "__main__":
    main()
