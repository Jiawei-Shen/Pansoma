"""E10b / E16 SNV summary with the COLO829T Illumina HPRC v1.1 (d9) arm at epoch 6 (the user's decision of 2026-10-20)
in place of the epoch-8 run of strata_paired.tsv (GIAB) and ../E16_bootstrap/strata_paired.tsv (SMaHT).

Reads strata_paired.tsv, ../E16_bootstrap/strata_paired.tsv and colo_ill_d9_ep6.json (colo_ill_d9_ep6.py); writes
summary_ep6.md: per stratum, the 12 SNV pairs (HPRC - linear graph): dF1 range, > 0, CI above 0, CI below 0.
usage: /wanglab/jshen/anaconda3/bin/python summary_ep6.py
"""
import csv
import json
from pathlib import Path

HERE = Path(__file__).resolve().parent
LABEL = {n: l for n, _, l in __import__("giab_strata").STRATA}


def rows(path):
    return list(csv.DictReader(open(path), delimiter="\t"))


def main():
    ep6 = {r["stratum"]: r for r in json.loads((HERE / "colo_ill_d9_ep6.json").read_text()) if r["arm"] == "d9_ep6"}
    out = []
    for src, name_of in ((HERE.parent / "E16_bootstrap" / "strata_paired.tsv", lambda s: f"SMaHT {s}"),
                         (HERE / "strata_paired.tsv", lambda s: LABEL[s])):
        groups = {}
        for r in rows(src):
            if r.get("kind", "SNV") != "SNV":
                continue
            lab = name_of(r["stratum"])
            d, lo, hi = float(r["dF1"]), float(r["dF1_lo"]), float(r["dF1_hi"])
            if r["set"] == "COLO829T_Illumina" and r["arm"] == "d9":
                e = ep6[lab]
                d, lo, hi = e["dF1"], e["dF1_lo"], e["dF1_hi"]
            if d != d:  # nan
                continue
            groups.setdefault(lab, []).append((d, lo, hi))
        for lab, v in groups.items():
            out.append(dict(stratum=lab, pairs=len(v), dmin=min(x[0] for x in v), dmax=max(x[0] for x in v),
                            pos=sum(x[0] > 0 for x in v), above=sum(x[1] > 0 for x in v), below=sum(x[2] < 0 for x in v)))
    L = ["# SNV, HPRC - GRCh38-only (linear) graph per stratum, COLO829T Illumina v1.1 at epoch 6", "",
         "| stratum | pairs | dF1 range | dF1 > 0 | CI above 0 | CI below 0 |", "|---|---:|---|---:|---:|---:|"]
    for o in out:
        L.append(f"| {o['stratum']} | {o['pairs']} | {o['dmin']:+.3f} to {o['dmax']:+.3f} | {o['pos']} | {o['above']} | {o['below']} |")
    (HERE / "summary_ep6.md").write_text("\n".join(L) + "\n")
    print("\n".join(L))


if __name__ == "__main__":
    main()
