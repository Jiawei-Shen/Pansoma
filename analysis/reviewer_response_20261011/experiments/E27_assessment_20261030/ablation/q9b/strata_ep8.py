"""Q9: E10/E10b/E16 strata SNV summaries with COLO829T Illumina v1.1 at the validation-selected epoch 8 (the d9 rows of
strata_paired.tsv = COLO829_Illumina_SNV_7fc7396_ddp4_b1024) vs at epoch 6 (E10b colo_ill_d9_ep6.json), and E10 extreme calls."""
import csv, json, sys
from pathlib import Path
E = Path("/scratch/jshen/Github/Pansoma/analysis/reviewer_response_20261011/experiments")
sys.path.insert(0, str(E / "E10b_giab_strata"))
LABEL = {n: l for n, _, l in __import__("giab_strata").STRATA}
ep6 = {r["stratum"]: r for r in json.loads((E / "E10b_giab_strata/colo_ill_d9_ep6.json").read_text()) if r["arm"] == "d9_ep6"}
def rows(p):
    return list(csv.DictReader(open(p), delimiter="\t"))
out = {}
for variant in ("ep8", "ep6"):
    groups = {}
    for src, name_of in ((E / "E16_bootstrap/strata_paired.tsv", lambda s: f"SMaHT {s}"), (E / "E10b_giab_strata/strata_paired.tsv", lambda s: LABEL[s])):
        for r in rows(src):
            if r.get("kind", "SNV") != "SNV":
                continue
            lab = name_of(r["stratum"])
            d, lo, hi = float(r["dF1"]), float(r["dF1_lo"]), float(r["dF1_hi"])
            if variant == "ep6" and r["set"] == "COLO829T_Illumina" and r["arm"] == "d9":
                e = ep6[lab]; d, lo, hi = e["dF1"], e["dF1_lo"], e["dF1_hi"]
            if d != d:
                continue
            groups.setdefault(lab, []).append((r["set"], r["arm"], d, lo, hi))
    out[variant] = groups
print("| stratum | pairs | dF1 range (ep8) | CI above 0 (ep8) | CI below 0 (ep8) | dF1 range (ep6) | CI above 0 (ep6) | COLO829T Illumina v1.1 dF1 [CI] ep8 / ep6 |")
print("|---|---|---|---|---|---|---|---|")
for lab in out["ep8"]:
    a, b = out["ep8"][lab], out["ep6"][lab]
    c8 = [x for x in a if x[0] == "COLO829T_Illumina" and x[1] == "d9"]; c6 = [x for x in b if x[0] == "COLO829T_Illumina" and x[1] == "d9"]
    f = lambda x: f"{x[2]:+.3f} [{x[3]:+.3f}, {x[4]:+.3f}]" if x else "-"
    print(f"| {lab} | {len(a)} | {min(x[2] for x in a):+.3f} to {max(x[2] for x in a):+.3f} | {sum(x[3] > 0 for x in a)} | {sum(x[4] < 0 for x in a)} | "
          f"{min(x[2] for x in b):+.3f} to {max(x[2] for x in b):+.3f} | {sum(x[3] > 0 for x in b)} | {f(c8[0] if c8 else None)} / {f(c6[0] if c6 else None)} |")
