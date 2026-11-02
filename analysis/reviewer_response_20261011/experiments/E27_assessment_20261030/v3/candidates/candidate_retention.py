"""Candidates per reference (letter Supplementary Table R1.2): for HG008T and COLO829T, three platforms, on the
GRCh38-only, HPRC v1.1 (d9) and HPRC v2.1 (d46) graphs, how many matched-normal germline alleles and how many somatic
truth SNVs become a candidate tensor whose representative allele is the truth allele (truth_recall.json
status.tensor_representative), chr1-22, all truth records. Writes candidate_retention.md / .json next to this script."""
import json
from pathlib import Path
D = Path("/scratch/jshen/data/pansoma_v2_tensors")
HERE = Path(__file__).resolve().parent
SETS = [("HG008T", "Illumina", "Illumina"), ("HG008T", "PacBio", "PacBio HiFi"), ("HG008T", "ONT", "ONT"),
        ("COLO829T", "Illumina", "Illumina"), ("COLO829T", "fiberseq", "Fiber-seq"), ("COLO829T", "ONT", "ONT")]
REFS = [("linear", "_linear"), ("d9", ""), ("d46", "_d46")]
out = {}
for t, p, label in SETS:
    for ref, suf in REFS:
        j = json.load(open(D / f"{t}_{p}{suf}" / "tensors" / "truth_recall.json"))
        g, s = j["germline"], j["somatic"]
        gs = g["by_kind"]["SNP"]; ss = s["by_kind"]["SNP"]
        out[f"{t}|{label}|{ref}"] = dict(
            germline_alleles=g["alleles"], germline_cand=g["status"].get("tensor_representative", 0),
            germline_snv_alleles=sum(v for k, v in gs.items() if ":" not in k), germline_snv_cand=gs.get("tensor_representative", 0),
            somatic_snv_alleles=sum(v for k, v in ss.items() if ":" not in k), somatic_snv_cand=ss.get("tensor_representative", 0))
(HERE / "candidate_retention.json").write_text(json.dumps(out, indent=1))
L = ["| tumor, platform | reference | germline alleles that become candidates | change vs GRCh38-only | somatic truth SNVs with an exact candidate |",
     "|---|---|---|---|---|"]
for t, p, label in SETS:
    base = out[f"{t}|{label}|linear"]["germline_cand"]
    for ref, name in (("linear", "GRCh38-only"), ("d9", "HPRC v1.1"), ("d46", "HPRC v2.1")):
        x = out[f"{t}|{label}|{ref}"]
        ch = "–" if ref == "linear" else f"−{100 * (1 - x['germline_cand'] / base):.1f}%"
        L.append(f"| {t} {label} | {name} | {x['germline_cand']:,} of {x['germline_alleles']:,} ({100 * x['germline_cand'] / x['germline_alleles']:.1f}%) | {ch} | "
                 f"{x['somatic_snv_cand']:,} of {x['somatic_snv_alleles']:,} ({100 * x['somatic_snv_cand'] / x['somatic_snv_alleles']:.1f}%) |")
(HERE / "candidate_retention.md").write_text("\n".join(L) + "\n")
print("\n".join(L))
# derived ranges
import itertools
red9 = [1 - out[f"{t}|{l}|d9"]["germline_cand"] / out[f"{t}|{l}|linear"]["germline_cand"] for t, p, l in SETS]
red46 = [1 - out[f"{t}|{l}|d46"]["germline_cand"] / out[f"{t}|{l}|linear"]["germline_cand"] for t, p, l in SETS]
red46_9 = [1 - out[f"{t}|{l}|d46"]["germline_cand"] / out[f"{t}|{l}|d9"]["germline_cand"] for t, p, l in SETS]
print("v1.1 vs linear %.1f-%.1f%%, v2.1 vs linear %.1f-%.1f%%, v2.1 vs v1.1 %.1f-%.1f%%" % (100*min(red9),100*max(red9),100*min(red46),100*max(red46),100*min(red46_9),100*max(red46_9)))
for t, p, l in SETS:
    a = [out[f"{t}|{l}|{r}"]["somatic_snv_cand"] / out[f"{t}|{l}|{r}"]["somatic_snv_alleles"] for r in ("linear", "d9", "d46")]
    print(t, l, "somatic SNV exact-candidate share linear/d9/d46: " + " / ".join(f"{100*x:.1f}%" for x in a), "d46-d9 %+.1f pp" % (100*(a[2]-a[1])))
