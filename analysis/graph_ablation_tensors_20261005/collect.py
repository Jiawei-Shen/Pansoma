"""Graph ablation (GRCh38-only linear vs HPRC v1.1 d9 vs HPRC v2.1 d46): tensor and label counts, truth recall,
discovery and run figures of every finalized set, read from the sets' own manifests (pansoma_v2_tensors/<set>)."""
import json, os, sys
from pathlib import Path
T = Path("/scratch/jshen/data/pansoma_v2_tensors")
GRAPHS = (("linear", "_linear"), ("v1.1 d9", ""), ("v2.1 d46", "_d46"))
SETS = [f"{s}_{p}" for s in ("HG008T", "COLO829T") for p in ("PacBio", "fiberseq", "ONT", "Illumina")
        if not (s == "HG008T" and p == "fiberseq") and not (s == "COLO829T" and p == "PacBio")]
def load(p):
    return json.loads(Path(p).read_text()) if Path(p).exists() else None
rows = []
for base in SETS:
    for graph, suffix in GRAPHS:
        s = base + suffix; root = T / s
        st = load(root / "run/status.json") or {}
        if st.get("status") != "finalized" or not st.get("labeled", True) or not (root / "tensors/truth_recall.json").exists():
            rows.append(dict(set=base, graph=graph, state=st.get("status") or "not built")); continue
        q = load(root / "run/queue_status.json") or {}
        cfg = load(root / "run/config.json")
        disc = load(root / "discovery/discovery_report.json") or {}
        rec = load(root / "tensors/truth_recall.json")
        r = dict(set=base, graph=graph, state="finalized", H=cfg["builder"]["haplotypes"], tasks=cfg["tasks"],
                 nodes_observed=disc.get("nodes_observed"), nodes_selected=disc.get("nodes_selected"),
                 tensors=st.get("tensors"), SNV=st["tensors_by_type"]["SNV"], INDEL=st["tensors_by_type"]["INDEL"],
                 task_hours=round(q.get("elapsed_seconds", 0) / 3600, 1), processes=q.get("processes"),
                 peak_gib=round((st.get("peak_sampled_rss_kib") or 0) / 1048576, 1))
        for kind in ("SNV", "INDEL"):
            tot = load(root / f"tensors/{kind}/labels.manifest.json")["totals"]
            for lab in ("somatic", "germline", "non", "ignore"):
                r[f"{kind}_{lab}"] = tot.get(lab, 0)
        for truth in ("somatic", "germline"):
            t = rec[truth]; n = t["alleles"]; s_ = t["status"]
            r[f"{truth}_alleles"] = n
            r[f"{truth}_rep"] = s_.get("tensor_representative", 0)
            r[f"{truth}_any"] = s_.get("tensor_representative", 0) + s_.get("tensor_non_representative_allele", 0)
            for k in ("SNP", "INS", "DEL"):
                b = t["by_kind"].get(k, {}); kn = sum(v for kk, v in b.items() if ":" not in kk)
                r[f"{truth}_{k}_rep"] = round(100 * b.get("tensor_representative", 0) / kn, 2) if kn else None
        rows.append(r)
out = Path(sys.argv[1]) if len(sys.argv) > 1 else Path(__file__).with_name("summary.json")
out.write_text(json.dumps(rows, indent=1) + "\n")
print(json.dumps(rows))
