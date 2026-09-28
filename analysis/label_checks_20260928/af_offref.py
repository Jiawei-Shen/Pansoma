"""(a) SNV AF thresholds: labels kept per platform; (b) off-reference (anchored) tensors by label, per platform."""
import json, glob, collections
R = "/scratch/jshen/data/pansoma_v2_tensors"
P = {"PacBio": f"{R}/Liss_lab_PacBio_Revio_20240125/v3_tensors", "ONT": f"{R}/Liss_lab_Northeastern-ONT-UL-20241216/v3_tensors",
     "Illumina": f"{R}/Liss_lab_BCM_Illumina-WGS_20240313/v3_tensors"}
TH = (0.06, 0.07, 0.08, 0.10, 0.15, 0.20, 0.30)
for n, d in P.items():
    af = collections.defaultdict(list); off = collections.Counter()
    for kind in ("SNV", "INDEL"):
        for f in sorted(glob.glob(f"{d}/{kind}/chr*_labels.ndjson")):
            summ = open(f.replace("_labels.ndjson", "_variant_summary.ndjson"))
            for a, b in zip(open(f), summ):
                l = json.loads(a)
                if l["grch38"] is None and l.get("anchor"):
                    off[(kind, l["label"], l["reason"])] += 1
                if kind == "SNV" and l["label"] in (0, 1, 2):
                    i = b.find('"af": '); af[l["label"]].append(float(b[i + 6:b.index(",", i)]))
    print(f"\n=== {n}  SNV tensors kept at AF >= threshold (label 1 somatic / 2 germline / 0 non)")
    print("  threshold " + "  ".join(f"{t:>6}" for t in TH))
    for lab in (1, 2, 0):
        v = af[lab]
        print(f"  label {lab:>2}  " + "  ".join(f"{sum(x >= t for x in v):>6}" if n != "Illumina" or lab != 0 else f"{sum(x >= t for x in v)//1000:>5}k" for t in TH))
    print(f"  off-reference anchored tensors by (kind, label, reason):")
    for k, c in sorted(off.items()): print("    ", k, c)
