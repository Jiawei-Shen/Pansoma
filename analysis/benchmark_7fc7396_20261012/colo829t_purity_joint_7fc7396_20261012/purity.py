"""What each COLO829T purity dilution gives the model on chr1, before any calling. Standard library only:

    python3 purity.py          # or as an LSF job: ris/submit.sh -g 0 -n 2 -m 8 NAME py.sh $PWD/purity.py

Per platform and purity (P100 = the undiluted set COLO829T_<platform>, P75 / P50 / P25 = COLO829T_<platform>_P<x>):

- truth: the chr1 somatic truth alleles that pass and lie in the BED (<set>/tensors/somatic.recall.tsv), SNV (SNP) and
  INDEL (INS + DEL), and how many have a tensor (status tensor_representative or tensor_non_representative_allele):
  the recall ceiling of any model on that set. The other statuses: filtered (the candidate failed the build's AF 0.08 /
  ALT-read floor) and no_candidate.
- AF of the labelled chr1 tensors (<set>/tensors/<KIND>/chr1_variant_summary.ndjson `af`, the representative allele's
  read AF, joined line by line to chr1_labels.ndjson): median and quartiles of the somatic (label 1) and germline
  (label 2) tensors. The somatic median relative to P100's is the measured purity; it is biased upward at low purity,
  because tensors exist only at AF >= 0.08.

Writes purity.md and purity.json next to this script.
"""
import json
import re
import statistics
from pathlib import Path

DATA = Path("/storage2/fs1/epigenome/Active/jiawei.s/Dataset/pansoma_v2_tensors")
HERE = Path(__file__).resolve().parent
PLATFORMS = ("Illumina", "fiberseq", "ONT")
PURITIES = ("P100", "P75", "P50", "P25")
CID = re.compile(r'"candidate_id": "([^"]+)"')
AF = re.compile(r'"af": ([0-9.eE+-]+)')
LABEL = re.compile(r'"label": (-?\d+)')


def set_name(platform, purity):
    return f"COLO829T_{platform}" + ("" if purity == "P100" else f"_{purity}")


def truth(t):
    """chr1 truth alleles (passed, in BED) per kind, and those with a tensor, from somatic.recall.tsv."""
    out = {k: {"truth": 0, "with_tensor": 0, "filtered": 0, "no_candidate": 0} for k in ("SNV", "INDEL")}
    with open(t / "somatic.recall.tsv") as f:
        head = f.readline().rstrip("\n").split("\t")
        ix = {c: i for i, c in enumerate(head)}
        for line in f:
            r = line.rstrip("\n").split("\t")
            if r[ix["chrom"]] != "chr1" or r[ix["passed"]] != "True" or r[ix["in_bed"]] != "True":
                continue
            k = "SNV" if r[ix["kind"]] == "SNP" else "INDEL"
            o, s = out[k], r[ix["status"]]
            o["truth"] += 1
            if s.startswith("tensor_"):
                o["with_tensor"] += 1
            elif s in ("filtered", "no_candidate"):
                o[s] += 1
    for o in out.values():
        o["ceiling"] = round(o["with_tensor"] / o["truth"], 4) if o["truth"] else None
    return out


def quartiles(v):
    if not v:
        return None
    q = statistics.quantiles(v, n=4) if len(v) > 1 else [v[0]] * 3
    return {"n": len(v), "q1": round(q[0], 4), "median": round(statistics.median(v), 4), "q3": round(q[2], 4)}


def afs(t, kind):
    """AF quartiles of the chr1 tensors labelled somatic (1) and germline (2)."""
    som, ger = [], []
    with open(t / kind / "chr1_variant_summary.ndjson") as fs, open(t / kind / "chr1_labels.ndjson") as fl:
        for n, (s, l) in enumerate(zip(fs, fl)):
            cs, cl = CID.search(s).group(1), CID.search(l).group(1)
            if cs != cl:
                raise SystemExit(f"{t}/{kind} chr1 line {n}: summary {cs} vs labels {cl}")
            lab = int(LABEL.search(l).group(1))
            if lab in (1, 2):
                (som if lab == 1 else ger).append(float(AF.search(s).group(1)))
    return {"somatic": quartiles(som), "germline": quartiles(ger)}


def main():
    res = {}
    for p in PLATFORMS:
        for q in PURITIES:
            t = DATA / set_name(p, q) / "tensors"
            if not all((t / k / "labels.manifest.json").is_file() for k in ("SNV", "INDEL")):
                print(f"{set_name(p, q)}: not ready, skipped"); continue
            res[f"{p} {q}"] = {"set": set_name(p, q), "truth": truth(t), "af": {k: afs(t, k) for k in ("SNV", "INDEL")}}
            print(f"{set_name(p, q)} done", flush=True)
    (HERE / "purity.json").write_text(json.dumps(res, indent=1) + "\n")
    L = ["# COLO829T dilutions on chr1, before calling (`purity.py`)\n",
         "Truth: chr1 somatic truth alleles that pass and lie in the BED; with a tensor = recall ceiling of any model. "
         "AF: median (quartiles) of the read AF of the chr1 tensors labelled somatic / germline; measured purity = the "
         "somatic SNV median relative to the platform's P100 (biased upward at low purity: tensors need AF >= 0.08).\n",
         "| platform | purity | SNV truth | SNV with a tensor (ceiling) | SNV filtered / no candidate | INDEL truth | "
         "INDEL with a tensor (ceiling) | somatic SNV tensors: AF median (q1-q3) | measured purity | germline SNV AF median | "
         "somatic INDEL tensors: AF median |",
         "|---|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|"]
    for p in PLATFORMS:
        base = (res.get(f"{p} P100") or {}).get("af", {}).get("SNV", {}).get("somatic") or {}
        for q in PURITIES:
            r = res.get(f"{p} {q}")
            if not r:
                L.append(f"| {p} | {q} | not ready | | | | | | | | |"); continue
            ts, ti = r["truth"]["SNV"], r["truth"]["INDEL"]
            s, g, si = r["af"]["SNV"]["somatic"], r["af"]["SNV"]["germline"], r["af"]["INDEL"]["somatic"]
            mp = f"{s['median'] / base['median']:.2f}" if s and base else "-"
            L.append(f"| {p} | {q} | {ts['truth']} | {ts['with_tensor']} ({ts['ceiling']:.3f}) | {ts['filtered']} / "
                     f"{ts['no_candidate']} | {ti['truth']} | {ti['with_tensor']} ({ti['ceiling']:.3f}) | "
                     f"{s['median']:.3f} ({s['q1']:.3f}-{s['q3']:.3f}), n {s['n']:,} | {mp} | "
                     f"{g['median']:.3f} | {si['median']:.3f} (n {si['n']:,}) |" if s and g and si else
                     f"| {p} | {q} | {ts['truth']} | {ts['with_tensor']} | | {ti['truth']} | {ti['with_tensor']} | - | - | - | - |")
    (HERE / "purity.md").write_text("\n".join(L) + "\n")
    print("\n".join(L))


if __name__ == "__main__":
    main()
