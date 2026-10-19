"""Where the false calls of the COLO829T Illumina SNV models come from (chr1, in the BED, after the PoN, PASS calls:
rtg's <run>/vcf_chr1/eval_bed/vcfeval/pon_calls/fp.vcf.gz). Standard library only (login node):

    PANSOMA_RUNS=/storage2/fs1/epigenome/Active/jiawei.s/pansoma_net_v2_runs python3 d9_fp_check.py

1. Each false call's tensor(s) (the VCF ID: candidate ids) looked up in the run's chr1 predictions: label and reason
   (outside_confident_region = label -1, never trained on), and how many lie in chr1:121-126 or 143-144 Mb.
2. The false calls of the 4-GPU d9 runs in those windows against the d46 set: is there a d46 tensor with the same
   GRCh38 allele (chr1_labels.ndjson grch38 pos0 / ref / alt), and does each d46 run call it (p_somatic >= its PASS
   threshold)?
Writes d9_fp_check.md next to this script.
"""
import gzip
import json
import os
import re
from collections import Counter
from pathlib import Path

RUNS = Path(os.environ.get("PANSOMA_RUNS", "/storage2/fs1/epigenome/Active/jiawei.s/pansoma_net_v2_runs"))
DATA = Path("/storage2/fs1/epigenome/Active/jiawei.s/Dataset/pansoma_v2_tensors")
HERE = Path(__file__).resolve().parent
D9 = ["COLO829_Illumina_SNV_7fc7396_ddp4_b1024", "COLO829_Illumina_SNV_7fc7396_ddp4_b1024_s1",
      "COLO829_Illumina_SNV_7fc7396_ddp4_b1024_s2"]
H90 = ["COLO829_Illumina_SNV_h90", "COLO829_Illumina_SNV_h90_s1", "COLO829_Illumina_SNV_h90_s2"]
D46 = ["COLO829_Illumina_SNV_d46_7fc7396", "COLO829_Illumina_SNV_d46_7fc7396_s1", "COLO829_Illumina_SNV_d46_7fc7396_s2"]
WINDOWS = [(121_000_000, 126_000_000), (143_000_000, 144_000_000)]


def in_window(pos):
    return any(a <= pos < b for a, b in WINDOWS)


def false_calls(run):
    """[(pos, ref, alt, [candidate ids])] of the run's PASS false calls in the BED after the PoN."""
    f = RUNS / run / "vcf_chr1" / "eval_bed" / "vcfeval" / "pon_calls" / "fp.vcf.gz"
    if not f.exists():
        return None
    out = []
    for line in gzip.open(f, "rt"):
        if line[0] == "#":
            continue
        c = line.split("\t", 6)
        out.append((int(c[1]), c[3], c[4], re.split("[;,]", c[2])))
    return out


def predictions(run):
    f = next((RUNS / run / "test_chr1").glob("*.predictions.ndjson.gz"))
    return {d["candidate_id"]: d for d in map(json.loads, gzip.open(f, "rt"))}


def threshold(run):
    return json.load(open(RUNS / run / "vcf_chr1" / "eval_bed" / "report.json"))["threshold"]


def main():
    L = ["# False calls of the COLO829T Illumina SNV models (chr1, in the BED, after the PoN, PASS)\n",
         "(`d9_fp_check.py`; windows chr1:121-126 and 143-144 Mb; reason = the label reason of the call's tensor.)\n",
         "\n## 1. Reasons of the false calls\n",
         "| run | recipe | false calls | in the windows | outside_confident_region (label -1) | of them in the windows | other reasons |",
         "|---|---|---:|---:|---:|---:|---|"]
    d9_window = {}
    for run, recipe in [(r, "d9, 4 GPUs") for r in D9] + [(r, "d9, 1 GPU") for r in H90] + [(r, "d46, 4 GPUs") for r in D46]:
        fps = false_calls(run)
        if fps is None:
            L.append(f"| {run} | {recipe} | not evaluated | | | | |")
            continue
        pred = predictions(run)
        reasons, out_win, n_win = Counter(), 0, 0
        for pos, ref, alt, ids in fps:
            r = next((pred[i]["reason"] for i in ids if i in pred), "no tensor found")
            reasons[r] += 1
            n_win += in_window(pos)
            out_win += in_window(pos) and r == "outside_confident_region"
            if run in D9 and in_window(pos):
                d9_window.setdefault((pos, ref, alt), set()).add(run)
        other = ", ".join(f"{k} {v}" for k, v in reasons.most_common() if k != "outside_confident_region")
        L.append(f"| {run} | {recipe} | {len(fps):,} | {n_win:,} | {reasons['outside_confident_region']:,} | "
                 f"{out_win:,} | {other} |")

    # 2. the d9 window false calls on d46
    alleles = {}
    pos_re = re.compile(r'"grch38": \{"chrom": "chr1", "pos0": (\d+), "ref": "([ACGTN]+)", "alt": "([ACGTN]+)"')
    cid_re = re.compile(r'"candidate_id": "([^"]+)"')
    for line in open(DATA / "COLO829T_Illumina_d46" / "tensors" / "SNV" / "chr1_labels.ndjson"):
        m = pos_re.search(line)
        if m:
            alleles.setdefault((int(m.group(1)) + 1, m.group(2), m.group(3)), []).append(cid_re.search(line).group(1))
    d46 = {r: (predictions(r), threshold(r)) for r in D46 if (RUNS / r / "vcf_chr1" / "eval_bed" / "report.json").exists()}
    keys = sorted(d9_window)
    found = [k for k in keys if k in alleles]
    L += ["\n## 2. The window false calls of the 4-GPU d9 runs on the d46 set\n",
          f"- distinct false-call alleles of the three 4-GPU d9 seeds in the windows: {len(keys):,} "
          f"(in all three seeds: {sum(len(v) == 3 for v in d9_window.values()):,})",
          f"- with a d46 tensor of the same GRCh38 allele: {len(found):,} ({len(found) / max(len(keys), 1):.0%})"]
    for r, (pred, thr) in d46.items():
        called = sum(any(pred.get(c, {}).get("p_somatic", 0) >= thr for c in alleles[k]) for k in found)
        L.append(f"- called by {r} (p_somatic >= {thr:.3f}): {called:,} of the {len(found):,}")
    (HERE / "d9_fp_check.md").write_text("\n".join(L) + "\n")
    print("\n".join(L))


if __name__ == "__main__":
    main()
