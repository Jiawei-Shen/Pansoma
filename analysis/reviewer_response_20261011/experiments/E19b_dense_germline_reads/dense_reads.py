"""E19b: what the reads look like at germline-dense sites on GRCh38-only and on the HPRC graphs (R2.5).

Two parts, chr1, inside the benchmark BED, SNV, germline density n50 = the normal's germline records within +-50 bp
(E19's strata 0 / 1 / 2+, read from E19's per-record files):

A. Reads at the truth somatic SNVs (model-free). For every truth SNV with an exact somatic tensor
   (reason representative_allele_is_somatic_truth) in every arm of the dataset (linear = GRCh38-only graph,
   d9 = HPRC v1.1, d46 = HPRC v2.1), per read row of the tensor, outside the site column (50):
     edits = mismatch columns (op X) + insertion / deletion / complex runs (op I, D, complex) in the 101-column window;
     MAPQ  = the row's MAPQ (channel 3).
   Per tensor: mean edits per read, share of reads with >= 2 edits, mean MAPQ; summarised per arm and stratum over the
   sites that have a tensor in all arms (paired).
B. The GRCh38-only arm's false positives (PASS after the PoN, in the BED; E19 / E4b per-record files) and what the
   graph arm does at the same allele: also a false positive / a tensor that is not called / no tensor at that
   position. With E4b's class of each linear FP (germline same allele, near germline, ...).

usage: /wanglab/jshen/anaconda3/bin/python dense_reads.py   (writes tables.md and summary.json next to this script)
"""
import collections
import csv
import gzip
import itertools
import json
import sys
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
CALL = Path("/scratch/jshen/data/reviewer_experiments_20261013/call_level")
TENSORS = Path("/scratch/jshen/data/pansoma_v2_tensors")
SITE_COL = 50
EXACT = "representative_allele_is_somatic_truth"

# dataset -> (tensor set prefix, {arm: call-set name or None})
DATASETS = {
    "HG008T Illumina": ("HG008T_Illumina", {"linear": "HG008_Illumina_SNV_linear_7fc7396__vcf_chr1",
                                            "d9": "HG008_Illumina_SNV_7fc7396_ddp4_b1024__vcf_chr1",
                                            "d46": "HG008_Illumina_SNV_d46_7fc7396__vcf_chr1"}),
    "HG008T PacBio": ("HG008T_PacBio", {"linear": "HG008_PacBio_SNV_linear_7fc7396__vcf_chr1",
                                        "d9": "HG008_PacBio_SNV_7fc7396__vcf_chr1",
                                        "d46": "HG008_PacBio_SNV_d46_7fc7396__vcf_chr1"}),
    "HG008T ONT": ("HG008T_ONT", {"linear": "HG008_ONT_SNV_linear_7fc7396__vcf_chr1",
                                  "d9": "HG008_ONT_SNV_7fc7396__vcf_chr1",
                                  "d46": "HG008_ONT_SNV_d46_7fc7396__vcf_chr1"}),
    "COLO829T Illumina": ("COLO829T_Illumina", {"linear": "COLO829_Illumina_SNV_linear_7fc7396__vcf_chr1", "d9": None,
                                                "d46": "COLO829_Illumina_SNV_d46_7fc7396__vcf_chr1"}),
    "COLO829T fiberseq": ("COLO829T_fiberseq", {"linear": "COLO829_fiberseq_SNV_linear_7fc7396__vcf_chr1", "d9": None,
                                                "d46": "COLO829_fiberseq_SNV_d46_7fc7396__vcf_chr1"}),
    "COLO829T ONT": ("COLO829T_ONT", {"linear": "COLO829_ONT_SNV_linear_7fc7396__vcf_chr1", "d9": None,
                                      "d46": "COLO829_ONT_SNV_d46_7fc7396__vcf_chr1"}),
}
ARMS = ("linear", "d9", "d46")
SUFFIX = {"linear": "_linear", "d9": "", "d46": "_d46"}


def stratum(n):
    n = int(n)
    return "2+" if n >= 2 else str(n)


def read_tsv(path):
    with gzip.open(path, "rt") as f:
        return list(csv.DictReader(f, delimiter="\t"))


def labels(set_name):
    """chr1 SNV labels: exact somatic tensors by (pos1, alt), every tensor's (pos1, alt) -> label, positions."""
    exact, by_allele, by_pos = {}, {}, collections.Counter()
    with open(TENSORS / set_name / "tensors" / "SNV" / "chr1_labels.ndjson") as f:
        for line in f:
            d = json.loads(line)
            g = d.get("grch38")
            if not g:
                continue
            key = (g["pos0"] + 1, g["alt"])
            by_allele[key] = (d["label"], d["reason"])
            by_pos[g["pos0"] + 1] += 1
            if d["label"] == 1 and d["reason"] == EXACT:
                exact[key] = (d["shard_file"], d["index_within_shard"])
    return exact, by_allele, by_pos


def row_metrics(t):
    """Per covered row of one tensor (8, 200, 101): edits outside the site column, and MAPQ."""
    ch0, ch3, ch4 = t[0], t[3], t[4]
    covered = (ch0 != 0).any(axis=1)
    keep = np.ones(t.shape[2], bool)
    keep[SITE_COL] = False
    op = np.where(keep[None, :], ch4, 0)
    x = (op == 2).sum(axis=1)
    runs = 0
    for code in (3, 4, 5):
        on = op == code
        start = on & ~np.concatenate([np.zeros((on.shape[0], 1), bool), on[:, :-1]], axis=1)
        runs = runs + start.sum(axis=1)
    edits = (x + runs)[covered]
    mapq = ch3.max(axis=1)[covered]
    return edits, mapq


def part_a(name, prefix, arms):
    sets = {a: prefix + SUFFIX[a] for a in arms}
    lab = {a: labels(s) for a, s in sets.items()}
    truth = [r for r in read_tsv(CALL / "E19" / f"{DATASETS[name][1]['linear']}.v2.tsv.gz")
             if r["part"] in ("tp-baseline", "fn")]
    per = collections.defaultdict(list)  # (arm, stratum) -> per-site tuples
    n_truth, n_paired = collections.Counter(), collections.Counter()
    wanted = []  # (arm, shard, index, stratum)
    for r in truth:
        key = (int(r["pos"]), r["alt"])
        s = stratum(r["n50"])
        n_truth[s] += 1
        locs = {a: lab[a][0].get(key) for a in arms}
        if any(v is None for v in locs.values()):
            continue
        n_paired[s] += 1
        wanted += [(a, loc[0], loc[1], s) for a, loc in locs.items()]
    # one read per tensor at its byte offset (np.load mmap faults page by page over BeeGFS), in file order
    for (a, shard), group in itertools.groupby(sorted(wanted), key=lambda w: (w[0], w[1])):
        with open(TENSORS / sets[a] / "tensors" / "SNV" / shard, "rb") as f:
            np.lib.format.read_magic(f)
            shape, fortran, dtype = np.lib.format.read_array_header_1_0(f)
            assert not fortran and dtype == np.int8, (shard, dtype)
            offset, size = f.tell(), int(np.prod(shape[1:]))
            for _, _, idx, s in group:
                f.seek(offset + idx * size)
                t = np.frombuffer(f.read(size), dtype=np.int8).reshape(shape[1:])
                e, q = row_metrics(t)
                if len(e) == 0:
                    continue
                per[(a, s)].append((e.mean(), (e >= 2).mean(), q.mean(), (q < 20).mean(), len(e)))
    out = {"truth": dict(n_truth), "paired": dict(n_paired), "arms": {}}
    for a in arms:
        out["arms"][a] = {}
        for s in ("0", "1", "2+"):
            v = np.array(per[(a, s)]) if per[(a, s)] else np.zeros((0, 5))
            out["arms"][a][s] = {
                "n": int(len(v)),
                "edits_mean": float(v[:, 0].mean()) if len(v) else None,
                "edits_median": float(np.median(v[:, 0])) if len(v) else None,
                "reads_ge2_edits": float(v[:, 1].mean()) if len(v) else None,
                "mapq_mean": float(v[:, 2].mean()) if len(v) else None,
                "mapq_lt20": float(v[:, 3].mean()) if len(v) else None,
                "reads": float(v[:, 4].mean()) if len(v) else None}
    return out, lab


def part_b(name, arms, lab):
    calls = DATASETS[name][1]
    lin = {(int(r["pos"]), r["alt"]): stratum(r["n50"])
           for r in read_tsv(CALL / "E19" / f"{calls['linear']}.v2.tsv.gz") if r["part"] == "fp"}
    cls = {(int(r["pos"]), r["alt"]): r["class"]
           for r in read_tsv(CALL / "E4b" / f"{calls['linear']}.eval_bed.fp.v2.tsv.gz")}
    out = {}
    for a in arms:
        if a == "linear" or not calls.get(a):
            continue
        gfp = {(int(r["pos"]), r["alt"]) for r in read_tsv(CALL / "E19" / f"{calls[a]}.v2.tsv.gz") if r["part"] == "fp"}
        exact_, by_allele, by_pos = lab[a]
        res = collections.defaultdict(collections.Counter)
        for key, s in lin.items():
            if key in gfp:
                fate = "also_fp"
            elif key in by_allele:
                fate = "tensor_not_called"
            elif key[0] in by_pos:
                fate = "other_allele_tensor"
            else:
                fate = "no_tensor"
            res[s][fate] += 1
            res[s]["class:" + cls.get(key, "?")] += 1
            res[s]["n"] += 1
            res[s][f"{fate}|{cls.get(key, '?')}"] += 1
        out[a] = {s: dict(c) for s, c in res.items()}
    return out


def fmt(x, d=2):
    return "-" if x is None else f"{x:.{d}f}"


def main():
    summary = {}
    for name, (prefix, calls) in DATASETS.items():
        print(name, file=sys.stderr, flush=True)
        a, lab = part_a(name, prefix, ARMS)
        b = part_b(name, ARMS, lab)
        summary[name] = {"A": a, "B": b}
    (HERE / "summary.json").write_text(json.dumps(summary, indent=1) + "\n")
    L = ["# E19b tables (generated by dense_reads.py)", "",
         "## A. Reads at chr1 truth somatic SNVs (paired: exact tensor in every arm), by germline density n50", "",
         "edits = mismatches + indel runs per read outside the site column (101-column window); "
         ">=2 = share of reads with at least two edits; MAPQ<20 = share of reads.", "",
         "| dataset | stratum | sites (truth / paired) | edits per read linear / d9 / d46 | reads with >=2 edits linear / d9 / d46 "
         "| mean MAPQ linear / d9 / d46 | MAPQ<20 linear / d9 / d46 |",
         "|---|---|---|---|---|---|---|"]
    for name, v in summary.items():
        A = v["A"]
        for s in ("0", "1", "2+"):
            g = lambda k, d=2: " / ".join(fmt(A["arms"][a][s][k], d) for a in ARMS)
            L.append(f"| {name} | {s} | {A['truth'].get(s, 0)} / {A['paired'].get(s, 0)} | {g('edits_mean')} | "
                     f"{g('reads_ge2_edits')} | {g('mapq_mean', 1)} | {g('mapq_lt20')} |")
    L += ["", "## B. GRCh38-only false positives (PASS, after the PoN, in the BED) and the graph arm at the same allele", "",
          "| dataset | graph arm | stratum | linear FP | also FP on graph | tensor, not called | other allele's tensor at the position "
          "| no tensor | linear FP class: germline same / overlap / near / HPRC / other |", "|---|---|---|---|---|---|---|---|---|"]
    for name, v in summary.items():
        for a, st in v["B"].items():
            for s in ("0", "1", "2+"):
                c = st.get(s, {})
                n = c.get("n", 0)
                L.append(f"| {name} | {a} | {s} | {n} | {c.get('also_fp', 0)} | {c.get('tensor_not_called', 0)} | "
                         f"{c.get('other_allele_tensor', 0)} | {c.get('no_tensor', 0)} | "
                         f"{c.get('class:germline_same', 0)} / {c.get('class:germline_overlap', 0)} / "
                         f"{c.get('class:near_germline', 0)} / {c.get('class:hprc', 0)} / {c.get('class:other', 0)} |")
    (HERE / "tables.md").write_text("\n".join(L) + "\n")
    print("\n".join(L))


if __name__ == "__main__":
    main()
