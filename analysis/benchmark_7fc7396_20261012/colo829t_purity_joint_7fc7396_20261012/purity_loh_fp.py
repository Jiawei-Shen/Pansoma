"""Why P25 scores above P75 / P50 on every platform: the germline false calls of the joint models, per purity, by read
AF against the AFs that tumor-LOH germline variants take in a dilution. Standard library only (login node):

    PANSOMA_RUNS=/storage2/fs1/epigenome/Active/jiawei.s/pansoma_net_v2_runs python3 purity_loh_fp.py

Models: Illumina the per-read BQ normalized joint runs (ris/tools/readnorm.py; the BQ read-source artefact removed),
fiberseq and ONT the joint runs as they are. A heterozygous germline variant of the normal in a tumor LOH region
appears at AF (1 - p) / 2 where the tumor lost its ALT copy ("lost") and (1 + p) / 2 where it kept it ("kept"), p = the
tumor read fraction (0.75 / 0.50 / 0.25). For each chr1 test tensor labelled germline (2) and called somatic (p_somatic
>= the threshold: SNV the PASS threshold, INDEL the chr1 curve's best point), its AF (the representative allele's
read AF), whether the PoN tags its GRCh38 record (<run>/vcf_chr1_<SET>/eval_bed/*.linear.pon.vcf.gz FILTER
PanelOfNormals; no GRCh38 record: "unplaced"), and its class: within 0.06 of "lost", of "kept", or "other". Writes
loh_fp.md and loh_fp.json next to this script.
"""
import glob
import gzip
import json
import os
import re
import statistics
from pathlib import Path

RUNS = Path(os.environ.get("PANSOMA_RUNS", "/storage2/fs1/epigenome/Active/jiawei.s/pansoma_net_v2_runs"))
DATA = Path("/storage2/fs1/epigenome/Active/jiawei.s/Dataset/pansoma_v2_tensors")
HERE = Path(__file__).resolve().parent
PURITY = {"P100": 1.0, "P75": 0.75, "P50": 0.5, "P25": 0.25}
MODELS = {  # (platform, kind): run
    ("Illumina", "SNV"): "COLO829_Illumina_SNV_joint_readnorm_gpu8",
    ("fiberseq", "SNV"): "COLO829_fiberseq_SNV_joint_7fc7396_gpu4",
    ("ONT", "SNV"): "COLO829_ONT_SNV_joint_7fc7396_gpu4",
    ("Illumina", "INDEL"): "COLO829_Illumina_INDEL_joint_readnorm",
    ("fiberseq", "INDEL"): "COLO829_fiberseq_INDEL_joint_7fc7396",
    ("ONT", "INDEL"): "COLO829_ONT_INDEL_joint_7fc7396",
}
WINDOW = 0.06
AF = re.compile(r'"af": ([0-9.eE+-]+)')
CID = re.compile(r'"candidate_id": "([^"]+)"')
BINS = [0, 0.1, 0.2, 0.3, 0.4, 0.6, 0.8, 1.01]


def set_name(platform, q):
    return f"COLO829T_{platform}" + ("" if q == "P100" else f"_{q}")


def bin_of(a):
    return next(f"{lo:.1f}-{hi:.1f}" for lo, hi in zip(BINS, BINS[1:]) if a < hi)


def pon_status(vcf):
    """{candidate id: True if the PoN tags its GRCh38 record}."""
    out = {}
    with gzip.open(vcf, "rt") as f:
        for line in f:
            if line.startswith("#"):
                continue
            c = line.split("\t", 8)
            tagged = "PanelOfNormals" in c[6]
            for cid in re.split("[;,]", c[2]):
                out[cid] = tagged
    return out


def one(platform, kind, run, q):
    s = set_name(platform, q); p = PURITY[q]
    d = RUNS / run
    rep = json.load(open(d / f"vcf_chr1_{s}" / "eval_bed" / "report.json"))
    pon = rep["pon"]
    thr = rep["threshold"] if kind == "SNV" else pon["best"]["p_somatic"]
    af = {}
    with open(DATA / s / "tensors" / kind / "chr1_variant_summary.ndjson") as f:
        for line in f:
            af[CID.search(line).group(1)] = float(AF.search(line).group(1))
    status = pon_status(glob.glob(str(d / f"vcf_chr1_{s}" / "eval_bed" / "*.linear.pon.vcf.gz"))[0])
    som_af, germ = [], []
    with gzip.open(glob.glob(str(d / f"test_chr1_{s}" / "*.predictions.ndjson.gz"))[0], "rt") as f:
        for line in f:
            r = json.loads(line)
            tl = r.get("test_label")
            if tl == 1:
                som_af.append(af[r["candidate_id"]])
            elif tl == 2 and r["p_somatic"] >= thr:
                a = af[r["candidate_id"]]
                st = status.get(r["candidate_id"])
                cls = "other"
                if p < 1 and abs(a - (1 - p) / 2) <= WINDOW:
                    cls = "lost"
                elif abs(a - (1 + p) / 2) <= WINDOW:
                    cls = "kept"
                germ.append({"af": a, "class": cls, "pon": "unplaced" if st is None else ("tagged" if st else "kept")})
    qs = statistics.quantiles(som_af, n=4) if len(som_af) > 1 else [float("nan")] * 3
    count = lambda sel: sum(1 for g in germ if sel(g))  # noqa: E731
    after = lambda g: g["pon"] == "kept"  # noqa: E731
    res = {"set": s, "threshold": thr, "pass_f1": pon["at_threshold"]["f1"], "pass_p": pon["at_threshold"]["precision"],
           "pass_r": pon["at_threshold"]["recall"], "best_f1": pon["best"]["f1"], "best_p": pon["best"]["precision"],
           "best_r": pon["best"]["recall"], "pon_ceiling": pon["ceiling"],
           "somatic_af": {"n": len(som_af), "q1": qs[0], "median": qs[1], "q3": qs[2]},
           "expected": {"lost": (1 - p) / 2 if p < 1 else None, "kept": (1 + p) / 2},
           "germline_called": {"all": len(germ), "after_pon": count(after), "unplaced": count(lambda g: g["pon"] == "unplaced")},
           "by_class": {c: {"before_pon": count(lambda g, c=c: g["class"] == c),
                            "after_pon": count(lambda g, c=c: g["class"] == c and after(g))} for c in ("lost", "kept", "other")},
           "by_af_bin": {b: {"before_pon": count(lambda g, b=b: bin_of(g["af"]) == b),
                             "after_pon": count(lambda g, b=b: bin_of(g["af"]) == b and after(g))} for b in
                         [f"{lo:.1f}-{hi:.1f}" for lo, hi in zip(BINS, BINS[1:])]}}
    return res


def main():
    res = {}
    for (platform, kind), run in MODELS.items():
        for q in PURITY:
            res[f"{platform} {kind} {q}"] = dict(run=run, **one(platform, kind, run, q))
            print(f"{platform} {kind} {q} done", flush=True)
    (HERE / "loh_fp.json").write_text(json.dumps(res, indent=1) + "\n")
    f3 = lambda x: "-" if x is None or x != x else f"{x:.3f}"  # noqa: E731
    L = ["# Germline false calls of the joint models per purity: tumor-LOH germline at AF (1 - p) / 2 and (1 + p) / 2\n",
         "(`purity_loh_fp.py`; chr1, in the BED. Illumina: the per-read BQ normalized models; fiberseq / ONT: the joint "
         "models as they are. SNV at the PASS threshold, INDEL at the chr1 curve's best point. 'lost' / 'kept': germline "
         f"tensors called somatic within {WINDOW} of (1 - p) / 2 / (1 + p) / 2; after the PoN = not tagged by it and "
         "placed on GRCh38.)\n"]
    for kind in ("SNV", "INDEL"):
        L.append(f"\n## {kind}\n")
        L.append("| platform | purity | F1 after PoN (P / R) | best F1 | somatic AF q1 / median / q3 | expected lost / kept AF | "
                 "germline called: all (after PoN) | lost (after PoN) | kept (after PoN) | other (after PoN) |")
        L.append("|---|---|---|---:|---|---|---:|---:|---:|---:|")
        for key, r in res.items():
            platform, k, q = key.split()
            if k != kind:
                continue
            s, e, g, c = r["somatic_af"], r["expected"], r["germline_called"], r["by_class"]
            fl = f"{r['pass_f1']:.3f} ({r['pass_p']:.3f} / {r['pass_r']:.3f})" if kind == "SNV" else \
                f"{r['best_f1']:.3f} ({r['best_p']:.3f} / {r['best_r']:.3f}) at best"
            L.append(f"| {platform} | {q} | {fl} | {r['best_f1']:.3f} | {f3(s['q1'])} / {f3(s['median'])} / {f3(s['q3'])} | "
                     f"{f3(e['lost'])} / {f3(e['kept'])} | {g['all']:,} ({g['after_pon']:,}) | "
                     f"{c['lost']['before_pon']:,} ({c['lost']['after_pon']:,}) | {c['kept']['before_pon']:,} ({c['kept']['after_pon']:,}) | "
                     f"{c['other']['before_pon']:,} ({c['other']['after_pon']:,}) |")
    (HERE / "loh_fp.md").write_text("\n".join(L) + "\n")
    print("\n".join(L))


if __name__ == "__main__":
    main()
