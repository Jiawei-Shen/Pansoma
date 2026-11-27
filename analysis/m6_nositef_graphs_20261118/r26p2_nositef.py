#!/usr/bin/env python3
"""r26p2_nositef.md: Reviewer 2 comment 6 point 2 (what happens to the germline SNV candidates the graph does not
absorb) on the six-tumor models of this batch, HPRC v1.1 and v2.1 without the population site features as a model
input (user 2026-11-24). The letter's sentence is computed on the HPRC models with the features
(../reviewer_response_20261011/ver4_20261123/r26p2/r26p2.py -> r26p2.md, cells.json).

Same chain and functions as r26p2.py (imported: load_predictions, file_pon, chain, QUANT): per cell, from the run's
chr1 predictions and call set,
  G    germline candidates = test tensors labelled germline (model independent)
  tag  the population filter's SNV tag of the tensor, from the set's site-feature file (these models' predictions carry
       no tag: predict reads no feature file)
  C    G with p_somatic >= the PASS threshold;  K = C on the GRCh38 path and not tagged
  rtg  FP records of raw_calls / pon_calls whose allele the matched normal's germline truth carries
Checks (asserted): the test set, G, tagged G, test tensors and somatic tensors of every cell equal those of the cell of
the same tumor, platform and graph in r26p2's cells.json (same tensors; only the model differs); PASS count = the
graph VCF's; rtg FP = report.json; sha256 of the germline VCFs = the manifests'.
usage: /home/jiawei.s/miniconda3/envs/convnextv2/bin/python r26p2_nositef.py   (login node, a few minutes; numpy)
"""
import json
import sys
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
V4 = HERE.parent / "reviewer_response_20261011/ver4_20261123/r26p2"
sys.path.insert(0, str(HERE))
sys.path.insert(0, str(V4))
import collect as C  # noqa: E402
import r26p2 as R  # noqa: E402
import strata as S  # noqa: E402  (chrom_seq without pysam, sha256; sets calls_lib.chrom_seq)

L = R.L
DATA = S.DATA
MOUNTS = (("/scratch/jshen/data/", DATA), ("/scratch/qfu/COLO829BL_DSA/", DATA / "COLO829BL_DSA"))  # ris/submit.sh
_GERM = {}


def germline(set_name):
    m = json.loads((S.TENSORS / set_name / "tensors" / "SNV" / "labels.manifest.json").read_text())["truth"]["germline"]
    if m["vcf"] not in _GERM:
        vcf = next(d / m["vcf"][len(a):] for a, d in MOUNTS if m["vcf"].startswith(a))
        assert S.sha256(vcf) == m["vcf_sha256"], vcf
        _GERM[m["vcf"]] = L.germline(str(vcf), S.FASTA)
    return _GERM[m["vcf"]]


def germline_same(fp_vcf, germ, seq):
    n = same = 0
    for rec in L.read_records(fp_vcf, with_format=False):
        keys = [L.normalize(seq, rec.pos, rec.ref, alt) for alt in rec.alts]
        assert any(k is not None for k in keys), (fp_vcf, rec.pos)
        n += 1
        same += any(k is not None and k in germ.keys for k in keys)
    return n, same


def cell(g, tumor, p, ref):
    set_name = C.test_set(tumor, p, g)
    run = C.RUNS / f"multicancer6_{p}_SNV_{g}_nositef"
    v = run / f"vcf_chr1_{set_name}"
    st = json.loads((v / f"{set_name}.tensors.SNV.graph.stats.json").read_text())
    rep = json.loads((v / "eval_bed/report.json").read_text())
    args = json.loads((run / "args.json").read_text())
    thr = st["threshold"]
    assert R.RULE in st["call"] and st["offref_rescored"] == 0 and rep["threshold"] == thr
    assert args["seed"] == 0 and args["call_recall"] == 0.95 and args["best_by"] == "postpon" and args["site_features_pon_only"]
    pred = run / f"test_chr1_{set_name}" / f"{set_name}.tensors.SNV.predictions.ndjson.gz"
    assert R.header_value(v / f"{set_name}.tensors.SNV.graph.vcf.gz", "pansoma_predictions") == str(pred)
    d = R.load_predictions(pred)
    assert not d["has_pon"] and len(d["p"]) == st["predictions"]
    tag, _ = R.file_pon(set_name, d, args["site_features"])
    assert int(((d["p"] >= thr) & ~d["exc"]).sum()) == st["called"]
    t = d["in_test"]
    c = R.chain(d["lab"][t], d["off"][t], tag[t] == 1, d["p"][t], thr)
    assert ref["set"] == set_name and all(c[k] == ref[k] for k in ("test_tensors", "somatic", "G", "G_tag", "G_untag_onref")), set_name
    out = dict(graph=g, tumor=tumor, platform=C.plat(tumor, p), set=set_name, threshold=thr, **c)
    germ, seq = germline(set_name), S.chrom_seq(S.FASTA)
    for stage in ("raw", "pon"):
        a = rep[stage]["at_threshold"]
        n, same = germline_same(v / f"eval_bed/vcfeval/{stage}_calls/fp.vcf.gz", germ, seq)
        assert n == a["fp"], (v, stage)
        out[stage] = dict(fp=a["fp"], fp_germline_same=same, precision=a["precision"], recall=a["recall"], f1=a["f1"])
    R.log(f"{g} {tumor} {p}: G {c['G']:,} C {c['C']:,} K {c['K']:,} | rtg FP after {out['pon']['fp']:,} germline {out['pon']['fp_germline_same']:,}")
    return out


def main():
    ref = {(c["graph"], c["tumor"], c["platform"]): c for c in json.loads((V4 / "cells.json").read_text())["cells"]}
    cells = [cell(g, t, p, ref[g, t, C.plat(t, p)]) for g, _ in C.GRAPHS[::-1] for t in ("HG008T", "COLO829T", "HCC1395T") for p in C.PLATFORMS]
    (HERE / "r26p2_nositef.json").write_text(json.dumps(cells, indent=1) + "\n")
    pct = R.pct
    md = ["# Reviewer 2 comment 6 point 2 on the HPRC six-tumor models without the population site features (chr1)", "",
          "Generated by `r26p2_nositef.py`; definitions as `ver4_20261123/r26p2/r26p2.md` (G = germline SNV candidates of "
          "chr1 in the benchmark regions, the germline variants the graph did not absorb; tagged = by the population "
          "filter's SNV rule; C = G called somatic at the PASS threshold; K = C left after the filter; rtg FP = false "
          "calls after the filter). G and the tagged share are those of r26p2.md (same tensors). \"with\" columns = the "
          "HPRC six-tumor models with the site features (r26p2's cells.json, the letter's numbers).", ""]
    for g, gl in C.GRAPHS[::-1]:
        md += [f"## {gl}", "",
               "| tumor | platform | G | tagged | C, with → without | rejected by the model, with → without | K, with → without | "
               "removed by the filter (1 − K/C), with → without | remaining K/G, with → without | rtg FP after the filter, with → without | "
               "germline allele among them, with → without |", "|" + "---|" * 11]
        for c in cells:
            if c["graph"] != g:
                continue
            o = ref[g, c["tumor"], c["platform"]]
            md.append(f"| {c['tumor']} | {c['platform']} | {c['G']:,} | {pct(c['G_tag'], c['G'])} | {o['C']:,} → {c['C']:,} | "
                      f"{pct(o['G'] - o['C'], o['G'])} → {pct(c['G'] - c['C'], c['G'])} | {o['K']:,} → {c['K']:,} | "
                      f"{pct(o['C'] - o['K'], o['C'])} → {pct(c['C'] - c['K'], c['C'])} | {pct(o['K'], o['G'])} → {pct(c['K'], c['G'])} | "
                      f"{o['pon']['fp']:,} → {c['pon']['fp']:,} | {pct(o['pon']['fp_germline_same'], o['pon']['fp'])} → "
                      f"{pct(c['pon']['fp_germline_same'], c['pon']['fp'])} |")
        md.append("")
    md += ["## Ranges for the sentence (HG008T chr1, both graphs; and the nine datasets, both graphs)", "",
           "| quantity | HG008T, with site features (the letter) | HG008T, without | nine datasets, with | nine datasets, without |",
           "|---|---|---|---|---|"]
    for q, label in (("rejects", "the model rejects (1 − C/G)"), ("removes", "the filter removes of those called somatic (1 − K/C)"),
                     ("remain", "remain after classification and filtering (K/G)"),
                     ("share_rtg", "germline alleles among the false calls after filtering (rtg FP)"),
                     ("tag_called", "tagged germline candidates called somatic"),
                     ("untag_called", "untagged germline candidates called somatic")):
        def r(cs, k=2 if q == "tag_called" else 1):
            return R.rng([R.frac(c, q) for c in cs], k)
        hg = [c for c in cells if c["tumor"] == "HG008T"]
        md.append(f"| {label} | {r([ref[c['graph'], c['tumor'], c['platform']] for c in hg])} | {r(hg)} | "
                  f"{r([ref[c['graph'], c['tumor'], c['platform']] for c in cells])} | {r(cells)} |")
    hg = [c for c in cells if c["tumor"] == "HG008T"]
    md += ["", f"Filter removals among the germline candidates called somatic, HG008T (C − K): with the site features "
           f"{min(ref[c['graph'], c['tumor'], c['platform']]['C'] - ref[c['graph'], c['tumor'], c['platform']]['K'] for c in hg)}–"
           f"{max(ref[c['graph'], c['tumor'], c['platform']]['C'] - ref[c['graph'], c['tumor'], c['platform']]['K'] for c in hg)} per call set, "
           f"without {min(c['C'] - c['K'] for c in hg):,}–{max(c['C'] - c['K'] for c in hg):,}.", ""]
    (HERE / "r26p2_nositef.md").write_text("\n".join(md) + "\n")
    print("\n".join(md))


if __name__ == "__main__":
    main()
