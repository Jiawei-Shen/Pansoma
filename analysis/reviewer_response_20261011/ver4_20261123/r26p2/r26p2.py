#!/usr/bin/env python3
"""ver4, R2.6 point 2: what happens to the germline SNV candidates the graph does not absorb (chr1, benchmark regions),
for the current six-tumor SNV models (seed 0) on HPRC v2.1, HPRC v1.1 and, as a contrast, the GRCh38-only graph.

Replaces the per-tensor chain of E1.1 (experiments/E1.1_ch6_counterfactual/ch6_counterfactual.py `summarize`, base mode)
and the call-level share of E4b (experiments/E4b_fp_decomposition/fp_decomposition.py, class "germline, same allele"),
which the letter's sentences were computed with on the October HG008T single-tumor models. Nothing is predicted here:
every input is an existing file.

Per cell (reference x tumor x platform), from the run's chr1 predictions (<run>/test_chr1/*.predictions.ndjson.gz, the
file named in the graph VCF header) and the call set (<run>/<vcf dir>):
  t      PASS threshold of the call set (<SET>.SNV.graph.stats.json `threshold`: rule 2, the validation threshold at
         recall 0.95 of the callable positives after the population filter)
  G      germline candidates: test tensors (in_test: labelled 0 / 1 / 2 inside the confident region; the SNV partial
         somatic matches left out) whose test_label is 2 (the candidate's allele is in the matched normal's germline
         truth). Model independent.
  G_tag  those the population filter's SNV rule tags (pon_tag 1; site_features.py = scripts/filter_panel_of_normals.py)
  C      G with p_somatic >= t                      -> "the model rejects" = 1 - C / G
  K      C on the GRCh38 path and not tagged       -> "the filter removes" = 1 - K / C ; "remain" = K / G
         (E1.1: keep = ~off_reference & ~pon; an off-reference call has no GRCh38 position and is never scored)
  FPt    somatic calls after the filter on test tensors labelled 0 or 2 (E1.1 truth_after_pon.fp) -> K / FPt
  call level (E4b): the rtg vcfeval FP records of the PASS calls after the filter in the benchmark regions
         (<vcf dir>/eval_bed/vcfeval/pon_calls/fp.vcf.gz) whose allele is carried by the normal's germline truth
         (calls_lib.Germline: every allele of a genotype, any FILTER, normalized; MNPs also as single bases)
         -> share of the false calls after filtering; the same on raw_calls/fp.vcf.gz (before the filter).

usage:  python r26p2.py repro     # reproduces the October numbers (E1.1 report.json, E4b summary.json) with this code
        python r26p2.py           # the 27 cells -> cells.json, r26p2.md
"""
import gzip
import json
import sys
import time
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
EXP = HERE.parents[1] / "experiments"
sys.path.insert(0, str(EXP / "call_level_common"))
import calls_lib as L  # noqa: E402  (read-only helpers: Germline, normalize, read_records, manifest, chrom_seq)

R = Path("/scratch/jshen/data/pansoma_net_v2_runs")
N = R / "snv_newpipe_20261109"
G = R / "m4m6_graphs_20261110"
E11_DATA = Path("/scratch/jshen/data/reviewer_experiments_20261013/E1.1_ch6_counterfactual")
E11_SUMMARY = EXP / "E1.1_ch6_counterfactual/results/summary.json"
E4B_SUMMARY = EXP / "E4b_fp_decomposition/summary.json"
SITEF = Path("/scratch/jshen/data/snv_af05_20261103/site_features")
TUMORS = ["HG008T", "COLO829T", "HCC1395T"]
PLATS = ["Illumina", "PacBio", "ONT"]
GRAPHS = ["d46", "d9", "linear_nositef"]
GRAPH_NAME = {"d46": "HPRC v2.1", "d9": "HPRC v1.1", "linear_nositef": "GRCh38-only"}
RULE = "validation post-PoN recall 0.95 of the callable positives"


def log(*a):
    print(f"[{time.strftime('%H:%M:%S')}]", *a, flush=True)


def ds(tumor, p):
    return "COLO829T_fiberseq" if (tumor, p) == ("COLO829T", "PacBio") else f"{tumor}_{p}"


def pname(tumor, p):
    return {"Illumina": "Illumina", "ONT": "ONT", "PacBio": "Fiber-seq" if tumor == "COLO829T" else "PacBio HiFi"}[p]


def vcf_dir(graph, tumor, p):
    if graph == "d46":
        return N / f"multicancer6_{p}_SNV_newpipe" / f"vcf_chr1_{ds(tumor, p)}_d46"
    return G / f"multicancer6_{p}_SNV_{graph}" / f"vcf_chr1_{tumor}"


def header_value(vcf, key):
    with gzip.open(vcf, "rt") as f:
        for line in f:
            if not line.startswith("##"):
                break
            if line.startswith(f"##{key}="):
                return line.rstrip("\n").split("=", 1)[1]
    raise ValueError(f"{vcf}: no ##{key}")


def load_predictions(path):
    """Arrays over the records of a chr1 predictions file (index order)."""
    in_test, lab, off, pon, p, exc, sub = [], [], [], [], [], [], []
    rescored = 0
    with gzip.open(path, "rt") as f:
        for line in f:
            r = json.loads(line)
            assert r["chrom"] == "chr1", r["chrom"]
            in_test.append(r["in_test"])
            lab.append(-1 if r["test_label"] is None else r["test_label"])
            off.append(r["off_reference"])
            pon.append(r.get("pon_tag", -9))
            p.append(r["p_somatic"])
            exc.append(bool(r.get("excluded")))
            sub.append(r["candidate_id"].rsplit(":", 1)[1])
            rescored += "offref_call" in r
    d = dict(in_test=np.array(in_test, bool), lab=np.array(lab, np.int8), off=np.array(off, bool),
             pon=np.array(pon, np.int8), p=np.array(p, np.float64), exc=np.array(exc, bool), sub=np.array(sub))
    assert rescored == 0, f"{path}: {rescored} off-reference rescored records"
    assert np.array_equal(d["in_test"], d["lab"] >= 0)
    d["has_pon"] = bool((d["pon"] != -9).any())
    assert not d["has_pon"] or np.array_equal(d["off"], d["pon"] == -1), f"{path}: off_reference and pon_tag -1 disagree"
    return d


def file_pon(set_name, d, sitef_dir):
    """pon_tag of the chr1 tensors from the set's site-feature file (chr1 is the first block of the index order).
    Returns (tags, share of on-reference tensors whose candidate REF>ALT equals the file's GRCh38 REF>ALT)."""
    z = np.load(Path(sitef_dir) / f"{set_name}.SNV.npz", allow_pickle=True)
    n = len(d["p"])
    assert str(z["chroms"][0]) == "chr1" and int(z["chrom_counts"][0]) == n, (set_name, n)
    tag = np.asarray(z["pon_tag"][:n], np.int8)
    assert np.array_equal(d["off"], tag == -1), f"{set_name}: off_reference and the file's pon_tag -1 disagree"
    on = ~d["off"]
    ra = np.char.add(np.char.add(np.asarray(z["vcf_ref"][:n]).astype(str), ">"), np.asarray(z["vcf_alt"][:n]).astype(str))
    return tag, float((ra[on] == d["sub"][on]).mean())


def chain(lab, off, pon, p, thr):
    """The per-tensor chain on test tensors (arrays restricted to in_test). pon: bool, tagged."""
    g = lab == 2
    call = p >= thr
    keep = ~off & ~pon
    return dict(
        test_tensors=int(len(lab)), somatic=int((lab == 1).sum()), non=int((lab == 0).sum()),
        G=int(g.sum()), G_offref=int((g & off).sum()), G_tag=int((g & pon).sum()),
        G_untag_onref=int((g & keep).sum()),
        C=int((g & call).sum()), C_tag=int((g & call & pon).sum()), C_offref=int((g & call & off).sum()),
        K=int((g & call & keep).sum()),
        non_called=int(((lab == 0) & call).sum()), non_called_after=int(((lab == 0) & call & keep).sum()),
        somatic_called=int(((lab == 1) & call).sum()), somatic_called_after=int(((lab == 1) & call & keep).sum()),
        FPt=int(((lab != 1) & call & keep).sum()))


_GERM = {}


def germline_same(fp_vcf, set_name):
    """(FP records, those whose allele the normal's germline truth carries) of an rtg fp.vcf.gz."""
    m = L.manifest(set_name, "SNV")
    key = (m["germline_vcf"], m["fasta"])
    if key not in _GERM:
        _GERM[key] = L.germline(*key)
    germ, seq = _GERM[key], L.chrom_seq(m["fasta"])
    n = same = 0
    for rec in L.read_records(fp_vcf, with_format=False):
        keys = [L.normalize(seq, rec.pos, rec.ref, alt) for alt in rec.alts]
        assert any(k is not None for k in keys), (fp_vcf, rec.pos)
        n += 1
        same += any(k is not None and k in germ.keys for k in keys)
    return n, same, m


# ---------------------------------------------------------------------------------------------------- reproduction
def repro():
    """This code on the October single-tumor HG008T models: must give E1.1's and E4b's numbers."""
    e11 = json.loads(E11_SUMMARY.read_text())["runs"]
    out = []
    # (a) per-tensor chain. E1.1 kept its own probabilities (probs_base.npy over positions.npy = index positions of the
    # set); labels, off_reference and the PoN tag come from a current predictions file of the same tensor set.
    for run, pred in (("HG008_PacBio_SNV_d46_7fc7396", N / "multicancer6_PacBio_SNV_newpipe/test_chr1/HG008T_PacBio_d46.tensors.SNV.predictions.ndjson.gz"),
                      ("HG008_ONT_SNV_d46_7fc7396", N / "multicancer6_ONT_SNV_newpipe/test_chr1/HG008T_ONT_d46.tensors.SNV.predictions.ndjson.gz"),
                      ("HG008_PacBio_SNV_7fc7396", G / "multicancer6_PacBio_SNV_d9/test_chr1/HG008T_PacBio.tensors.SNV.predictions.ndjson.gz"),
                      ("HG008_ONT_SNV_7fc7396", G / "multicancer6_ONT_SNV_d9/test_chr1/HG008T_ONT.tensors.SNV.predictions.ndjson.gz")):
        old = e11[run]
        set_name = old["set"]
        assert f"/{set_name}.tensors." in str(pred), (set_name, pred)
        z = np.load(SITEF / f"{set_name}.SNV.npz", allow_pickle=True)
        assert str(z["chroms"][0]) == "chr1"  # chr1 is the first block of the index: index position = chr1 line
        d = load_predictions(pred)
        assert len(d["p"]) == int(z["chrom_counts"][0]) == old["chr1_tensors"]
        pos = np.load(E11_DATA / run / "positions.npy")
        probs = np.load(E11_DATA / run / "probs_base.npy")
        assert len(pos) == len(probs) == old["evaluated"]
        chr1 = pos < len(d["p"])  # --assembly runs also scored a few keys of other chromosomes
        test = np.zeros(len(pos), bool)
        test[chr1] = d["in_test"][pos[chr1]]
        assert int(test.sum()) == old["test_tensors"] == int(d["in_test"].sum()), (run, int(test.sum()))
        tp = pos[test]
        tags = json.loads((E11_DATA / run / "pon_tags.json").read_text())
        old_pon = np.array([bool(tags.get(str(int(q)))) for q in tp])
        new_pon = d["pon"][tp] == 1
        diff = int((old_pon != new_pon).sum())
        c = chain(d["lab"][tp], d["off"][tp], new_pon, probs[test][:, 1].astype(np.float64),
                  old["threshold_val_recall_0_9"])
        b = old["results"]["base"]
        want = dict(G=b["germline"]["n"], C=b["germline"]["ge_thr"], K=b["germline_after_pon_ge_thr"],
                    FPt=b["truth_after_pon"]["fp"])
        got = {k: c[k] for k in want}
        ok = got == want
        out.append(f"E1.1 {run} ({set_name}): this code {got} | E1.1 summary.json {want} | "
                   f"{'IDENTICAL' if ok else 'DIFFERENT'}; PoN tag (site-feature file vs E1.1 pon_tags.json) differs on "
                   f"{diff} of {len(tp):,} test tensors")
        log(out[-1])
        assert ok and diff == 0, run
    # (b) call level: E4b's "germline, same allele" on the October call sets
    for cs in json.loads(E4B_SUMMARY.read_text())["call_sets"]:
        if not (cs["eval"] == "eval_bed" and cs["kind"] == "SNV" and cs["group"] == "insample"
                and cs["sample"] == "HG008T" and cs["variant"] == "main"):
            continue
        fp = Path(cs["vcf_dir"]) / "eval_bed/vcfeval/pon_calls/fp.vcf.gz"
        n, same, _ = germline_same(fp, cs["set"])
        want = (cs["fp_rows"], cs["classes"].get("germline_same", 0))
        ok = (n, same) == want
        out.append(f"E4b {cs['id']} ({cs['arm']}): this code {same} / {n} | E4b summary.json {want[1]} / {want[0]} | "
                   f"{'IDENTICAL' if ok else 'DIFFERENT'}")
        log(out[-1])
        assert ok, cs["id"]
    (HERE / "repro.txt").write_text("\n".join(out) + "\n")
    log(f"reproduction: {len(out)} comparisons identical -> repro.txt")


# ---------------------------------------------------------------------------------------------------- the cells
def cell(graph, tumor, p):
    v = vcf_dir(graph, tumor, p)
    gstats = [x for x in v.glob("*.SNV.graph.stats.json")]
    assert len(gstats) == 1, v
    st = json.loads(gstats[0].read_text())
    set_name = gstats[0].name.split(".SNV.")[0]
    gvcf = v / f"{set_name}.SNV.graph.vcf.gz"
    pred = Path(header_value(gvcf, "pansoma_predictions"))
    assert pred.parent.parent == v.parent and f"/{set_name}.tensors.SNV." in str(pred), (pred, v)
    assert Path(header_value(gvcf, "pansoma_tensors")).parent.name == set_name
    thr = st["threshold"]
    assert RULE in st["call"] and st["offref_rescored"] == 0 and st["chromosomes"] == ["chr1"], st["call"]
    rep = json.loads((v / "eval_bed/report.json").read_text())
    assert rep["threshold"] == thr and rep["kind"] == "SNV" and rep["chromosomes"] == ["chr1"]
    met = json.loads(Path(str(pred).replace(".predictions.ndjson.gz", ".metrics.json")).read_text())
    assert Path(met["checkpoint"]) == v.parent / "best.pth" and Path(met["directory"]).parent.parent.name == set_name
    args = json.loads((v.parent / "args.json").read_text())
    assert args["seed"] == 0 and args["call_recall"] == 0.95 and args["best_by"] == "postpon"
    assert bool(args.get("site_features_pon_only")) == (graph == "linear_nositef")

    d = load_predictions(pred)
    assert len(d["p"]) == st["predictions"] == met["all_tensors"]
    tag, same_sub = file_pon(set_name, d, args["site_features"])
    assert d["has_pon"] == (graph != "linear_nositef")
    if d["has_pon"]:  # the predictions carry the tag: it must be the file's, tensor by tensor
        assert np.array_equal(tag, d["pon"]), f"{set_name}: pon_tag of the predictions != site-feature file"
    else:  # GRCh38-only models do not read the features: the tag comes from the file; every node is forward
        assert same_sub == 1.0 and not d["off"].any(), (set_name, same_sub)
        d["pon"] = tag
    assert int(((d["p"] >= thr) & ~d["exc"]).sum()) == st["called"], "PASS count of graph_vcf"
    t = d["in_test"]
    assert int(t.sum()) == met["tensors"]
    c = chain(d["lab"][t], d["off"][t], d["pon"][t] == 1, d["p"][t], thr)
    conf = met["argmax"]["confusion"]  # rows = truth non / somatic / germline
    assert [sum(r) for r in conf] == [c["non"], c["somatic"], c["G"]], "class supports of predict's metrics.json"
    if graph != "linear_nositef":
        # predict's own post-filter rule at this threshold (a tagged tensor is never called; off-reference ones are):
        assert met["rules"]["val_postpon_recall"]["threshold"] == thr
        conf2 = met["rules"]["val_postpon_recall"]["confusion"]
        assert conf2[2][1] == c["K"] + c["C_offref"] and conf2[0][1] == int(((d["lab"][t] == 0) & (d["p"][t] >= thr) & (d["pon"][t] != 1)).sum())
        assert met["pon_tagged"] == int((d["pon"] == 1).sum())  # predict counts the tagged tensors of all of chr1

    out = dict(graph=graph, reference=GRAPH_NAME[graph], tumor=tumor, platform=pname(tumor, p), set=set_name,
               run=str(v.parent), vcf_dir=str(v), predictions=str(pred), threshold=thr, epoch=met_epoch(v.parent),
               site_features_input=graph != "linear_nositef", candidate_ref_alt_equals_grch38=same_sub, **c)
    for stage in ("raw", "pon"):
        a = rep[stage]["at_threshold"]
        n, same, m = germline_same(v / f"eval_bed/vcfeval/{stage}_calls/fp.vcf.gz", set_name)
        assert n == a["fp"], (v, stage, n, a["fp"])
        assert a["tp_baseline"] + a["fn"] == rep[stage]["baseline"]
        out[stage] = dict(truth=a["tp_baseline"] + a["fn"], tp=a["tp_baseline"], fp=a["fp"], fn=a["fn"],
                          precision=a["precision"], recall=a["recall"], f1=a["f1"], fp_germline_same=same)
    assert out["raw"]["truth"] == out["pon"]["truth"]
    out.update(unplaced_pass=rep["unplaced_pass"], germline_vcf=m["germline_vcf"], truth_vcf=m["truth_vcf"],
               truth_bed=m["truth_bed"])
    assert rep["truth"] == m["truth_vcf"] and rep["bed"] == m["truth_bed"]
    log(f"{graph} {tumor} {p} ({set_name}): G {c['G']:,} tagged {c['G_tag']:,} C {c['C']:,} K {c['K']:,} FPt {c['FPt']:,} | "
        f"rtg FP after {out['pon']['fp']:,} germline same {out['pon']['fp_germline_same']:,}")
    return out


def met_epoch(run):
    import re
    m = re.findall(r"^done; best epoch (\d+)", (run / "train.log").read_text(), re.M)
    return int(m[-1])


def pct(n, d, k=1):
    return "-" if not d else f"{100 * n / d:.{k}f}%"


def rng(vals, k=1):
    lo, hi = min(vals), max(vals)
    return f"{100 * lo:.{k}f}%" if f"{100 * lo:.{k}f}" == f"{100 * hi:.{k}f}" else f"{100 * lo:.{k}f}–{100 * hi:.{k}f}%"


QUANT = {  # name -> (numerator, denominator) functions of a cell
    "rejects": lambda c: (c["G"] - c["C"], c["G"]),
    "tagged": lambda c: (c["G_tag"], c["G"]),
    "called": lambda c: (c["C"], c["G"]),
    "removes": lambda c: (c["C"] - c["K"], c["C"]),
    "remain": lambda c: (c["K"], c["G"]),
    "untag_called": lambda c: (c["K"], c["G_untag_onref"]),
    "tag_called": lambda c: (c["C_tag"], c["G_tag"]),
    "share_rtg": lambda c: (c["pon"]["fp_germline_same"], c["pon"]["fp"]),
    "share_tensor": lambda c: (c["K"], c["FPt"]),
    "share_rtg_raw": lambda c: (c["raw"]["fp_germline_same"], c["raw"]["fp"]),
    "removes_rtg": lambda c: (c["raw"]["fp_germline_same"] - c["pon"]["fp_germline_same"], c["raw"]["fp_germline_same"]),
}


def frac(c, q):
    n, d = QUANT[q](c)
    return n / d


def table(cells, graph):
    head = ["tumor", "platform", "germline SNV candidates, chr1 (G)", "tagged by the filter's rule",
            "called somatic at the PASS threshold (C)", "rejected by the model, 1 − C/G",
            "left after the filter (K)", "removed by the filter, 1 − K/C", "remaining, K/G",
            "false calls after the filter, chr1 (rtg FP)", "of them the normal's germline allele", "share"]
    rows = ["| " + " | ".join(head) + " |", "|" + "---|" * len(head)]
    for c in cells:
        if c["graph"] != graph:
            continue
        rows.append("| " + " | ".join([
            c["tumor"], c["platform"], f"{c['G']:,}", f"{c['G_tag']:,} ({pct(c['G_tag'], c['G'])})",
            f"{c['C']:,}", pct(c["G"] - c["C"], c["G"]), f"{c['K']:,}", pct(c["C"] - c["K"], c["C"]),
            pct(c["K"], c["G"]), f"{c['pon']['fp']:,}", f"{c['pon']['fp_germline_same']:,}",
            pct(c["pon"]["fp_germline_same"], c["pon"]["fp"])]) + " |")
    return rows


def table2(cells, graph):
    head = ["tumor", "platform", "PASS threshold (epoch)", "tagged candidates called somatic",
            "untagged candidates (on GRCh38 path)", "of them called somatic (= K)",
            "per-tensor false calls after the filter (FPt)", "germline share K/FPt",
            "rtg FP before the filter, chr1", "of them germline allele", "germline FP records removed by the filter (rtg)",
            "chr1 truth SNVs", "chr1 P / R / F1 after the filter"]
    rows = ["| " + " | ".join(head) + " |", "|" + "---|" * len(head)]
    for c in cells:
        if c["graph"] != graph:
            continue
        a, r = c["pon"], c["raw"]
        rows.append("| " + " | ".join([
            c["tumor"], c["platform"], f"{c['threshold']:.4f} ({c['epoch']})",
            f"{c['C_tag']:,} of {c['G_tag']:,} ({pct(c['C_tag'], c['G_tag'], 2)})", f"{c['G_untag_onref']:,}",
            f"{c['K']:,} ({pct(c['K'], c['G_untag_onref'])})", f"{c['FPt']:,}", pct(c["K"], c["FPt"]),
            f"{r['fp']:,}", f"{r['fp_germline_same']:,} ({pct(r['fp_germline_same'], r['fp'])})",
            f"{r['fp_germline_same'] - a['fp_germline_same']:,} ({pct(r['fp_germline_same'] - a['fp_germline_same'], r['fp_germline_same'])})",
            f"{a['truth']:,}", f"{a['precision']:.3f} / {a['recall']:.3f} / {a['f1']:.3f}"]) + " |")
    return rows


def old_rows():
    """The October HG008T single-tumor numbers the letter's sentences were written from (copied from the old files)."""
    e11 = json.loads(E11_SUMMARY.read_text())["runs"]
    e4b = {(c["arm"], c["platform"]): c for c in json.loads(E4B_SUMMARY.read_text())["call_sets"]
           if c["eval"] == "eval_bed" and c["kind"] == "SNV" and c["group"] == "insample" and c["sample"] == "HG008T"
           and c["variant"] == "main"}
    runs = {("d9", "Illumina"): "HG008_Illumina_SNV_7fc7396_ddp4_b1024", ("d9", "PacBio"): "HG008_PacBio_SNV_7fc7396",
            ("d9", "ONT"): "HG008_ONT_SNV_7fc7396", ("d46", "Illumina"): "HG008_Illumina_SNV_d46_7fc7396",
            ("d46", "PacBio"): "HG008_PacBio_SNV_d46_7fc7396", ("d46", "ONT"): "HG008_ONT_SNV_d46_7fc7396"}
    out = []
    for (g, p), run in runs.items():
        b = e11[run]["results"]["base"]
        f = e4b[(g, p)]
        out.append(dict(graph=g, reference=GRAPH_NAME[g], tumor="HG008T", platform=pname("HG008T", p), run=run,
                        G=b["germline"]["n"], C=b["germline"]["ge_thr"], K=b["germline_after_pon_ge_thr"],
                        FPt=b["truth_after_pon"]["fp"], sampled="sampling" in e11[run],
                        pon=dict(fp=f["fp_rows"], fp_germline_same=f["classes"].get("germline_same", 0))))
    return out


def main():
    cells = [cell(g, t, p) for g in GRAPHS for t in TUMORS for p in PLATS]
    # truth counts: equal across the three references of a tumor x platform, before and after the filter
    for t in TUMORS:
        for p in PLATS:
            n = {c["raw"]["truth"] for c in cells if (c["tumor"], c["platform"]) == (t, pname(t, p))}
            assert len(n) == 1, (t, p, n)
    old = old_rows()
    by = lambda g, t=None: [c for c in cells if c["graph"] == g and (t is None or c["tumor"] == t)]  # noqa: E731
    oby = lambda g: [c for c in old if c["graph"] == g]  # noqa: E731
    ranges = {}
    for g in GRAPHS:
        for scope, sel in (("HG008T", by(g, "HG008T")), ("COLO829T", by(g, "COLO829T")), ("HCC1395T", by(g, "HCC1395T")),
                           ("all nine", by(g))):
            ranges[f"{g}|{scope}"] = {q: [min(frac(c, q) for c in sel), max(frac(c, q) for c in sel)] for q in QUANT}
    for scope, sel in (("HG008T", by("d46", "HG008T") + by("d9", "HG008T")), ("all nine", by("d46") + by("d9"))):
        ranges[f"d9+d46|{scope}"] = {q: [min(frac(c, q) for c in sel), max(frac(c, q) for c in sel)] for q in QUANT}
    oq = ("rejects", "removes", "remain", "share_rtg", "share_tensor")
    old_ranges = {g: {q: [min(frac(c, q) for c in oby(g)), max(frac(c, q) for c in oby(g))] for q in oq} for g in ("d9", "d46")}
    old_ranges["d9+d46"] = {q: [min(frac(c, q) for c in old), max(frac(c, q) for c in old)] for q in oq}
    (HERE / "cells.json").write_text(json.dumps(dict(
        generated=time.strftime("%Y-%m-%d %H:%M"), script=str(Path(__file__).resolve()), cells=cells, ranges=ranges,
        october_single_tumor=old, october_ranges=old_ranges), indent=1) + "\n")

    def r(key, q, k=1):
        return rng(ranges[key][q], k)

    def ro(g, q, k=1):
        return rng(old_ranges[g][q], k)

    md = ["# ver4, R2.6 point 2: germline SNV candidates the graph does not absorb, six-tumor models (chr1)", "",
          f"Generated by `r26p2.py` ({time.strftime('%Y-%m-%d %H:%M')}); every number is read from files (cells.json holds "
          "all of them). Models: the current six-tumor SNV models, seed 0 "
          "(HPRC v2.1 `snv_newpipe_20261109/multicancer6_<P>_SNV_newpipe`; HPRC v1.1 and GRCh38-only "
          "`m4m6_graphs_20261110/multicancer6_<P>_SNV_{d9,linear_nositef}`); chr1 of HG008T, COLO829T and HCC1395T is "
          "held out of them. PASS = rule 2 (p_somatic >= the validation threshold at recall 0.95 of the callable "
          "positives after the population filter). chr1, benchmark regions.", "",
          "Definitions. **G**: chr1 SNV candidates (tensors) inside the benchmark's confident regions whose allele is in "
          "the matched normal's germline truth (label germline; HG008-N and COLO829BL dipcall, HCC1395BL DeepVariant); "
          "these are the germline variants the graph did not absorb. Model independent. **tagged**: G whose GRCh38 allele "
          "the population filter's SNV rule tags (gnomAD or CoLoRSdb AF >= 1e-4, dbSNP non-somatic, 1000G panel). "
          "**C**: G with p_somatic >= the PASS threshold. **K**: C that the filter does not tag (and that have a GRCh38 "
          "position). **rtg FP**: false-positive records of rtg vcfeval for the PASS calls after the filter, benchmark "
          "regions; **germline allele**: the call's allele is carried by the matched normal's germline truth (any "
          "FILTER). Tensor sets: HG008T and all ONT sets at build AF 0.08; COLO829T and HCC1395T Illumina / PacBio HiFi / "
          "Fiber-seq at AF 0.05 (`_snvaf05`), so their G is larger than on the AF 0.08 sets.", ""]
    for i, g in enumerate(GRAPHS, 1):
        note = {"d46": "population site features are a model input",
                "d9": "population site features are a model input",
                "linear_nositef": "contrast only: no population site features as a model input; on GRCh38-only every "
                                  "germline variant is a candidate, so G here is not \"absent from the graph\""}[g]
        md += [f"## {i}. {GRAPH_NAME[g]} six-tumor models, chr1 ({note})", ""] + table(cells, g) + [""]
    md += ["## 4. Detail per cell, chr1 (threshold, tagged / untagged split, per-tensor share, before the filter)", "",
           "`tagged candidates called somatic`: C among the tagged G. `untagged candidates`: G on the GRCh38 path that "
           "the filter does not tag (rare germline alleles absent from the graph and from the resources); K is the part "
           "of them the model calls somatic. `FPt`: PASS calls after the filter on test tensors labelled germline or "
           "non-variant (the E1.1 per-tensor count). `rtg FP before the filter`: `raw_calls`.", ""]
    for g in GRAPHS:
        md += [f"### {GRAPH_NAME[g]}, chr1", ""] + table2(cells, g) + [""]

    md += ["## 5. October single-tumor HG008T models (what the ver3 sentences were computed from), chr1", "",
           "Copied from `experiments/E1.1_ch6_counterfactual/results/summary.json` (base mode; PASS then = the validation "
           "recall-0.9 threshold, no site features as input; Illumina = weighted estimate of a class-stratified sample) "
           "and `experiments/E4b_fp_decomposition/summary.json` (eval_bed, main call sets).", "",
           "| reference | platform | G | C | rejected by the model | K | removed by the filter | remaining, K/G | rtg FP after the filter | germline allele | share |",
           "|---|---|---|---|---|---|---|---|---|---|---|"]
    for c in old:
        md.append(f"| {c['reference']} | {c['platform']} | {c['G']:,} | {c['C']:,} | {pct(c['G'] - c['C'], c['G'])} | "
                  f"{c['K']:,} | {pct(c['C'] - c['K'], c['C'])} | {pct(c['K'], c['G'])} | {c['pon']['fp']:,} | "
                  f"{c['pon']['fp_germline_same']:,} | {pct(c['pon']['fp_germline_same'], c['pon']['fp'])} |")
    md += ["", "ver3 ranges, HG008T chr1, v1.1 / v2.1 / both: rejected by the model "
           f"{ro('d9', 'rejects')} / {ro('d46', 'rejects')} / {ro('d9+d46', 'rejects')}; removed by the filter "
           f"{ro('d9', 'removes')} / {ro('d46', 'removes')} / {ro('d9+d46', 'removes')}; remaining "
           f"{ro('d9', 'remain')} / {ro('d46', 'remain')} / {ro('d9+d46', 'remain')}; share of the false calls after "
           f"filtering (rtg) {ro('d9', 'share_rtg')} / {ro('d46', 'share_rtg')} / {ro('d9+d46', 'share_rtg')}. "
           "Each ver3 range (79–88%, 74–89%, 2.4–3.2%, 83–98%) spans both graphs.", ""]

    md += ["## 6. Ranges for the sentences (six-tumor models, chr1)", "",
           "| quantity | HG008T v2.1, chr1 | HG008T v1.1, chr1 | HG008T both graphs, chr1 | nine datasets v2.1, chr1 | nine datasets v1.1, chr1 | COLO829T v2.1, chr1 | HCC1395T v2.1, chr1 | ver3 (October HG008T, both graphs), chr1 |",
           "|---|---|---|---|---|---|---|---|---|"]
    names = [("rejects", "the model rejects (1 − C/G)", 1), ("removes", "the filter removes of those called somatic (1 − K/C)", 1),
             ("remain", "remain after classification and filtering (K/G)", 1),
             ("share_rtg", "germline alleles among the false calls after filtering (rtg FP)", 1),
             ("share_tensor", "the same per tensor (K/FPt)", 1),
             ("tagged", "germline candidates the filter's rule tags", 1),
             ("tag_called", "tagged germline candidates called somatic", 2),
             ("untag_called", "untagged germline candidates called somatic", 1),
             ("removes_rtg", "germline FP records removed by the filter (rtg, call level)", 1)]
    for q, label, k in names:
        md.append(f"| {label} | {r('d46|HG008T', q, k)} | {r('d9|HG008T', q, k)} | {r('d9+d46|HG008T', q, k)} | "
                  f"{r('d46|all nine', q, k)} | {r('d9|all nine', q, k)} | {r('d46|COLO829T', q, k)} | "
                  f"{r('d46|HCC1395T', q, k)} | {ro('d9+d46', q, k) if q in oq else '-'} |")
    md += ["", "GRCh38-only six-tumor models (no population site features as input), chr1, HG008T / nine datasets: the "
           f"model rejects {r('linear_nositef|HG008T', 'rejects')} / {r('linear_nositef|all nine', 'rejects')} of the "
           f"germline candidates, the filter removes {r('linear_nositef|HG008T', 'removes')} / "
           f"{r('linear_nositef|all nine', 'removes')} of those called somatic.", ""]
    md += ["## 7. Candidate level (model independent): germline SNV candidates per reference, chr1", "",
           "G of sections 1-3 side by side: the chr1 germline SNV candidates in the benchmark regions on each reference "
           "(same reads, same AF floor per row), and the share of the GRCh38-only count that is still a candidate on "
           "the graph.", "",
           "| tumor | platform | GRCh38-only, chr1 | HPRC v1.1, chr1 (of GRCh38-only) | HPRC v2.1, chr1 (of GRCh38-only) |",
           "|---|---|---|---|---|"]
    for t in TUMORS:
        for pl in PLATS:
            g = {c["graph"]: c["G"] for c in cells if (c["tumor"], c["platform"]) == (t, pname(t, pl))}
            md.append(f"| {t} | {pname(t, pl)} | {g['linear_nositef']:,} | {g['d9']:,} ({pct(g['d9'], g['linear_nositef'])}) | "
                      f"{g['d46']:,} ({pct(g['d46'], g['linear_nositef'])}) |")
    md.append("")
    (HERE / "r26p2.md").write_text("\n".join(md) + "\n")
    log("wrote cells.json, r26p2.md")


if __name__ == "__main__":
    if len(sys.argv) > 1 and sys.argv[1] == "repro":
        repro()
    else:
        main()
