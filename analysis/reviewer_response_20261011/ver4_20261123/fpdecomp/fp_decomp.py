#!/usr/bin/env python3
"""ver4 (2026-11-23) of Response Tables 3.8 (= 2.6) and 3.9: the GRCh38-only model's PASS false positives split by
whether the allele is a germline allele of the matched normal and whether the HPRC v2.1 graph carries it, and what
happens to the same allele on HPRC v2.1 with the same reads. Modified copy of
tmp/r32_prepon_20261104/r32_prepon.py (kept unchanged next to this file as r32_prepon_20261104_original.py).

Changes against the original:
- the call sets are given as VCF directories (not L.RUNS/<run>/vcf_chr1); the tensor set of each arm is read from the
  directory's <set>.<KIND>.graph.vcf.gz, so the v2.1 "candidate retained" lookup uses the test set of the v2.1 model
  (the _d46_snvaf05 sets for SNV on COLO829T Illumina / Fiber-seq and HCC1395T Illumina / PacBio HiFi);
- --mode ver4: the six-tumor current models, HG008T / COLO829T / HCC1395T x Illumina / PacBio HiFi (COLO829T
  Fiber-seq) / ONT x SNV / INDEL, seed 0:
    GRCh38-only  R/m4m6_graphs_20261110/multicancer6_<P>_<KIND>_linear_nositef/vcf_[ppr8_|pp_]chr1_<TUMOR>
    HPRC v2.1    SNV R/snv_newpipe_20261109/multicancer6_<P>_SNV_newpipe/vcf_chr1_<SET>_d46,
                 INDEL R/release_s0_<P>_INDEL/vcf_[ppr8_|pp_]chr1_<SET>_d46   (ONT INDEL: pp)
- the decomposition is done twice: on the false positives before the population filter (raw_calls/fp.vcf.gz) and
  after it (pon_calls/fp.vcf.gz);
- --mode october: cells of the 2026-10-12 single-tumor runs (reviewer_cpu_20261013/<S>_<P>_<K>_{linear,d46}_7fc7396),
  compared field by field with the old r32_prepon.tsv (the reproduction check);
- checks: truth count (TP + FN) equal in the two arms, equal to report.json "baseline", and equal to the records of
  tp-baseline.vcf.gz + fn.vcf.gz; FP records of fp.vcf.gz equal to report.json fp; same truth VCF / BED / germline VCF
  in the two arms;
- only the v2.1 tensor keys that are GRCh38-only false positives are kept in memory.

Unchanged: chr1, benchmark BED (eval_bed), alleles normalized with calls_lib.normalize, germline = an allele carried
by a genotype of the set's germline truth VCF (calls_lib.Germline keys; HG008-N and COLO829BL dipcall, HCC1395BL
DeepVariant QUAL >= 10), HPRC v2.1 carriers among the 463 haplotypes (E7 v2.1 hprc_alleles chr1, n_carriers), >= 46 =
"in the HPRC v2.1 graph" (the d46 threshold), candidate on v2.1 = a chr1 tensor of the v2.1 set whose GRCh38
projection is the same allele.

Writes into --out (default: this directory): fp_decomp.tsv (long), fp_decomp.json (per cell), tables.md, ranges.md."""
import argparse
import collections
import gzip
import json
import sys
from pathlib import Path

sys.path.insert(0, "/scratch/jshen/Github/Pansoma/analysis/reviewer_response_20261011/experiments/call_level_common")
import calls_lib as L

HERE = Path(__file__).resolve().parent
R = Path("/scratch/jshen/data/pansoma_net_v2_runs")
V21 = Path("/scratch/jshen/data/reviewer_experiments_20261013/E7_graph_composition/v2.1/hprc_alleles/chr1.tsv.gz")
OLD_TSV = Path("/scratch/jshen/Github/Pansoma/tmp/r32_prepon_20261104/r32_prepon.tsv")
TUMORS = ("HG008T", "COLO829T", "HCC1395T")
PLATFORMS = ("Illumina", "PacBio", "ONT")
KINDS = ("SNV", "INDEL")
STAGES = {"raw": "raw_calls", "pon": "pon_calls"}
FIELDS = ("n", "pon_removes", "d46_tensor", "d46_raw_fp", "d46_pon_fp")


def ds(tumor, p):
    return "fiberseq" if (tumor, p) == ("COLO829T", "PacBio") else p


def pname(tumor, p):
    return {"Illumina": "Illumina", "ONT": "ONT", "PacBio": "Fiber-seq" if tumor == "COLO829T" else "PacBio HiFi",
            "fiberseq": "Fiber-seq"}[p]


def ver4_cells():
    out = []
    for kind in KINDS:
        for tumor in TUMORS:
            for p in PLATFORMS:
                r = "" if kind == "SNV" else ("pp_" if p == "ONT" else "ppr8_")
                lin = R / "m4m6_graphs_20261110" / f"multicancer6_{p}_{kind}_linear_nositef" / f"vcf_{r}chr1_{tumor}"
                if kind == "SNV":
                    d46 = R / "snv_newpipe_20261109" / f"multicancer6_{p}_SNV_newpipe" / f"vcf_chr1_{tumor}_{ds(tumor, p)}_d46"
                else:
                    d46 = R / f"release_s0_{p}_INDEL" / f"vcf_{r}chr1_{tumor}_{ds(tumor, p)}_d46"
                out.append(dict(tumor=tumor, platform=p, kind=kind, lin_dir=lin, d46_dir=d46))
    return out


def october_cells(spec):
    out = []
    for item in spec.split(","):
        tumor, p, kind = item.split(":")
        s = tumor[:-1]
        out.append(dict(tumor=tumor, platform=p, kind=kind,
                        lin_dir=L.RUNS / f"{s}_{p}_{kind}_linear_7fc7396" / "vcf_chr1",
                        d46_dir=L.RUNS / f"{s}_{p}_{kind}_d46_7fc7396" / "vcf_chr1"))
    return out


def set_of(vcf_dir, kind):
    g = sorted(Path(vcf_dir).glob(f"*.{kind}.graph.vcf.gz"))
    assert len(g) == 1, (vcf_dir, g)
    return g[0].name[: -len(f".{kind}.graph.vcf.gz")]


def n_records(path):
    n = 0
    with gzip.open(path, "rt") as f:
        for line in f:
            n += line[0] != "#"
    return n


def fp_keys(vcf_dir, which, seq):
    """(normalized allele keys, records) of one rtg run's fp.vcf.gz (as the original fp_keys)."""
    f = Path(vcf_dir) / "eval_bed" / "vcfeval" / which / "fp.vcf.gz"
    keys = set()
    recs = L.read_records(f, with_format=False)
    for r in recs:
        for a in r.alts:
            n = L.normalize(seq, r.pos, r.ref, a)
            if n:
                keys.add(n)
    return keys, len(recs)


def tensor_keys(set_name, kind, seq, keep=None):
    """Normalized GRCh38 alleles of the set's chr1 tensors (as the original); with `keep`, only the keys in it."""
    keys = set()
    n = 0
    for line in open(L.TENSORS / set_name / "tensors" / kind / "chr1_labels.ndjson"):
        g = json.loads(line).get("grch38")
        if not g or g["chrom"] != "chr1":
            continue
        n += 1
        p0, ref, alt = g["pos0"], g["ref"], g["alt"]
        if len(ref) == len(alt):
            k = L.normalize(seq, p0 + 1, ref, alt)
        else:  # INDEL: pad with the base before pos0
            b = seq[p0 - 1]
            k = L.normalize(seq, p0, b + ref, b + alt)
        if k and (keep is None or k in keep):
            keys.add(k)
    return keys, n


def check_reports(c):
    """Truth counts and FP counts of the two arms against report.json and the rtg files."""
    out = {}
    reps = {arm: json.loads((c[f"{arm}_dir"] / "eval_bed" / "report.json").read_text()) for arm in ("lin", "d46")}
    for key in ("truth", "bed", "truth_records"):
        assert reps["lin"][key] == reps["d46"][key], (c["lin_dir"], key, reps["lin"][key], reps["d46"][key])
    assert reps["lin"]["pon_files"] == reps["d46"]["pon_files"], c["lin_dir"]
    for stage, which in STAGES.items():
        truth = set()
        for arm in ("lin", "d46"):
            a = reps[arm][stage]["at_threshold"]
            base = reps[arm][stage]["baseline"]
            assert a["tp_baseline"] + a["fn"] == base, (c[f"{arm}_dir"], stage, a, base)
            d = c[f"{arm}_dir"] / "eval_bed" / "vcfeval" / which
            assert (d / "done").exists(), d
            n_tp, n_fn = n_records(d / "tp-baseline.vcf.gz"), n_records(d / "fn.vcf.gz")
            assert (n_tp, n_fn) == (a["tp_baseline"], a["fn"]), (d, n_tp, n_fn, a)
            truth.add(base)
            out[f"{arm}_{stage}"] = dict(tp=a["tp_baseline"], fn=a["fn"], fp=a["fp"], truth=base,
                                         precision=a["precision"], recall=a["recall"], f1=a["f1"])
        assert len(truth) == 1, (c["lin_dir"], stage, truth)
    for arm in ("lin", "d46"):
        out[f"{arm}_call_rule"] = reps[arm].get("steps", {}).get("graph_vcf", {}).get("call")
        out[f"{arm}_threshold"] = reps[arm].get("threshold")
        out[f"{arm}_checkpoint"] = reps[arm].get("steps", {}).get("graph_vcf", {}).get("checkpoint")
    out["truth_vcf"], out["truth_bed"] = reps["lin"]["truth"], reps["lin"]["bed"]
    return out


def compute(cells):
    need = set()
    for c in cells:
        t, p, kind = c["tumor"], c["platform"], c["kind"]
        c["lin_set"], c["d46_set"] = set_of(c["lin_dir"], kind), set_of(c["d46_dir"], kind)
        assert c["lin_set"].startswith(f"{t}_{ds(t, p)}_linear") and c["lin_set"].replace("_linear", "_d46") == c["d46_set"], c
        c["rep"] = check_reports(c)
        m, m46 = L.manifest(c["lin_set"], kind), L.manifest(c["d46_set"], kind)
        for key in ("truth_vcf", "truth_bed", "germline_vcf", "fasta"):
            assert m[key] == m46[key], (c["lin_set"], key, m[key], m46[key])
        assert m["truth_vcf"] == c["rep"]["truth_vcf"], (m["truth_vcf"], c["rep"]["truth_vcf"])
        c["germline_vcf"] = m["germline_vcf"]
        seq = L.chrom_seq(m["fasta"])
        c["g"] = L.germline(m["germline_vcf"], m["fasta"])
        for arm in ("lin", "d46"):
            for stage, which in STAGES.items():
                keys, recs = fp_keys(c[f"{arm}_dir"], which, seq)
                assert recs == c["rep"][f"{arm}_{stage}"]["fp"], (c[f"{arm}_dir"], stage, recs, c["rep"][f"{arm}_{stage}"])
                c[f"{arm}_{stage}"], c[f"{arm}_{stage}_records"] = keys, recs
        own = c["lin_raw"] | c["lin_pon"]
        c["tk"], c["d46_tensors_chr1"] = tensor_keys(c["d46_set"], kind, seq, keep=own)
        need |= own
        L.log(t, p, kind, c["lin_set"], c["d46_set"], "GRCh38-only FP raw / pon", len(c["lin_raw"]), len(c["lin_pon"]),
              "v2.1 FP raw / pon", len(c["d46_raw"]), len(c["d46_pon"]), "v2.1 tensors", c["d46_tensors_chr1"],
              "of them GRCh38-only FP alleles", len(c["tk"]))
    car = {}
    with gzip.open(V21, "rt") as fh:
        head = fh.readline().rstrip("\n").split("\t")
        assert head[0] == "key" and head[6] == "n_carriers", head
        for line in fh:
            k, rest = line.split("\t", 1)
            parts = k.split(":")
            if len(parts) != 4:  # symbolic / breakend ALTs
                continue
            key = (int(parts[1]), parts[2], parts[3])
            if key in need:
                car[key] = int(rest.split("\t")[5])
    L.log("v2.1 carrier keys found", len(car), "of", len(need))
    rows = []
    for c in cells:
        for stage in STAGES:
            cnt = collections.Counter()
            for k in c[f"lin_{stage}"]:
                germ = "germline" if k in c["g"].keys else "other"
                n = car.get(k, 0)
                grp = "ge46" if n >= 46 else "1to45" if n > 0 else "0"
                for cls in (f"{germ}|{grp}", f"{germ}|all", "all|all"):
                    cnt[cls, "n"] += 1
                    cnt[cls, "pon_removes"] += k not in c["lin_pon"]
                    cnt[cls, "d46_tensor"] += k in c["tk"]
                    cnt[cls, "d46_raw_fp"] += k in c["d46_raw"]
                    cnt[cls, "d46_pon_fp"] += k in c["d46_pon"]
            for cls in sorted({x for x, _ in cnt}):
                rows.append(dict(tumor=c["tumor"], platform=c["platform"], kind=c["kind"], stage=stage, cls=cls,
                                 **{f: cnt[cls, f] for f in FIELDS},
                                 linear_raw_fp=len(c["lin_raw"]), linear_pon_fp=len(c["lin_pon"]),
                                 d46_raw_fp_all=len(c["d46_raw"]), d46_pon_fp_all=len(c["d46_pon"])))
    return rows


def cells_json(cells, rows):
    by = collections.defaultdict(dict)
    for r in rows:
        by[(r["tumor"], r["platform"], r["kind"], r["stage"])][r["cls"]] = r
    out = []
    for c in cells:
        e = dict(tumor=c["tumor"], platform=c["platform"], dataset=f"{c['tumor']} {pname(c['tumor'], c['platform'])}",
                 kind=c["kind"], lin_dir=str(c["lin_dir"]), d46_dir=str(c["d46_dir"]), lin_set=c["lin_set"],
                 d46_set=c["d46_set"], germline_vcf=c["germline_vcf"], d46_tensors_chr1=c["d46_tensors_chr1"],
                 report=c["rep"], stages={})
        for stage in STAGES:
            d = by[(c["tumor"], c["platform"], c["kind"], stage)]
            z = lambda k: {f: d.get(k, {}).get(f, 0) for f in FIELDS}
            add = lambda a, b: {f: a[f] + b[f] for f in FIELDS}
            e["stages"][stage] = dict(
                lin_fp_records=c[f"lin_{stage}_records"], lin_fp_keys=len(c[f"lin_{stage}"]),
                d46_fp_records=c[f"d46_{stage}_records"], d46_fp_keys=len(c[f"d46_{stage}"]),
                truth=c["rep"][f"lin_{stage}"]["truth"],
                germline_in_graph=z("germline|ge46"),
                germline_not_in_graph=add(z("germline|1to45"), z("germline|0")),
                germline_all=z("germline|all"), non_germline=z("other|all"),
                non_germline_in_graph=z("other|ge46"), all=z("all|all"))
        out.append(e)
    return out


# ------------------------------------------------------------------ tables
def p_share(n, d):
    if not d:
        return "-"
    x = 100 * n / d
    if n == 0:
        return "0%"
    return f"{x:.1f}%" if round(x, 1) < 1 else f"{x:.0f}%"


SMALL = 20  # a class with fewer alleles gets "k of n" instead of a percentage and is left out of the ranges


def p_ret(n, d, small=True):
    if not d:
        return "-"
    if small and d < SMALL:
        return f"{n} of {d}"
    x = 100 * n / d
    if n == 0:
        return "0%"
    if x < 0.05:
        return "<0.1%"
    return f"{x:.1f}%" if round(x, 1) < 10 else f"{x:.0f}%"


STAGE_TEXT = {"raw": "before population filtering", "pon": "after population filtering"}


def letter_table(J, kind, stage):
    h = ["Dataset", f"GRCh38-graph FP {STAGE_TEXT[stage]}",
         "Matched-normal germline in HPRC v2.1 graph (share; candidates retained on HPRC v2.1)",
         "Matched-normal germline not in HPRC v2.1 graph (share; candidates retained on HPRC v2.1)",
         "Non-germline FP (% of GRCh38-graph FP)", f"HPRC v2.1 FP {STAGE_TEXT[stage]}"]
    md = ["| " + " | ".join(h) + " |", "|" + "|".join("---" for _ in h) + "|"]
    for e in J:
        if e["kind"] != kind:
            continue
        s = e["stages"][stage]
        n = s["all"]["n"]
        a, b, o = s["germline_in_graph"], s["germline_not_in_graph"], s["non_germline"]
        assert a["n"] + b["n"] + o["n"] == n == s["lin_fp_keys"]
        md.append(f"| {e['dataset']} | {s['lin_fp_records']:,} | {p_share(a['n'], n)}; {p_ret(a['d46_tensor'], a['n'])} | "
                  f"{p_share(b['n'], n)}; {p_ret(b['d46_tensor'], b['n'])} | {p_share(o['n'], n)} | {s['d46_fp_records']:,} |")
    return md


def count_table(J, kind, stage):
    fp_same = "d46_raw_fp" if stage == "raw" else "d46_pon_fp"
    h = ["Dataset", "truth (TP + FN)", "GRCh38-only FP records (allele keys)",
         "germline, >= 46 haplotypes: n; candidate on v2.1; FP on v2.1", "germline, < 46 haplotypes: n; candidate on v2.1; FP on v2.1",
         "not germline: n; candidate on v2.1; FP on v2.1", "HPRC v2.1 FP records (allele keys)",
         "v2.1 FP that are also GRCh38-only FP"]
    md = ["| " + " | ".join(h) + " |", "|" + "|".join("---" for _ in h) + "|"]
    for e in J:
        if e["kind"] != kind:
            continue
        s = e["stages"][stage]
        f = lambda x: f"{x['n']:,}; {x['d46_tensor']:,}; {x[fp_same]:,}"
        md.append(f"| {e['dataset']} | {s['truth']:,} | {s['lin_fp_records']:,} ({s['lin_fp_keys']:,}) | {f(s['germline_in_graph'])} | "
                  f"{f(s['germline_not_in_graph'])} | {f(s['non_germline'])} | {s['d46_fp_records']:,} ({s['d46_fp_keys']:,}) | "
                  f"{s['all'][fp_same]:,} |")
    return md


def rng(vals, fmt, min_d=1):
    """min-max of the (n, d) pairs with d >= min_d, formatted with fmt; the smallest denominator used and the number
    of rows left out."""
    used = [(n, d) for n, d in vals if d >= min_d]
    if not used:
        return f"- ({len(vals)} rows under {min_d} alleles)" if vals else "-"
    lo, hi = min(used, key=lambda v: v[0] / v[1]), max(used, key=lambda v: v[0] / v[1])
    a, b = fmt(*lo), fmt(*hi)
    note = f"smallest n {min(d for _, d in used):,}"
    if len(used) < len(vals):
        note += f"; {len(vals) - len(used)} of {len(vals)} rows left out, class under {min_d} alleles"
    return (a if a == b else f"{a[:-1]}-{b}") + f" ({note})"


def ranges(J):
    md = []
    groups = [("HG008T + COLO829T (six rows)", ("HG008T", "COLO829T")), ("all nine rows", TUMORS),
              ("HCC1395T (three rows)", ("HCC1395T",))]
    for stage in STAGES:
        md += [f"### FP {STAGE_TEXT[stage]}", ""]
        h = ["quantity (sentence it replaces)"] + [g for g, _ in groups]
        md += ["| " + " | ".join(h) + " |", "|" + "|".join("---" for _ in h) + "|"]

        def row(label, kind, cls, what, plats=PLATFORMS):
            cells = []
            for _, tumors in groups:
                vals = []
                for e in J:
                    if e["kind"] != kind or e["tumor"] not in tumors or e["platform"] not in plats:
                        continue
                    s = e["stages"][stage]
                    if what == "share":
                        vals.append((s[cls]["n"], s["all"]["n"]))
                    elif what == "ratio":
                        vals.append((s["d46_fp_records"], s["lin_fp_records"]))
                    else:
                        vals.append((s[cls]["d46_tensor"], s[cls]["n"]))
                cells.append(rng(vals, p_share) if what in ("share", "ratio") else rng(vals, p_ret, SMALL))
            md.append("| " + " | ".join([label] + cells) + " |")

        row('SNV: graph-represented germline, share of GRCh38-graph FP (R3.2 / Table 3.8 caption "71-86%"; R2.5 "over 70%")', "SNV", "germline_in_graph", "share")
        row('SNV: graph-represented germline alleles that keep a candidate on v2.1 ("0.2-0.6%")', "SNV", "germline_in_graph", "ret")
        row('SNV: graph-absent germline alleles that keep a candidate on v2.1 ("83-88%")', "SNV", "germline_not_in_graph", "ret")
        row("SNV: graph-absent germline, share of GRCh38-graph FP", "SNV", "germline_not_in_graph", "share")
        row("SNV: non-germline, share of GRCh38-graph FP", "SNV", "non_germline", "share")
        row('INDEL: graph-represented germline alleles that keep a candidate on v2.1, Illumina ("0.1-0.3%")', "INDEL", "germline_in_graph", "ret", ("Illumina",))
        row('INDEL: graph-represented germline alleles that keep a candidate on v2.1, long reads ("15-22%")', "INDEL", "germline_in_graph", "ret", ("PacBio", "ONT"))
        row("INDEL: graph-represented germline, share of GRCh38-graph FP", "INDEL", "germline_in_graph", "share")
        row("INDEL: graph-absent germline alleles that keep a candidate on v2.1", "INDEL", "germline_not_in_graph", "ret")
        row('INDEL: non-germline, share of GRCh38-graph FP, ONT ("43-57%")', "INDEL", "non_germline", "share", ("ONT",))
        row("INDEL: non-germline, share of GRCh38-graph FP, Illumina + PacBio HiFi / Fiber-seq", "INDEL", "non_germline", "share", ("Illumina", "PacBio"))
        row("SNV: HPRC v2.1 FP as % of GRCh38-graph FP (last column / first count column)", "SNV", None, "ratio")
        row("INDEL: HPRC v2.1 FP as % of GRCh38-graph FP", "INDEL", None, "ratio")
        md.append("")
    return md


def arm_table(J):
    """PASS recall / precision of the two arms at both stages (report.json), so the FP totals can be read with them."""
    h = ["kind", "Dataset", "truth", "GRCh38-only before: TP / FP / recall / precision", "HPRC v2.1 before: TP / FP / recall / precision",
         "GRCh38-only after: TP / FP / recall / precision", "HPRC v2.1 after: TP / FP / recall / precision"]
    md = ["| " + " | ".join(h) + " |", "|" + "|".join("---" for _ in h) + "|"]
    for e in J:
        r = e["report"]
        f = lambda x: f"{x['tp']:,} / {x['fp']:,} / {x['tp'] / x['truth']:.3f} / {x['tp'] / (x['tp'] + x['fp']):.3f}"
        md.append(f"| {e['kind']} | {e['dataset']} | {r['lin_raw']['truth']:,} | {f(r['lin_raw'])} | {f(r['d46_raw'])} | {f(r['lin_pon'])} | {f(r['d46_pon'])} |")
    return md


def write_tables(J, out):
    md = ["# False-positive decomposition, six-tumor current models (ver4, 2026-11-23)", "",
          "chr1, benchmark regions, PASS calls, seed 0. GRCh38-only = six-tumor model without population site features as a model "
          "input (m4m6_graphs_20261110/multicancer6_<P>_<KIND>_linear_nositef); HPRC v2.1 = six-tumor current model "
          "(SNV snv_newpipe_20261109/multicancer6_<P>_SNV_newpipe, INDEL release_s0_<P>_INDEL). Share = % of the GRCh38-graph FP "
          "alleles; candidates retained = % of the alleles of that class with a tensor in the v2.1 test set of the same reads. "
          "In the HPRC v2.1 graph = carried by >= 46 of 463 HPRC v2.1 haplotypes. A class with fewer than 20 alleles shows "
          "\"k of n\" instead of the retained percentage; \"-\" = the class is empty; \"<0.1%\" = above zero and under 0.05 %.", ""]
    n = {("SNV", "raw"): "Response Table 3.8 (= 2.6), SNV, FP before population filtering",
         ("INDEL", "raw"): "Response Table 3.9, INDEL, FP before population filtering",
         ("SNV", "pon"): "Second pair: SNV, FP after population filtering",
         ("INDEL", "pon"): "Second pair: INDEL, FP after population filtering"}
    for (kind, stage), title in n.items():
        md += [f"## {title}", ""] + letter_table(J, kind, stage) + [""]
    md += ["## Counts behind the four tables", "",
           "n = GRCh38-only FP alleles of the class; candidate on v2.1 = of them, alleles with a tensor in the v2.1 set; FP on v2.1 = of "
           "them, alleles that are also a PASS false positive of the v2.1 model at the same stage.", ""]
    for (kind, stage), title in n.items():
        md += [f"### {kind}, FP {STAGE_TEXT[stage]}", ""] + count_table(J, kind, stage) + [""]
    md += ["## PASS recall and precision of the two models (report.json)", "",
           "The two models are not at the same recall: the FP totals of the first and last column go with these recalls.", ""] + arm_table(J) + [""]
    md += ["## Sets and PASS rules", "", "| Dataset | kind | GRCh38-only set | PASS rule (GRCh38-only) | HPRC v2.1 set | PASS rule (v2.1) |", "|---|---|---|---|---|---|"]
    for e in J:
        md.append(f"| {e['dataset']} | {e['kind']} | {e['lin_set']} | {Path(e['lin_dir']).name}: {e['report']['lin_call_rule']} | {e['d46_set']} | "
                  f"{Path(e['d46_dir']).name}: {e['report']['d46_call_rule']} |")
    (out / "tables.md").write_text("\n".join(md) + "\n")
    for (kind, stage), title in n.items():
        (out / f"table_{kind}_{'before' if stage == 'raw' else 'after'}_filter.md").write_text("\n".join(letter_table(J, kind, stage)) + "\n")
    r = ["# Ranges for the R3.2 / R2.5 sentences, recomputed on the ver4 tables", "",
         "Range = min-max over the rows of the group, from the unrounded values, in the cell format of the tables. "
         "\"smallest n\" = the smallest denominator among the rows of the range. Retained percentages of classes with fewer than "
         "20 alleles are left out of a range (stated in the cell).", ""] + ranges(J)
    (out / "ranges.md").write_text("\n".join(r) + "\n")
    return md, r


def october_check(rows):
    """Field-by-field comparison of the raw-stage rows with the old r32_prepon.tsv."""
    old = {}
    with open(OLD_TSV) as f:
        head = f.readline().rstrip("\n").split("\t")
        for line in f:
            d = dict(zip(head, line.rstrip("\n").split("\t")))
            old[(d["tumor"], d["platform"], d["kind"], d["cls"])] = d
    cells = {(r["tumor"], r["platform"], r["kind"]) for r in rows}
    n = bad = 0
    lines = []
    for key, d in old.items():
        if key[:3] not in cells:
            continue
        new = [r for r in rows if r["stage"] == "raw" and (r["tumor"], r["platform"], r["kind"], r["cls"]) == key]
        assert len(new) == 1, key
        for f in head[4:]:
            n += 1
            ok = int(d[f]) == new[0][f]
            bad += not ok
            if not ok or (key[3] in ("germline|ge46", "all|all") and f in ("n", "d46_tensor", "linear_raw_fp", "d46_raw_fp_all")):
                lines.append(f"{' '.join(key)} {f}: old {d[f]} new {new[0][f]} {'OK' if ok else 'DIFFERENT'}")
    new_keys = {(r["tumor"], r["platform"], r["kind"], r["cls"]) for r in rows if r["stage"] == "raw"}
    extra = new_keys - set(old)
    lines.append(f"compared {n} numbers of {len(cells)} cell(s) with {OLD_TSV}: {n - bad} equal, {bad} different; "
                 f"classes only in the new output: {sorted(extra)}")
    return lines, bad + len(extra)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--mode", choices=("ver4", "october"), default="ver4")
    ap.add_argument("--cells", default="HG008T:Illumina:SNV", help="october: TUMOR:PLATFORM:KIND,...")
    ap.add_argument("--out", default=str(HERE))
    ap.add_argument("--tables-only", action="store_true", help="rewrite tables.md / ranges.md from fp_decomp.json")
    a = ap.parse_args()
    out = Path(a.out)
    out.mkdir(parents=True, exist_ok=True)
    if a.tables_only:
        J = json.loads((out / "fp_decomp.json").read_text())
    else:
        cells = ver4_cells() if a.mode == "ver4" else october_cells(a.cells)
        rows = compute(cells)
        with open(out / "fp_decomp.tsv", "w") as f:
            f.write("\t".join(rows[0]) + "\n")
            for r in rows:
                f.write("\t".join(str(v) for v in r.values()) + "\n")
        J = cells_json(cells, rows)
        (out / "fp_decomp.json").write_text(json.dumps(J, indent=1) + "\n")
        if a.mode == "october":
            lines, bad = october_check(rows)
            (out / "october_check.txt").write_text("\n".join(lines) + "\n")
            print("\n".join(lines))
            if bad:
                raise SystemExit("october check FAILED")
    md, r = write_tables(J, out)
    print("\n".join(md))
    print("\n".join(r))


if __name__ == "__main__":
    main()
