"""E1.6: recall funnel of the chr1 somatic truth, per in-sample call set (R1.6 "low sensitivity", R2.3 "explain the
lower recall").

Universe: every truth record that rtg vcfeval evaluated in the benchmark BED on chr1, i.e. tp-baseline + fn of
<run>/vcf_chr1*/eval_bed/vcfeval/pon_calls (PASS truth records of the kind, chr1, in the BED; see
machine_learning/pansoma_net_v2/vcfeval.py). Each record gets the FIRST stage at which it is lost:

  4 called           rtg TP after the PoN (pon_calls tp-baseline).
  3 pon              rtg TP before the PoN (raw_calls tp-baseline) but not after: a PASS call matched it and the
                     population filter (scripts/filter_panel_of_normals.py, as vcf.sh ran it) removed that call.
  2 below_threshold  not matched by any PASS call, and some tensor's representative allele (A1) is exactly the truth
                     allele (its candidate_id is one of the truth allele's node keys in <set>/truth/somatic.graph.tsv).
                     Sub-reasons: below_threshold (no exact tensor is called: p_somatic < the run's threshold, the
                     checkpoint's validation recall-0.9 threshold; SNV off-reference tensors: offref_call), pass_unplaced
                     (called, but linear_vcf found no GRCh38 position), pass_not_matched (called and placed, rtg FN).
  1 no_tensor        no exact representative tensor. Reason from <set>/tensors/somatic.recall.tsv (truth_recall.json
                     per allele): af_floor (filtered, reasons include min_af: tumour AF upper bound < the build AF floor
                     0.08), min_variants (filtered for read support only), absorbed (no_candidate and the arm's graph
                     spells the truth ALT: E5b spell_{d9,d46}; never on the GRCh38-only arm), no_candidate (no candidate,
                     ALT not spelled), other_allele (tensor_non_representative_allele: the site has a tensor whose A1 is
                     another allele), other (complex_allele, no_unique_grch38_node_key). Flag partial_tensor: a tensor
                     labelled somatic for this truth by a partial match (predictions truth_ids, candidate not a key).

Splits: variant kind (SNV / INS / DEL, truth_labels.split_alleles); truth VAF bin [0, 0.08), [0.08, 0.2),
[0.2, 0.4), [0.4, 1] with the source of the VAF recorded per record, in this order:
  tensor          the tumour read AF of the truth allele in the call set's own tensor set (exact match,
                  representative tensor first; E2.4's cache DATA/E2.4/truth_af, vaf_bins.truth_af);
  build_min_af    own set, recall status filtered with min_af: AF upper bound < 0.08 (bin 0-0.08, no value);
  other_arm       the same AF in the same reads' tensor set of another reference arm (vaf_bins.other_arm_sets);
  truth_vcf       COLO829T truth VCF VAF_Ill (Illumina) / VAF_PB (fiberseq, ONT);
  other_platform  HG008T only: median tensor AF over the other platforms' sets (any arm) of the same tumour;
  unknown.
COLO829T also gets a uniform truth-VCF VAF bin for every record (vcf_vaf_bin). Repeat context: calls_lib.Repeats
(indexed_gam_pipeline_v4 truth_labels RepeatContext for INS / DEL; its tract() for SNVs), grouped HP>=7 / STR / other.

Check: per call set, stage 4 = report.json pon.at_threshold.tp_baseline, stages 1-3 = pon fn, stage 3 = raw tp_baseline
- pon tp_baseline, stages 1-2 = raw fn.

Outputs: DATA/E1.6_recall_funnel/per_truth.tsv.gz (one row per call set x truth record); here summary.json and
tables.md. usage: python e16_recall_funnel.py
"""
import gzip
import json
import statistics
import sys
import time
from collections import Counter, defaultdict
from pathlib import Path

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parent / "call_level_common"))
sys.path.insert(0, str(HERE.parent / "E2.4_vaf_bins"))
import calls_lib as L  # noqa: E402
import vaf_bins as V  # noqa: E402
from indexed_gam_pipeline_v4.tensor_postprocessing.truth_labels import split_alleles  # noqa: E402

OUT = Path("/scratch/jshen/data/reviewer_experiments_20261013/E1.6_recall_funnel")
E5B = Path("/scratch/jshen/data/reviewer_experiments_20261013/E5b_absorbed_d46/per_truth.tsv.gz")
STAGES = {1: "no_tensor", 2: "below_threshold", 3: "pon", 4: "called"}
STAGE_LABEL = {1: "no tensor", 2: "tensor, not PASS", 3: "PASS, PoN-tagged", 4: "called"}
REASONS1 = ("af_floor", "min_variants", "absorbed", "no_candidate", "other_allele", "other")
REASONS2 = ("below_threshold", "pass_unplaced", "pass_not_matched")
VAF_BINS = ((0.0, 0.08, "0-0.08"), (0.08, 0.2, "0.08-0.2"), (0.2, 0.4, "0.2-0.4"), (0.4, 1.01, ">0.4"))
VAF_NAMES = [b[2] for b in VAF_BINS] + ["unknown"]
VAF_SOURCES = ("tensor", "build_min_af", "other_arm", "truth_vcf", "other_platform", "unknown")
KIND_NAMES = ("SNV", "INS", "DEL", "COMPLEX")
REPEAT_NAMES = ("HP>=7", "STR", "other")
ARM_NAME = {"linear": "GRCh38-only", "d9": "HPRC v1.1 d9", "d46": "HPRC v2.1 d46"}
STATUS_RANK = {"tensor_representative": 0, "tensor_non_representative_allele": 1, "filtered": 2, "no_candidate": 3,
               "complex_allele": 4, "no_unique_grch38_node_key": 5, "missing": 9}
PLATFORMS = {"HG008T": ("Illumina", "PacBio", "ONT"), "COLO829T": ("Illumina", "fiberseq", "ONT")}


def vaf_bin(x):
    if x is None:
        return "unknown"
    for lo, hi, name in VAF_BINS:
        if lo <= x < hi:
            return name
    return ">0.4"


def call_sets():
    """In-sample main call sets (linear / d9 / d46) and the COLO829T fiberseq joint d9 model at P100 (supplementary)."""
    out = [cs for cs in L.discover(groups=("insample", "joint"), variants=("main",)) if cs.purity == "P100"]
    return sorted(out, key=L.sort_key)


# ------------------------------------------------------------------ per tensor set tables (cached in memory)
_SET = {}


def set_tables(set_name):
    """truth alleles of chr1: keys (vcf_pos, vcf_ref, vcf_alt) -> [truth_id]; truth_id -> dict(keys, kind, status,
    detail)."""
    if set_name in _SET:
        return _SET[set_name]
    by_rec, tid_info = defaultdict(list), {}
    with open(L.TENSORS / set_name / "truth" / "somatic.graph.tsv") as f:
        cols = f.readline().rstrip("\n").split("\t")
        for line in f:
            c = dict(zip(cols, line.rstrip("\n").split("\t")))
            if c["chrom"] != L.CHROM:
                continue
            tid = int(c["truth_id"])
            by_rec[(int(c["vcf_pos"]), c["vcf_ref"], c["vcf_alt"])].append(tid)
            tid_info[tid] = dict(keys=set(c["keys"].split(",")) if c["keys"] else set(), kind=c["kind"],
                                 vcf=(int(c["vcf_pos"]), c["vcf_ref"], c["vcf_alt"]))
    with open(L.TENSORS / set_name / "tensors" / "somatic.recall.tsv") as f:
        cols = f.readline().rstrip("\n").split("\t")
        for line in f:
            c = dict(zip(cols, line.rstrip("\n").split("\t")))
            if c["chrom"] != L.CHROM:
                continue
            t = tid_info[int(c["truth_id"])]
            if t["vcf"][0] != int(c["vcf_pos"]):
                raise SystemExit(f"{set_name}: truth {c['truth_id']} position differs between graph.tsv and recall.tsv")
            t["status"], t["detail"] = c["status"], c["detail"]
    _SET[set_name] = (by_rec, tid_info)
    return _SET[set_name]


_E5B = None


def spelled():
    """(sample, truth_id) -> {d9: 0/1, d46: 0/1, vcf: (pos, ref, alt)} from E5b (chr1 rows)."""
    global _E5B
    if _E5B is None:
        _E5B = {}
        with gzip.open(E5B, "rt") as f:
            cols = f.readline().rstrip("\n").split("\t")
            for line in f:
                c = dict(zip(cols, line.rstrip("\n").split("\t")))
                if c["chrom"] != L.CHROM:
                    continue
                _E5B[(c["sample"], int(c["truth_id"]))] = dict(d9=c["d9_spelled"] == "1", d46=c["d46_spelled"] == "1",
                                                               vcf=(int(c["vcf_pos"]), c["vcf_ref"], c["vcf_alt"]))
    return _E5B


_AF = {}


def truth_af(set_name, kind):
    if (set_name, kind) not in _AF:
        _AF[(set_name, kind)] = V.truth_af(set_name, kind)
    return _AF[(set_name, kind)]


def other_platform_sets(cs):
    if cs.sample != "HG008T":
        return []
    out = []
    for p in PLATFORMS[cs.sample]:
        if p == cs.platform:
            continue
        for suffix in ("_linear", "", "_d46"):
            s = f"{cs.sample}_{p}{suffix}"
            if (L.TENSORS / s / "tensors" / cs.kind / f"{L.CHROM}_labels.ndjson").exists():
                out.append(s)
    return out


# ------------------------------------------------------------------ per call set inputs
def predictions(cs):
    """truth_id -> [dict(cid, p, called, off_reference)] for the tensors labelled somatic (predictions truth_ids:
    representative and partial matches), and the threshold."""
    rep = cs.report("eval_bed")
    thr = rep["threshold"]
    path = L.RUNS / cs.run / "test_chr1" / f"{cs.set}.tensors.{cs.kind}.predictions.ndjson.gz"
    by_tid = defaultdict(list)
    n = excluded = 0
    with gzip.open(path, "rt") as f:
        for line in f:
            n += 1
            if '"truth_ids": []' in line and '"excluded": true' not in line:
                continue
            p = json.loads(line)
            if p.get("excluded"):
                excluded += 1
                continue
            called = p["offref_call"] if "offref_call" in p else p["p_somatic"] >= thr
            for tid in p["truth_ids"]:
                by_tid[tid].append(dict(cid=p["candidate_id"], p=p["p_somatic"], called=bool(called),
                                        off_reference=p.get("off_reference")))
    return by_tid, thr, dict(predictions=n, excluded=excluded, path=str(path))


def call_records(cs):
    """candidate_id -> (placed: bool, pass_before_pon: bool, pon_tag or '') from the PoN-tagged linear VCF and the
    unplaced VCF."""
    out = {}
    pon_vcf = cs.vcf_dir / "eval_bed" / f"{cs.prefix}.linear.pon.vcf.gz"
    for path, placed in ((pon_vcf, True), (cs.unplaced_vcf, False)):
        with gzip.open(path, "rt") as f:
            for line in f:
                if line[0] == "#":
                    continue
                c = line.split("\t", 8)
                filters = set(c[6].split(";"))
                passed = not (filters - {"PASS", "PanelOfNormals"})
                tag = ""
                if "PanelOfNormals" in filters:
                    for item in c[7].split(";"):
                        if item.startswith("PANSOMA_PON="):
                            tag = item[len("PANSOMA_PON="):]
                    tag = tag or "PanelOfNormals"
                for cid in c[2].split(";"):
                    out[cid] = (placed, passed, tag)
    return out


def truth_pon_tags(cs):
    """(pos, ref, alts) -> PoN tag of the truth records (eval_bed/truth.<KIND>.pon.vcf.gz)."""
    out = {}
    for r in L.read_records(cs.vcf_dir / "eval_bed" / f"truth.{cs.kind}.pon.vcf.gz", with_format=False):
        tag = r.info_get("PANSOMA_PON") or ("PanelOfNormals" if "PanelOfNormals" in r.filt else "")
        out[(r.pos, r.ref, tuple(r.alts))] = tag
    return out


def var_kind(rec):
    kinds = []
    for alt in rec.alts:
        for _, _, _, kind, _ in split_alleles(rec.pos, rec.ref, [alt]):
            kinds.append("SNV" if kind == "SNP" else kind)
    kinds = sorted(set(kinds))
    return kinds[0] if len(kinds) == 1 else "COMPLEX"


def vcf_vaf(cs, rec):
    if cs.sample != "COLO829T":
        return None
    v = rec.info_get("VAF_Ill" if cs.platform == "Illumina" else "VAF_PB")
    try:
        return float(v.split(",")[0])
    except (AttributeError, ValueError):
        return None


# ------------------------------------------------------------------ one call set
def process(cs, repeats):
    t0 = time.time()
    rep = cs.report("eval_bed")
    by_rec, tid_info = set_tables(cs.set)
    by_tid, thr, pred_stats = predictions(cs)
    calls = call_records(cs)
    tpon = truth_pon_tags(cs)
    e5b = spelled()
    own_af = truth_af(cs.set, cs.kind)
    arm_afs = [truth_af(s, cs.kind) for s in V.other_arm_sets(cs)]
    plat_afs = [truth_af(s, cs.kind) for s in other_platform_sets(cs)]
    pon = L.rtg_parts(cs, "eval_bed", "pon_calls")
    raw = L.rtg_parts(cs, "eval_bed", "raw_calls")
    key = lambda r: (r.pos, r.ref, tuple(r.alts))  # noqa: E731
    raw_tp = {key(r) for r in raw["tp-baseline"]}
    raw_fn = {key(r) for r in raw["fn"]}
    rows = []
    for part in ("tp-baseline", "fn"):
        for rec in pon[part]:
            k = key(rec)
            if k not in raw_tp and k not in raw_fn:
                raise SystemExit(f"{cs.id}: {k} is in pon_calls but not in raw_calls")
            tids = sorted({t for alt in rec.alts for t in by_rec.get((rec.pos, rec.ref, alt), [])})
            statuses = [(STATUS_RANK.get(tid_info[t].get("status", "missing"), 9), tid_info[t].get("status", "missing"),
                         tid_info[t].get("detail", "")) for t in tids] or [(9, "missing", "")]
            _, status, detail = min(statuses)
            exact, partial = [], []
            for t in tids:
                for x in by_tid.get(t, []):
                    (exact if x["cid"] in tid_info[t]["keys"] else partial).append(x)
            called_exact = [x for x in exact if x["called"]]
            placed = [x for x in called_exact if calls.get(x["cid"], (False,))[0]]
            # stage
            if part == "tp-baseline":
                stage, reason = 4, "called"
            elif k in raw_tp:
                stage, reason = 3, "pon"
            elif exact:
                stage = 2
                reason = ("below_threshold" if not called_exact else "pass_not_matched" if placed else "pass_unplaced")
            else:
                stage = 1
                sp = [e5b.get((cs.sample, t)) for t in tids]
                for t, s in zip(tids, sp):
                    if s is not None and s["vcf"] != tid_info[t]["vcf"]:
                        raise SystemExit(f"{cs.set}: E5b truth {t} differs from somatic.graph.tsv")
                absorbed = cs.arm in ("d9", "d46") and any(s and s[cs.arm] for s in sp)
                if status == "filtered":
                    reason = "af_floor" if "min_af" in detail else "min_variants"
                elif status == "no_candidate":
                    reason = "absorbed" if absorbed else "no_candidate"
                elif status == "tensor_non_representative_allele":
                    reason = "other_allele"
                elif status == "tensor_representative":
                    reason = "inconsistent"  # a representative tensor in the labels but none in the predictions
                else:
                    reason = "other"
            # PoN tags
            call_tag = ""
            for x in called_exact:
                c = calls.get(x["cid"])
                if c and c[2]:
                    call_tag = c[2]
            # VAF
            af, src = None, "unknown"
            a = V.lookup(own_af, rec)
            if a is not None:
                af, src = a, "tensor"
            elif status == "filtered" and "min_af" in detail:
                src = "build_min_af"
            else:
                for m in arm_afs:
                    a = V.lookup(m, rec)
                    if a is not None:
                        af, src = a, "other_arm"
                        break
                if src == "unknown":
                    a = vcf_vaf(cs, rec)
                    if a is not None:
                        af, src = a, "truth_vcf"
                if src == "unknown" and plat_afs:
                    vals = [x for x in (V.lookup(m, rec) for m in plat_afs) if x is not None]
                    if vals:
                        af, src = statistics.median(vals), "other_platform"
            vbin = "0-0.08" if src == "build_min_af" else vaf_bin(af)
            vv = vcf_vaf(cs, rec)
            rclass = repeats.classify(rec.pos, rec.ref, rec.alts[0])
            rows.append(dict(
                run=cs.run, vcf_dir=cs.vcf_dir.name, group=cs.group, set=cs.set, sample=cs.sample, platform=cs.platform,
                arm=cs.arm, kind=cs.kind, chrom=L.CHROM, pos=rec.pos, ref=rec.ref, alt=",".join(rec.alts),
                truth_ids=",".join(map(str, tids)), var_kind=var_kind(rec), stage=stage, stage_name=STAGES[stage],
                reason=reason, recall_status=status, recall_detail=detail,
                n_exact=len(exact), n_partial=len(partial), partial_pass=sum(x["called"] for x in partial),
                max_p_exact=max((x["p"] for x in exact), default=None), threshold=thr,
                exact_called=len(called_exact), exact_called_placed=len(placed),
                truth_pon=tpon.get(k, ""), call_pon=call_tag, vaf=af, vaf_source=src, vaf_bin=vbin,
                vcf_vaf=vv, vcf_vaf_bin=vaf_bin(vv) if cs.sample == "COLO829T" else "",
                repeat_class=rclass, repeat_group=L.repeat_group(rclass)))
    # check against report.json
    st = Counter(r["stage"] for r in rows)
    p_at, r_at = rep["pon"]["at_threshold"], rep["raw"]["at_threshold"]
    check = dict(pon_tp=p_at["tp_baseline"], pon_fn=p_at["fn"], raw_tp=r_at["tp_baseline"], raw_fn=r_at["fn"],
                 stage4=st[4], stage123=st[1] + st[2] + st[3], stage3=st[3], stage12=st[1] + st[2])
    check["exact"] = (check["stage4"] == check["pon_tp"] and check["stage123"] == check["pon_fn"]
                      and check["stage3"] == check["raw_tp"] - check["pon_tp"] and check["stage12"] == check["raw_fn"])
    L.log(f"{cs.id}: {len(rows)} truth records, stages {dict(sorted(st.items()))}, check {check['exact']} "
          f"({time.time() - t0:.0f} s)")
    return rows, dict(threshold=thr, predictions=pred_stats, check=check,
                      truth_pon_tags=rep.get("truth_pon_tags"), unplaced_pass=rep.get("unplaced_pass"))


# ------------------------------------------------------------------ summaries and tables
def funnel(rows):
    n = len(rows)
    st = Counter(r["stage"] for r in rows)
    return dict(n=n, **{STAGES[s]: st[s] for s in STAGES})


def pct(k, n):
    return "-" if not n else f"{100 * k / n:.1f}"


def cell(k, n):
    return "-" if not n else f"{k} ({100 * k / n:.1f} %)"


def run_label(cs_meta):
    s = f"{cs_meta['sample']} {cs_meta['platform']} {cs_meta['kind']}"
    return s


def arm_label(m):
    a = ARM_NAME[m["arm"]]
    return a + (" (joint model)" if m["group"] == "joint" else "")


def funnel_row(label, arm, rows, extra=()):
    f = funnel(rows)
    n = f["n"]
    return [label, arm, str(n)] + [cell(f[STAGES[s]], n) for s in (1, 2, 3, 4)] + list(extra)


def split_rows(meta, by_run, field, levels, label_fn=lambda x: x):
    out = []
    for rid, m in meta.items():
        rows = by_run[rid]
        for lev in levels:
            sub = [r for r in rows if r[field] == lev]
            if not sub:
                continue
            f = funnel(sub)
            out.append([run_label(m), arm_label(m), label_fn(lev), str(f["n"])]
                       + [pct(f[STAGES[s]], f["n"]) for s in (1, 2, 3, 4)])
    return out


def md(header, rows):
    out = ["| " + " | ".join(header) + " |", "|" + "|".join("---" for _ in header) + "|"]
    out += ["| " + " | ".join(str(x) for x in r) + " |" for r in rows]
    return out


def main():
    t0 = time.time()
    OUT.mkdir(parents=True, exist_ok=True)
    sets = call_sets()
    L.log(f"{len(sets)} call sets")
    repeats_by_fasta = {}
    all_rows, meta, by_run = [], {}, {}
    for cs in sets:
        fasta = L.manifest(cs.set, cs.kind)["fasta"]
        if fasta not in repeats_by_fasta:
            repeats_by_fasta[fasta] = L.Repeats(fasta)
        rows, info = process(cs, repeats_by_fasta[fasta])
        all_rows += rows
        by_run[cs.id] = rows
        meta[cs.id] = dict(run=cs.run, vcf_dir=cs.vcf_dir.name, group=cs.group, set=cs.set, sample=cs.sample,
                           platform=cs.platform, arm=cs.arm, kind=cs.kind, **info)
    cols = list(all_rows[0])
    L.write_tsv_gz(OUT / "per_truth.tsv.gz", cols, [[r[c] for c in cols] for r in all_rows])

    # ---- summary.json
    summary = dict(created=time.strftime("%Y-%m-%d %H:%M"), per_truth=str(OUT / "per_truth.tsv.gz"), runs={})
    for rid, m in meta.items():
        rows = by_run[rid]
        s = dict(meta=m, funnel=funnel(rows),
                 reasons_stage1=dict(Counter(r["reason"] for r in rows if r["stage"] == 1)),
                 reasons_stage2=dict(Counter(r["reason"] for r in rows if r["stage"] == 2)),
                 stage1_partial_tensor=sum(r["n_partial"] > 0 for r in rows if r["stage"] == 1),
                 stage1_partial_pass=sum(r["partial_pass"] > 0 for r in rows if r["stage"] == 1),
                 truth_pon_tagged={STAGES[s]: sum(bool(r["truth_pon"]) for r in rows if r["stage"] == s) for s in STAGES},
                 stage3_call_pon_tagged=sum(bool(r["call_pon"]) for r in rows if r["stage"] == 3),
                 stage4_without_exact_tensor=sum(r["n_exact"] == 0 for r in rows if r["stage"] == 4),
                 stage3_without_exact_tensor=sum(r["n_exact"] == 0 for r in rows if r["stage"] == 3),
                 vaf_sources=dict(Counter(r["vaf_source"] for r in rows)),
                 by_kind={k: funnel([r for r in rows if r["var_kind"] == k]) for k in KIND_NAMES
                          if any(r["var_kind"] == k for r in rows)},
                 by_vaf={b: funnel([r for r in rows if r["vaf_bin"] == b]) for b in VAF_NAMES
                         if any(r["vaf_bin"] == b for r in rows)},
                 by_repeat={g: funnel([r for r in rows if r["repeat_group"] == g]) for g in REPEAT_NAMES
                            if any(r["repeat_group"] == g for r in rows)})
        if m["sample"] == "COLO829T":
            s["by_vcf_vaf"] = {b: funnel([r for r in rows if r["vcf_vaf_bin"] == b]) for b in VAF_NAMES
                               if any(r["vcf_vaf_bin"] == b for r in rows)}
        untagged = [r for r in rows if not r["truth_pon"]]
        s["untagged_funnel"] = funnel(untagged)
        summary["runs"][rid] = s
    # pooled per sample x kind x arm (in-sample runs; sum over the three platforms)
    pooled = defaultdict(list)
    for rid, m in meta.items():
        if m["group"] == "insample":
            pooled[(m["sample"], m["kind"], m["arm"])] += by_run[rid]
    def splits(rows, field, levels):
        return {x: funnel([r for r in rows if r[field] == x]) for x in levels if any(r[field] == x for r in rows)}

    summary["pooled_platforms"] = {
        f"{s}|{k}|{a}": dict(funnel=funnel(rows), by_kind=splits(rows, "var_kind", KIND_NAMES),
                             by_vaf=splits(rows, "vaf_bin", VAF_NAMES), by_repeat=splits(rows, "repeat_group", REPEAT_NAMES),
                             **({"by_vcf_vaf": splits(rows, "vcf_vaf_bin", VAF_NAMES)} if s == "COLO829T" else {}),
                             reasons_stage1=dict(Counter(r["reason"] for r in rows if r["stage"] == 1)),
                             stage1_partial_pass=sum(r["partial_pass"] > 0 for r in rows if r["stage"] == 1),
                             stage1_truth_pon_tagged=sum(bool(r["truth_pon"]) for r in rows if r["stage"] == 1),
                             truth_pon_tagged=sum(bool(r["truth_pon"]) for r in rows),
                             untagged_funnel=funnel([r for r in rows if not r["truth_pon"]]),
                             vaf_sources=dict(Counter(r["vaf_source"] for r in rows)))
        for (s, k, a), rows in pooled.items()}
    summary["all_checks_exact"] = all(m["check"]["exact"] for m in meta.values())
    (HERE / "summary.json").write_text(json.dumps(summary, indent=1) + "\n")

    # ---- tables.md
    T = [f"# E1.6 recall funnel: tables (written {summary['created']} by e16_recall_funnel.py)", "",
         "chr1, truth records in the benchmark BED (rtg baseline), PASS at each run's validation recall-0.9 threshold; "
         "stage = the first stage at which a truth record is lost. Arms: GRCh38-only = linear, HPRC v1.1 d9, HPRC v2.1 d46. "
         "`(joint model)` = COLO829T fiberseq P100 called by the joint P100+P75+P50+P25 d9 model (supplementary, not "
         "arm-matched). Per-record rows: `" + str(OUT / "per_truth.tsv.gz") + "`.", ""]
    order = list(meta)
    T += ["## 0. Headline: % of the chr1 truth in the BED lost at each stage, per arm", "",
          "Cell: `no tensor / tensor, not PASS / PASS, PoN-tagged / called` in % of the truth records (n in section 1). "
          "`n.l.` = no local single-set run.", ""]
    head0 = {}
    for rid in order:
        m = meta[rid]
        if m["group"] != "insample":
            continue
        f = summary["runs"][rid]["funnel"]
        head0[(m["sample"], m["platform"], m["kind"], m["arm"])] = " / ".join(pct(f[STAGES[x]], f["n"]) for x in (1, 2, 3, 4))
    rows0 = []
    for sample in ("HG008T", "COLO829T"):
        for kind in ("SNV", "INDEL"):
            for plat in PLATFORMS[sample]:
                rows0.append([f"{sample} {plat} {kind}"] + [head0.get((sample, plat, kind, a), "n.l.")
                                                            for a in ("linear", "d9", "d46")])
    T += md(["chr1 truth in BED", "GRCh38-only", "HPRC v1.1 d9", "HPRC v2.1 d46"], rows0)
    T += [""]
    T += ["## 1. Funnel per call set (chr1, in the BED; n (% of truth))", ""]
    rows1 = []
    for rid in order:
        m, rows = meta[rid], by_run[rid]
        c = m["check"]
        rows1.append(funnel_row(run_label(m), arm_label(m), rows,
                                extra=[f"{m['threshold']:.4f}", f"{c['pon_tp']} / {c['pon_fn']}",
                                       "yes" if c["exact"] else "NO"]))
    T += md(["call set (chr1)", "arm", "truth", "1 no tensor", "2 tensor, not PASS", "3 PASS, PoN-tagged", "4 called",
             "threshold", "rtg TP / FN after PoN (report.json)", "stages reproduce rtg"], rows1)
    T += ["", "Check: stage 4 = rtg TP-baseline after the PoN; 1+2+3 = FN after the PoN; 3 = TP before - TP after the PoN; "
          "1+2 = FN before the PoN (report.json `pon` / `raw` `at_threshold`).", ""]

    T += ["## 2. Stage 1 (no exact tensor) by reason (chr1; n)", ""]
    rows2 = []
    for rid in order:
        m, s = meta[rid], summary["runs"][rid]
        r1 = s["reasons_stage1"]
        rows2.append([run_label(m), arm_label(m), str(s["funnel"]["no_tensor"])] + [str(r1.get(x, 0)) for x in REASONS1]
                     + [str(r1.get("inconsistent", 0)), f"{s['stage1_partial_tensor']} ({s['stage1_partial_pass']})"])
    T += md(["call set (chr1)", "arm", "stage 1", "AF floor (min_af)", "min_variants", "absorbed (graph spells ALT)",
             "no candidate", "other allele at site", "other", "inconsistent", "with a partial-match tensor (PASS)"], rows2)
    T += ["", "Reasons from `<set>/tensors/somatic.recall.tsv`; absorbed = no_candidate and the arm's graph spells the ALT "
          "(E5b spell tables; never on GRCh38-only). Partial-match tensor: labelled somatic for this truth by an overlap "
          "match (another allele / representation); its call cannot be an rtg match.", ""]

    T += ["## 3. Stage 2 (exact tensor, no matching PASS call) by reason (chr1; n)", ""]
    rows3 = []
    for rid in order:
        m, s = meta[rid], summary["runs"][rid]
        r2 = s["reasons_stage2"]
        rows3.append([run_label(m), arm_label(m), str(s["funnel"]["below_threshold"])]
                     + [str(r2.get(x, 0)) for x in REASONS2])
    T += md(["call set (chr1)", "arm", "stage 2", "p_somatic below threshold", "PASS, no GRCh38 position",
             "PASS, placed, rtg FN"], rows3)
    T += [""]

    T += ["## 4. PoN tags of the truth records per stage (chr1; tagged / stage n)", ""]
    rows4 = []
    for rid in order:
        m, s = meta[rid], summary["runs"][rid]
        f, tg = s["funnel"], s["truth_pon_tagged"]
        u = s["untagged_funnel"]
        rows4.append([run_label(m), arm_label(m), f"{sum(tg.values())} / {f['n']} ({pct(sum(tg.values()), f['n'])} %)"]
                     + [f"{tg[STAGES[x]]} / {f[STAGES[x]]}" for x in (1, 2, 3, 4)]
                     + [f"{s['stage3_call_pon_tagged']} / {f['pon']}", str(u["n"])]
                     + [pct(u[STAGES[x]], u["n"]) for x in (1, 2, 3, 4)])
    T += md(["call set (chr1)", "arm", "truth tagged", "stage 1", "stage 2", "stage 3", "stage 4",
             "stage 3: exact tensor's call tagged", "untagged truth n", "untagged: 1 %", "2 %", "3 %", "4 %"], rows4)
    T += ["", "Truth tag: the PoN run by vcfeval over the truth records (`eval_bed/truth.<KIND>.pon.vcf.gz`, allele-exact). "
          "The untagged columns give the funnel of the truth the PoN does not tag (the share any caller can still report "
          "after the PoN).", ""]

    T += ["## 5. Funnel by variant kind (chr1; % of the kind's truth)", ""]
    T += md(["call set (chr1)", "arm", "kind", "n", "1 no tensor %", "2 not PASS %", "3 PoN %", "4 called %"],
            split_rows(meta, by_run, "var_kind", KIND_NAMES))
    T += [""]

    T += ["## 6. Funnel by truth VAF bin (chr1; % of the bin's truth)", "",
          "VAF source order: own tensor AF, own build filter min_af (< 0.08), other-arm tensor AF (same reads), "
          "COLO829T truth VCF (VAF_Ill / VAF_PB), HG008T other-platform tensor AF (median), unknown. Section 7 gives "
          "the sources.", ""]
    T += md(["call set (chr1)", "arm", "VAF bin", "n", "1 no tensor %", "2 not PASS %", "3 PoN %", "4 called %"],
            split_rows(meta, by_run, "vaf_bin", VAF_NAMES))
    T += ["", "## 7. VAF sources (chr1; truth records)", ""]
    rows7 = []
    for rid in order:
        m, s = meta[rid], summary["runs"][rid]
        rows7.append([run_label(m), arm_label(m), str(s["funnel"]["n"])]
                     + [str(s["vaf_sources"].get(x, 0)) for x in VAF_SOURCES])
    T += md(["call set (chr1)", "arm", "truth"] + list(VAF_SOURCES), rows7)
    T += ["", "## 8. COLO829T funnel by the truth VCF's VAF for every record (chr1; VAF_Ill for Illumina, VAF_PB for "
          "fiberseq / ONT; % of the bin)", ""]
    T += md(["call set (chr1)", "arm", "truth-VCF VAF bin", "n", "1 no tensor %", "2 not PASS %", "3 PoN %", "4 called %"],
            split_rows({k: v for k, v in meta.items() if v["sample"] == "COLO829T"}, by_run, "vcf_vaf_bin", VAF_NAMES))
    T += ["", "## 9. Funnel by repeat context (chr1; % of the context's truth)", "",
          "HP>=7: homopolymer tract >= 7 bp; STR: unit 2-6, >= 3 units; other: the rest (calls_lib.Repeats: pipeline "
          "RepeatContext for INS / DEL, its tract() for SNVs).", ""]
    T += md(["call set (chr1)", "arm", "repeat", "n", "1 no tensor %", "2 not PASS %", "3 PoN %", "4 called %"],
            split_rows(meta, by_run, "repeat_group", REPEAT_NAMES))
    T += ["", "## 10. Pooled over the three platforms per sample, kind and arm (in-sample runs; sum of per-run records; "
          "% of n)", ""]
    pkeys = sorted(summary["pooled_platforms"], key=lambda k: (k.split("|")[0], L.KIND_ORDER[k.split("|")[1]],
                                                                 L.ARM_ORDER[k.split("|")[2]]))

    def pooled_rows(field, only=None):
        out = []
        for key in pkeys:
            s, k, a = key.split("|")
            if only and s != only:
                continue
            p = summary["pooled_platforms"][key]
            items = [("all", p["funnel"])] if field is None else list(p[field].items())
            for lev, f in items:
                out.append([f"{s} {k}", ARM_NAME[a], lev, f"{f['n']} ({pct(f['n'], p['funnel']['n'])} %)"]
                           + [pct(f[STAGES[x]], f["n"]) for x in (1, 2, 3, 4)])
        return out

    head = ["1 no tensor %", "2 not PASS %", "3 PoN %", "4 called %"]
    T += md(["sample kind (chr1)", "arm", "", "n (3 platforms; % of all)"] + head, pooled_rows(None))
    T += ["", "### 10a. by variant kind", ""]
    T += md(["sample kind (chr1)", "arm", "kind", "n (3 platforms; % of all)"] + head, pooled_rows("by_kind"))
    T += ["", "### 10b. by truth VAF bin (sources as in section 6)", ""]
    T += md(["sample kind (chr1)", "arm", "VAF bin", "n (3 platforms; % of all)"] + head, pooled_rows("by_vaf"))
    T += ["", "### 10c. COLO829T by the truth VCF's VAF (every record)", ""]
    T += md(["sample kind (chr1)", "arm", "truth-VCF VAF bin", "n (3 platforms; % of all)"] + head,
            pooled_rows("by_vcf_vaf", only="COLO829T"))
    T += ["", "### 10d. by repeat context", ""]
    T += md(["sample kind (chr1)", "arm", "repeat", "n (3 platforms; % of all)"] + head, pooled_rows("by_repeat"))
    T += ["", "### 10e. stage 1 reasons (n (% of all truth n)), PoN tags and the funnel of the untagged truth", ""]
    rows10e = []
    for key in pkeys:
        s, k, a = key.split("|")
        p = summary["pooled_platforms"][key]
        r1, u, n1 = p["reasons_stage1"], p["untagged_funnel"], p["funnel"]["no_tensor"]
        nall = p["funnel"]["n"]
        rows10e.append([f"{s} {k}", ARM_NAME[a], f"{n1} ({pct(n1, nall)} %)"]
                       + [f"{r1.get(x, 0)} ({pct(r1.get(x, 0), nall)} %)" for x in REASONS1]
                       + [str(p["stage1_partial_pass"]), f"{p['stage1_truth_pon_tagged']} / {n1}",
                          f"{p['truth_pon_tagged']} / {p['funnel']['n']} ({pct(p['truth_pon_tagged'], p['funnel']['n'])} %)", str(u["n"])]
                       + [pct(u[STAGES[x]], u["n"]) for x in (1, 2, 3, 4)])
    T += md(["sample kind (chr1)", "arm", "stage 1"] + list(REASONS1)
            + ["stage 1 with a PASS partial-match tensor", "stage 1 PoN-tagged", "truth PoN-tagged", "untagged n",
               "untagged: 1 %", "2 %", "3 %", "4 %"], rows10e)
    R = summary["runs"]
    inv = dict(
        stage3_truth_all_tagged=sum(v["truth_pon_tagged"]["pon"] == v["funnel"]["pon"] for v in R.values()),
        stage4_truth_none_tagged=sum(v["truth_pon_tagged"]["called"] == 0 for v in R.values()),
        stage2_only_below_threshold=sum(set(v["reasons_stage2"]) <= {"below_threshold"} for v in R.values()),
        stage1_inconsistent=sum(v["reasons_stage1"].get("inconsistent", 0) for v in R.values()),
        stage34_without_exact_tensor={k: (v["stage3_without_exact_tensor"], v["stage4_without_exact_tensor"])
                                      for k, v in R.items()
                                      if v["stage3_without_exact_tensor"] or v["stage4_without_exact_tensor"]},
        stage3_exact_call_untagged={k: v["funnel"]["pon"] - v["stage3_call_pon_tagged"] for k, v in R.items()
                                    if v["stage3_call_pon_tagged"] != v["funnel"]["pon"]})
    summary["invariants"] = inv
    (HERE / "summary.json").write_text(json.dumps(summary, indent=1) + "\n")
    T += ["", "## 11. Invariants (all call sets)", "",
          f"- Stage-3 truth records all PoN-tagged themselves: {inv['stage3_truth_all_tagged']} of {len(R)} call sets.",
          f"- No stage-4 (called) truth record is PoN-tagged: {inv['stage4_truth_none_tagged']} of {len(R)} call sets.",
          f"- Stage 2 is only `below_threshold` (no PASS exact tensor unplaced or unmatched by rtg): "
          f"{inv['stage2_only_below_threshold']} of {len(R)} call sets.",
          f"- Stage-1 records whose recall status says `tensor_representative` but without an exact tensor in the "
          f"predictions: {inv['stage1_inconsistent']}.",
          f"- Stage 3 / 4 records without an exact tensor (rtg matched another tensor's call), per call set: "
          f"{inv['stage34_without_exact_tensor'] or 'none'}.",
          f"- Stage-3 records whose exact tensor's call is not PoN-tagged (the match came from another call): "
          f"{inv['stage3_exact_call_untagged'] or 'none'}."]
    T += ["", f"All {len(meta)} call sets reproduce rtg exactly: {summary['all_checks_exact']}. Runtime "
          f"{time.time() - t0:.0f} s.", ""]
    (HERE / "tables.md").write_text("\n".join(T))
    L.log(f"done: {len(all_rows)} rows, {len(meta)} call sets, all checks exact: {summary['all_checks_exact']} "
          f"({time.time() - t0:.0f} s)")


if __name__ == "__main__":
    main()
