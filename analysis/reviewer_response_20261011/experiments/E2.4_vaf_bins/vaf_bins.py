"""E2.4: call-level recall per truth VAF bin (R2.4 low-VAF detection).

Truth VAF = the tumor read AF of the truth allele in the call set's own tensor set: the `af` of the matched allele in
<set>/tensors/<KIND>/chr1_variant_summary.ndjson, for the tensor that chr1_labels.ndjson links to the truth record
(the representative tensor first; exact allele matches only, not partial ones). A truth allele without a tensor in
its own set takes, in this order: the AF in the same reads' tensor set of another reference arm (linear, d9, d46; same
sample, platform and purity; e.g. a truth the graph absorbed has no graph tensor but a linear one); the truth VCF's own
annotation, COLO829T VAF_Ill (Illumina sets) or VAF_PB (PacBio-derived; fiberseq and ONT sets), multiplied by the purity
for the P75 / P50 / P25 mixtures (an expected VAF); else "unknown" (HG008T's truth VCF has no read VAF). Bins: [0, 0.05), [0.05, 0.1), [0.1, 0.2), [0.2, 0.3), [0.3, 0.5), [0.5, 1].
Recall = rtg TP-baseline / truth of the bin in the BED (pon_calls: PASS after the PoN; raw_calls: before it).
Precision per bin of the call's own FORMAT AF (tumor reads at the call). Share of the truth below the tensor build AF
floor 0.08. Cache: DATA/E2.4/truth_af/<set>.<KIND>.chr1.tsv.gz. usage: python vaf_bins.py
"""
import gzip
import json
import os
import sys
import time
from collections import defaultdict
from pathlib import Path

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parent / "call_level_common"))
import calls_lib as L  # noqa: E402

OUT = L.DATA / "E2.4"
BINS = ((0.0, 0.05, "0-0.05"), (0.05, 0.1, "0.05-0.1"), (0.1, 0.2, "0.1-0.2"), (0.2, 0.3, "0.2-0.3"),
        (0.3, 0.5, "0.3-0.5"), (0.5, 1.01, ">0.5"))
BIN_NAMES = [b[2] for b in BINS] + ["unknown"]
PREFIX = '{"candidate_id": "'


def vaf_bin(x):
    if x is None:
        return "unknown"
    for lo, hi, name in BINS:
        if lo <= x < hi:
            return name
    return ">0.5"


def truth_af(set_name, kind):
    """{(vcf_pos, vcf_ref, vcf_alt): (af, coverage, alt_count, representative)} for chr1 truth alleles with a tensor."""
    d = L.TENSORS / set_name / "tensors" / kind
    lab, summ = d / f"{L.CHROM}_labels.ndjson", d / f"{L.CHROM}_variant_summary.ndjson"
    cache = OUT / "truth_af" / f"{set_name}.{kind}.{L.CHROM}.tsv.gz"
    if cache.exists() and cache.stat().st_mtime > max(lab.stat().st_mtime, summ.stat().st_mtime):
        out = {}
        with gzip.open(cache, "rt") as f:
            next(f)
            for line in f:
                p, r, a, af, cov, alt, rep = line.rstrip("\n").split("\t")
                out[(int(p), r, a)] = (float(af) if af else None, int(cov or 0), int(alt or 0), rep == "1")
        return out
    t = time.time()
    entries = []
    need = set()
    with open(lab) as f:
        for line in f:
            if '"somatic": [{' not in line:
                continue
            r = json.loads(line)
            for e in r["somatic"]:
                entries.append((int(e["vcf_pos"]), e["vcf_ref"], e["vcf_alt"], r["candidate_id"], e["matched_candidate_id"],
                                bool(e.get("representative"))))
                need.add(r["candidate_id"])
    alleles = {}
    n = len(PREFIX)
    with open(summ) as f:
        for line in f:
            cid = line[n:line.index('"', n)]
            if cid not in need:
                continue
            r = json.loads(line)
            per = {a["candidate_id"]: (a.get("af"), a.get("coverage"), a.get("alt_count")) for a in r.get("alleles", [])}
            per.setdefault(r["candidate_id"], (r.get("af"), r.get("coverage"), r.get("alt_count")))
            alleles[cid] = per
    out = {}
    for p, ref, alt, tcid, mcid, rep in entries:
        x = alleles.get(tcid, {}).get(mcid)
        if x is None or x[0] is None:
            continue
        k = (p, ref, alt)
        old = out.get(k)
        cand = (float(x[0]), int(x[1] or 0), int(x[2] or 0), rep)
        if old is None or (rep and not old[3]) or (rep == old[3] and cand[0] > old[0]):
            out[k] = cand
    L.write_tsv_gz(cache, ["vcf_pos", "vcf_ref", "vcf_alt", "af", "coverage", "alt_count", "representative"],
                   [[k[0], k[1], k[2], v[0], v[1], v[2], int(v[3])] for k, v in sorted(out.items())])
    L.log(f"truth AF {set_name} {kind}: {len(entries)} label entries, {len(out)} truth alleles with AF "
          f"({time.time() - t:.0f} s)")
    return out


def vcf_vaf(cs, rec):
    """Fallback: the truth VCF's own VAF (COLO829T only)."""
    if cs.sample != "COLO829T":
        return None
    key = "VAF_Ill" if cs.platform == "Illumina" else "VAF_PB"
    v = rec.info_get(key)
    try:
        x = float(v.split(",")[0])
    except (AttributeError, ValueError):
        return None
    return x * int(cs.purity[1:]) / 100


def other_arm_sets(cs):
    """The same reads' tensor sets of the other reference arms (same sample, platform and purity) that exist."""
    base = f"{cs.sample}_{cs.platform}" + (f"_{cs.purity}" if cs.purity != "P100" else "")
    out = []
    for s in (f"{base}_linear", base, f"{base}_d46"):
        if s != cs.set and (L.TENSORS / s / "tensors" / cs.kind / f"{L.CHROM}_labels.ndjson").exists():
            out.append(s)
    return out


def lookup(afmap, rec):
    af = None
    for alt in rec.alts:
        x = afmap.get((rec.pos, rec.ref, alt))
        if x and (af is None or x[0] > af):
            af = x[0]
    return af


def truth_rows(cs, which, afmap, others):
    """One row per truth record in the BED: found, VAF, its source (tensor / other_arm / truth_vcf / unknown)."""
    rows = []
    for part in ("tp-baseline", "fn"):
        for rec in L.read_records(cs.path("eval_bed", part, which)):
            af, src = lookup(afmap, rec), "tensor"
            if af is None:
                src = "other_arm"
                for m in others:
                    af = lookup(m, rec)
                    if af is not None:
                        break
            if af is None:
                af, src = vcf_vaf(cs, rec), "truth_vcf"
            if af is None:
                src = "unknown"
            rows.append(dict(found=part == "tp-baseline", af=af, src=src, bin=vaf_bin(af)))
    return rows


def calibration(cs, afmap):
    """Truth-VCF VAF against the own-tensor AF over the truth records that have both (COLO829T): n, Pearson r, medians,
    median ratio tensor / VCF, and how many VCF VAF < 0.08 truths do have a tensor."""
    if cs.sample != "COLO829T":
        return None
    import numpy as np
    pairs = []
    for part in ("tp-baseline", "fn"):
        for rec in L.read_records(cs.path("eval_bed", part, "pon_calls")):
            a, v = lookup(afmap, rec), vcf_vaf(cs, rec)
            if a is not None and v is not None:
                pairs.append((a, v))
    if len(pairs) < 3:
        return None
    x = np.array(pairs)
    low = x[x[:, 1] < L.BUILD_AF_FLOOR]
    return dict(n=len(x), pearson_r=float(np.corrcoef(x[:, 0], x[:, 1])[0, 1]), median_tensor_af=float(np.median(x[:, 0])),
                median_vcf_vaf=float(np.median(x[:, 1])), median_ratio=float(np.median(x[:, 0] / np.maximum(x[:, 1], 1e-3))),
                vcf_below_floor_with_tensor=len(low), their_median_tensor_af=float(np.median(low[:, 0])) if len(low) else None,
                vcf_field="VAF_Ill" if cs.platform == "Illumina" else "VAF_PB")


def recall_by_bin(rows):
    d = {b: dict(truth=0, tp=0, tensor=0, other_arm=0, truth_vcf=0, tp_tensor=0) for b in BIN_NAMES}
    for r in rows:
        x = d[r["bin"]]
        x["truth"] += 1
        x["tp"] += r["found"]
        if r["src"] != "unknown":
            x[r["src"]] += 1
        x["tp_tensor"] += r["found"] and r["src"] == "tensor"
    for x in d.values():
        x["recall"] = x["tp"] / x["truth"] if x["truth"] else None
        x["recall_with_tensor"] = x["tp_tensor"] / x["tensor"] if x["tensor"] else None
    return d


def precision_by_bin(cs):
    d = {b: dict(tp=0, fp=0) for b in BIN_NAMES}
    for part, k in (("tp", "tp"), ("fp", "fp")):
        for rec in L.read_records(cs.path("eval_bed", part, "pon_calls")):
            d[vaf_bin(L.call_af(rec))][k] += 1
    for x in d.values():
        x["precision"] = x["tp"] / (x["tp"] + x["fp"]) if x["tp"] + x["fp"] else None
    return d


def floor_share(cs, rows):
    n = len(rows)
    by = {k: [r for r in rows if r["src"] == k] for k in ("tensor", "other_arm", "truth_vcf", "unknown")}
    below = {k: sum(r["af"] < L.BUILD_AF_FLOOR for r in v) for k, v in by.items() if k != "unknown"}
    known_below = sum(below.values())
    return dict(truth=n, with_tensor=len(by["tensor"]), tensor_af_below_floor=below["tensor"],
                without_tensor=n - len(by["tensor"]), other_arm=len(by["other_arm"]), other_arm_below_floor=below["other_arm"],
                truth_vcf=len(by["truth_vcf"]), truth_vcf_below_floor=below["truth_vcf"], unknown=len(by["unknown"]),
                below_floor_known=known_below, share_without_tensor=(n - len(by["tensor"])) / n if n else None,
                share_below_floor_known=known_below / n if n else None)


def cells(d):
    return [("-" if d[b]["truth"] == 0 else f"{L.fmt(d[b]['recall'])} ({d[b]['tp']}/{d[b]['truth']})") for b in BIN_NAMES]


def main():
    t0 = time.time()
    sets = sorted(L.discover(groups=("insample", "joint"), variants=("main", "valbest")), key=L.sort_key)
    only = os.environ.get("ONLY")
    if only:
        sets = [cs for cs in sets if only in cs.id]
    L.log(f"{len(sets)} call sets")
    summ = []
    afmaps = {}
    for cs in sets:
        k = (cs.set, cs.kind)
        if k not in afmaps:
            afmaps[k] = truth_af(cs.set, cs.kind)
        others = []
        for o in other_arm_sets(cs):
            if (o, cs.kind) not in afmaps:
                afmaps[(o, cs.kind)] = truth_af(o, cs.kind)
            others.append(afmaps[(o, cs.kind)])
        pon = truth_rows(cs, "pon_calls", afmaps[k], others)
        raw = truth_rows(cs, "raw_calls", afmaps[k], others)
        src = defaultdict(int)
        for r in pon:
            src[r["src"]] += 1
        summ.append(dict(id=cs.id, label=cs.label, group=cs.group, variant=cs.variant, sample=cs.sample,
                         platform=cs.platform, kind=cs.kind, arm=cs.arm, purity=cs.purity, set=cs.set,
                         recall_pon=recall_by_bin(pon), recall_raw=recall_by_bin(raw), precision_call_af=precision_by_bin(cs),
                         vaf_source=dict(src), floor=floor_share(cs, pon), calibration=calibration(cs, afmaps[k]),
                         overall=dict(truth=len(pon), tp=sum(r["found"] for r in pon))))
        L.log(f"{cs.id}: {len(pon)} truth records")

    out = ["# E2.4: recall per truth VAF bin (chr1)", "",
           f"Generated by `vaf_bins.py` ({time.strftime('%Y-%m-%d %H:%M')}). Truth in the BED; cells `R (TP/truth)`. "
           "VAF: the tensor variant_summary AF of the truth allele in the call set's own tensors; without one, the AF in "
           "the same reads' tensors of another reference arm (linear, d9, d46: same sample, platform, purity); then the "
           "truth VCF's VAF (COLO829T: VAF_Ill for Illumina, VAF_PB otherwise, x purity for mixtures); else unknown "
           "(HG008T). Source counts per call set in the last column.", ""]
    head = ["call set"] + BIN_NAMES + ["all", "VAF source: own tensor / other-arm tensor / truth VCF / unknown"]

    def rrows(rs, key):
        out_rows = []
        for s in rs:
            tp = sum(x["tp"] for x in s[key].values())
            n = sum(x["truth"] for x in s[key].values())
            out_rows.append([s["label"]] + cells(s[key]) + [f"{L.fmt(tp / n if n else None)} ({tp}/{n})",
                            f"{s['vaf_source'].get('tensor', 0)} / {s['vaf_source'].get('other_arm', 0)} / {s['vaf_source'].get('truth_vcf', 0)} / "
                            f"{s['vaf_source'].get('unknown', 0)}"])
        return out_rows

    sec = 1
    for title, pred in (("in-sample", lambda s: s["group"] == "insample"),
                        ("COLO829T fiberseq joint model, per purity", lambda s: s["group"] == "joint")):
        for key, what in (("recall_pon", "PASS after the PoN"), ("recall_raw", "PASS before the PoN")):
            rs = [s for s in summ if pred(s)]
            out += [f"## {sec}. Recall per VAF bin, {what}: {title}", ""]
            sec += 1
            out += (L.md_table(head, rrows(rs, key)) if rs else ["(none finished yet)"]) + [""]
    out += [f"## {sec}. Recall per VAF bin among the truth with an own-set tensor (PASS after the PoN)", "",
            "The model's sensitivity once a tensor exists: cells `R (TP/truth with a tensor)`.", ""]
    sec += 1
    trow = []
    for s in summ:
        d = s["recall_pon"]
        trow.append([s["label"]] + [("-" if d[b]["tensor"] == 0 else
                                     f"{L.fmt(d[b]['recall_with_tensor'])} ({d[b]['tp_tensor']}/{d[b]['tensor']})")
                                    for b in BIN_NAMES[:-1]])
    out += (L.md_table(["call set"] + BIN_NAMES[:-1], trow) if trow else ["(none)"]) + [""]
    # floor
    out += [f"## {sec}. Truth below the tensor build AF floor ({L.BUILD_AF_FLOOR})", "",
            "Per tensor set and kind (the truth in the BED; one row per set). `no tensor`: the truth allele has no exact "
            "tensor in the set (below the build floor, absorbed by the graph, or not a candidate for another reason); "
            "for these the VAF comes from another arm's tensor of the same reads, else the truth VCF, else unknown "
            "(counts below 0.08 in brackets). `known below floor`: own-tensor AF, other-arm AF or truth-VCF VAF < 0.08.", ""]
    sec += 1
    seen = set()
    frows = []
    for s in summ:
        if (s["set"], s["kind"]) in seen:
            continue
        seen.add((s["set"], s["kind"]))
        f = s["floor"]
        frows.append([f"{s['set']} {s['kind']}", f["truth"], f["with_tensor"], f["tensor_af_below_floor"],
                      L.pct(f["without_tensor"], f["truth"]), f"{f['other_arm']} ({f['other_arm_below_floor']})",
                      f"{f['truth_vcf']} ({f['truth_vcf_below_floor']})", f["unknown"],
                      L.pct(f["below_floor_known"], f["truth"])])
    out += (L.md_table(["set", "truth", "with tensor", "tensor AF < 0.08", "no tensor", "no tensor: other-arm tensor AF (< 0.08)",
                        "no tensor: truth-VCF VAF (< 0.08)", "no tensor: VAF unknown", "known below floor"], frows)
            if frows else ["(none)"]) + [""]
    out += [f"## {sec}. Calibration of the truth-VCF fallback (COLO829T): truth-VCF VAF against the own-tensor AF", "",
            "Over the truth records in the BED that have both. VAF_Ill for Illumina sets, VAF_PB otherwise (x purity for "
            "mixtures).", ""]
    sec += 1
    crow = []
    for s in summ:
        c = s.get("calibration")
        if c:
            crow.append([s["label"], c["vcf_field"], c["n"], c["pearson_r"], c["median_tensor_af"], c["median_vcf_vaf"],
                         c["median_ratio"], c["vcf_below_floor_with_tensor"], c["their_median_tensor_af"]])
    out += (L.md_table(["call set", "VCF field", "n", "Pearson r", "median tensor AF", "median VCF VAF",
                        "median ratio tensor / VCF", "VCF VAF < 0.08 with a tensor", "their median tensor AF"], crow)
            if crow else ["(no COLO829T call set)"]) + [""]
    out += [f"## {sec}. Precision per bin of the call's own read AF (FORMAT AF; PASS after the PoN, in the BED)", ""]
    sec += 1
    prow = [[s["label"]] + [("-" if s["precision_call_af"][b]["tp"] + s["precision_call_af"][b]["fp"] == 0 else
                             f"{L.fmt(s['precision_call_af'][b]['precision'])} ({s['precision_call_af'][b]['tp']}/"
                             f"{s['precision_call_af'][b]['tp'] + s['precision_call_af'][b]['fp']})") for b in BIN_NAMES[:-1]]
            for s in summ]
    out += (L.md_table(["call set"] + BIN_NAMES[:-1], prow) if prow else ["(none)"]) + [""]
    (HERE / "tables.md").write_text("\n".join(out) + "\n")
    (HERE / "summary.json").write_text(json.dumps(dict(generated=time.strftime("%Y-%m-%d %H:%M"), call_sets=summ), indent=1) + "\n")
    L.log(f"done in {time.time() - t0:.0f} s")


if __name__ == "__main__":
    main()
