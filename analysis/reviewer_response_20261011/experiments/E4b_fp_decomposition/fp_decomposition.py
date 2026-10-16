"""E4b: what the false positives that survive the PoN are (R3.2, R2.5).

For every finished call set (calls_lib.discover), the rtg FP records of the PASS calls after the PoN
(<vcf_dir>/<eval_bed|eval_nobed>/vcfeval/pon_calls/fp.vcf.gz) get flags and one exclusive class, in this order:
  1 germline_same     the normalized allele is carried by the normal's germline truth (any FILTER)
  2 germline_overlap  a germline record of another allele intersects the call's reference interval
  3 near_germline     a germline record intersects the interval widened by 10 bp (labels' near_bp)
  4 hprc              the allele is carried by >= 1 HPRC v1.1 haplotype (AC > 0, normalized)
  5 pop_below_floor   gnomAD or CoLoRSdb carries the allele (normalized) with AF below the PoN's floor for it
                      (SNV 1e-4; INDEL 0.01 in homopolymers >= 7 bp, else 0.001; scripts/filter_panel_of_normals.py)
  6 pop_other_repr    gnomAD / CoLoRSdb carry it at AF >= the floor, but in another representation (the PoN matches
                      POS / REF / ALT exactly)
  7 repeat            HP>=7 or STR (calls_lib.Repeats: the pipeline's RepeatContext)
  8 other
Per-call rows: DATA/E4b/<call set>.<eval>.fp.v2.tsv.gz (recomputed when the rtg file is newer). Tables: tables.md,
summary.json. usage: python fp_decomposition.py
"""
import gzip
import json
import os
import statistics
import sys
import time
from collections import Counter, defaultdict
from pathlib import Path

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parent / "call_level_common"))
import calls_lib as L  # noqa: E402

OUT = L.DATA / "E4b"
CLASSES = ["germline_same", "germline_overlap", "near_germline", "hprc", "pop_below_floor", "pop_other_repr",
           "repeat", "other"]
CLASS_LABEL = {"germline_same": "germline, same allele", "germline_overlap": "germline, overlapping other allele",
               "near_germline": "near germline (<= 10 bp)", "hprc": "HPRC v1.1 allele",
               "pop_below_floor": "gnomAD / CoLoRSdb below PoN floor", "pop_other_repr": "gnomAD / CoLoRSdb >= floor, other representation",
               "repeat": "HP>=7 / STR", "other": "other"}
GERM_CLASSES = ("germline_same", "germline_overlap", "near_germline")
POP_CLASSES = ("hprc", "pop_below_floor", "pop_other_repr")
NEAR = 10
COLS = ["pos", "ref", "alt", "id", "p_somatic", "call_af", "npos", "nref", "nalt", "germline_same_filter",
        "germline_overlap", "germline_near10", "germline_nearest_bp", "hprc_ac", "gnomad_af", "colorsdb_af",
        "pon_floor", "pop_exact", "dbsnp", "dbsnp_somatic_only", "g1000", "repeat", "in_germline_bed", "class"]
EVALS = ("eval_bed", "eval_nobed")


class Ctx:
    """Per-sample resources (germline, repeats, population), loaded once."""

    def __init__(self):
        self.hprc = L.HPRC()
        self.pop = {}
        self.rep = {}
        self.memo = {}  # (germline VCF, pos, ref, alt) -> classification (the same FP recurs across evals and runs)

    def get(self, m):
        f = m["fasta"]
        if f not in self.pop:
            self.pop[f] = L.Population(f)
            self.rep[f] = L.Repeats(f)
        return L.germline(m["germline_vcf"], f), self.pop[f], self.rep[f], L.chrom_seq(f)


def classify_allele(rec, alt, germ, pop, rep, seq, hprc, gbed):
    key = L.normalize(seq, rec.pos, rec.ref, alt)
    if key is None:
        return None
    a, b = L.interval(*key)
    same = germ.keys.get(key)
    overlap = germ.overlapping(a, b) if same is None else 0
    near = germ.count_near(a, b, NEAR) if same is None and not overlap else 0
    nearest = germ.nearest(a, b)
    ac = hprc.lookup(key)
    hits = pop.lookup(rec.pos, rec.ref, alt)
    floor = pop.floor(rec.pos, rec.ref, alt)
    g = hits.get("gnomAD", {}).get("af")
    c = hits.get("CoLoRSdb", {}).get("af")
    afs = [x for x in (g, c) if x is not None]
    pop_hit = "gnomAD" in hits or "CoLoRSdb" in hits
    pop_exact = any(hits.get(k, {}).get("exact") for k in ("gnomAD", "CoLoRSdb"))
    maxaf = max(afs) if afs else None
    repeat = rep.classify(rec.pos, rec.ref, alt)
    if same is not None:
        cls = "germline_same"
    elif overlap:
        cls = "germline_overlap"
    elif near:
        cls = "near_germline"
    elif ac > 0:
        cls = "hprc"
    elif pop_hit and (maxaf is None or maxaf < floor * (1 - 1e-6)):
        cls = "pop_below_floor"
    elif pop_hit:
        cls = "pop_other_repr"
    elif repeat in ("HP>=7", "STR"):
        cls = "repeat"
    else:
        cls = "other"
    dbs = hits.get("dbSNP")
    return dict(npos=key[0], nref=key[1], nalt=key[2], germline_same_filter=same or "", germline_overlap=overlap,
                germline_near10=int(bool(near)), germline_nearest_bp=nearest, hprc_ac=ac, gnomad_af=g,
                colorsdb_af=c, pon_floor=floor, pop_exact=int(bool(pop_exact)), dbsnp=int(dbs is not None),
                dbsnp_somatic_only=int(bool(dbs and dbs.get("somatic_only"))), g1000=int("1000G" in hits),
                repeat=repeat, in_germline_bed=int(gbed.overlaps(a, b)), **{"class": cls})


def process(cs, ev, ctx):
    src = cs.path(ev, "fp")
    out = OUT / f"{cs.id}.{ev}.fp.v2.tsv.gz"
    if out.exists() and out.stat().st_mtime > src.stat().st_mtime:
        return read_rows(out)
    m = L.manifest(cs.set, cs.kind)
    germ, pop, rep, seq = ctx.get(m)
    rows = []
    for rec in L.read_records(src):
        best = None
        for alt in rec.alts:
            mk = (m["germline_vcf"], rec.pos, rec.ref, alt)
            if mk not in ctx.memo:
                ctx.memo[mk] = classify_allele(rec, alt, germ, pop, rep, seq, ctx.hprc, L.bed(m["germline_bed"]))
            d = ctx.memo[mk]
            if d and (best is None or CLASSES.index(d["class"]) < CLASSES.index(best["class"])):
                best = dict(d, alt=alt)
        if best is None:
            continue
        rows.append(dict(best, pos=rec.pos, ref=rec.ref, id=rec.id, p_somatic=rec.info_get("P_SOMATIC"),
                         call_af=L.call_af(rec)))
    L.write_tsv_gz(out, COLS, [[r.get(c) for c in COLS] for r in rows])
    return read_rows(out)


def read_rows(path):
    with gzip.open(path, "rt") as f:
        head = f.readline().rstrip("\n").split("\t")
        return [dict(zip(head, line.rstrip("\n").split("\t"))) for line in f]


def summarize(cs, ev, rows):
    rep = cs.report(ev)
    at = rep["pon"]["at_threshold"]
    cls = Counter(r["class"] for r in rows)
    flags = Counter()
    afs = defaultdict(list)
    for r in rows:
        flags["germline_same"] += r["germline_same_filter"] != ""
        flags["germline_same_pass"] += r["germline_same_filter"] in ("PASS", ".")
        flags["germline_same_hprc"] += r["class"] == "germline_same" and int(r["hprc_ac"] or 0) > 0
        flags["germline_any_hprc"] += r["class"] in GERM_CLASSES and int(r["hprc_ac"] or 0) > 0
        flags["germline_overlap"] += int(r["germline_overlap"] or 0) > 0
        flags["germline_near10"] += int(r["germline_near10"] or 0) > 0
        flags["germline_any_within10"] += (r["germline_same_filter"] != "" or int(r["germline_overlap"] or 0) > 0
                                           or (r["germline_nearest_bp"] not in ("", None) and int(r["germline_nearest_bp"]) <= NEAR))
        flags["hprc"] += int(r["hprc_ac"] or 0) > 0
        flags["gnomad_or_colorsdb"] += r["gnomad_af"] != "" or r["colorsdb_af"] != ""
        flags["dbsnp_or_1000g"] += (r["dbsnp"] == "1" and r["dbsnp_somatic_only"] != "1") or r["g1000"] == "1"
        flags["HP>=7"] += r["repeat"] == "HP>=7"
        flags["STR"] += r["repeat"] == "STR"
        flags["outside_germline_bed"] += r["in_germline_bed"] == "0"
        if r["in_germline_bed"] == "0":
            flags[f"outside_germline_bed:{r['class']}"] += 1
        if r["call_af"] not in ("", "None"):
            afs[r["class"]].append(float(r["call_af"]))
            afs["all"].append(float(r["call_af"]))
    med = {k: statistics.median(v) for k, v in afs.items() if v}
    return dict(id=cs.id, label=cs.label, run=cs.run, vcf_dir=str(cs.vcf_dir), group=cs.group, variant=cs.variant,
                sample=cs.sample, platform=cs.platform, kind=cs.kind, arm=cs.arm, purity=cs.purity, set=cs.set,
                model=cs.model, eval=ev, tp=at["tp_baseline"], tp_call=at["tp_call"], fp=at["fp"], fn=at["fn"],
                precision=at["precision"], recall=at["recall"], f1=at["f1"], fp_rows=len(rows), classes=dict(cls),
                flags=dict(flags), median_call_af=med, unplaced_pass=rep.get("unplaced_pass"))


def class_cells(s):
    n = s["fp_rows"]
    return [L.pct(s["classes"].get(c, 0), n) for c in CLASSES]


def tables(summ):
    out = ["# E4b: false positives that survive the PoN (chr1)", "",
           f"Generated by `fp_decomposition.py` ({time.strftime('%Y-%m-%d %H:%M')}). PASS calls at the run's threshold "
           "after the GRCh38 PoN; rtg FP records. Exclusive classes in the order of the columns (README). "
           "COLO829T's BED is the whole genome, so its eval_bed and eval_nobed rows are equal.", ""]
    head = ["call set", "TP", "FP", "P", "F1"] + [CLASS_LABEL[c] for c in CLASSES]

    def rows_of(pred, ev="eval_bed"):
        return [s for s in summ if s["eval"] == ev and pred(s)]

    sec = 1
    for title, pred in (("In-sample call sets", lambda s: s["group"] == "insample"),
                        ("Joint-model call sets (COLO829T fiberseq P100 / P75 / P50 / P25)", lambda s: s["group"] == "joint"),
                        ("Cross-sample call sets (model of one sample on another sample's tensors)", lambda s: s["group"] == "xs")):
        rs = rows_of(pred)
        out += [f"## {sec}. {title}, in the BED", ""]
        sec += 1
        if not rs:
            out += ["(none finished yet)", ""]
            continue
        out += L.md_table(head, [[s["label"], s["tp"], s["fp"], s["precision"], s["f1"]] + class_cells(s) for s in rs]) + [""]

    # arm comparison
    out += [f"## {sec}. Arms of the same sample / platform / kind (in-sample, main threshold, in the BED)", "",
            "FP counts per class; `lin - d9` and `lin - d46` are the linear arm's extra FPs; `germline share` and "
            "`population share` are the shares of the extra FPs (all classes) in the germline classes (1-3) and the "
            "population classes (4-6). The row in brackets splits class 1 by HPRC v1.1 membership (common germline "
            "alleles that a pangenome carries vs private / rare ones).", ""]
    sec += 1
    groups = defaultdict(dict)
    for s in summ:
        if s["eval"] == "eval_bed" and s["group"] == "insample" and s["variant"] == "main":
            groups[(s["kind"], s["sample"], s["platform"])][s["arm"]] = s
    comp = []
    any_rows = False
    for (kind, sample, platform), arms in sorted(groups.items(), key=lambda x: (L.KIND_ORDER[x[0][0]], x[0][1], L.PLATFORM_ORDER.get(x[0][2], 9))):
        hdr = ["class"] + [a for a in ("linear", "d9", "d46") if a in arms]
        diffs = [g for g in ("d9", "d46") if g in arms and "linear" in arms]
        hdr += [f"lin - {g}" for g in diffs]
        rows = []
        for c in CLASSES + ["germline_same_hprc", "total"]:
            if c == "total":
                vals = {a: s["fp_rows"] for a, s in arms.items()}
            elif c == "germline_same_hprc":
                vals = {a: s["flags"].get("germline_same_hprc", 0) for a, s in arms.items()}
            else:
                vals = {a: s["classes"].get(c, 0) for a, s in arms.items()}
            label = "(germline, same allele: of which HPRC v1.1 allele)" if c == "germline_same_hprc" else CLASS_LABEL.get(c, c)
            rows.append([label] + [vals[a] for a in ("linear", "d9", "d46") if a in arms]
                        + [vals["linear"] - vals[g] for g in diffs])
        share = []
        for g in diffs:
            tot = arms["linear"]["fp_rows"] - arms[g]["fp_rows"]
            ge = sum(arms["linear"]["classes"].get(c, 0) - arms[g]["classes"].get(c, 0) for c in GERM_CLASSES)
            po = sum(arms["linear"]["classes"].get(c, 0) - arms[g]["classes"].get(c, 0) for c in POP_CLASSES)
            gh = arms["linear"]["flags"].get("germline_any_hprc", 0) - arms[g]["flags"].get("germline_any_hprc", 0)
            share.append(dict(vs=g, extra_fp=tot, germline=ge, population=po, germline_hprc=gh,
                              germline_share=ge / tot if tot > 0 else None, population_share=po / tot if tot > 0 else None))
        comp.append(dict(kind=kind, sample=sample, platform=platform, arms=sorted(arms), shares=share,
                         counts={a: dict(s["classes"], total=s["fp_rows"]) for a, s in arms.items()}))
        out += [f"**{sample} {platform} {kind}** (arms: {', '.join(a for a in ('linear', 'd9', 'd46') if a in arms)})", ""]
        out += L.md_table(hdr, rows) + [""]
        if share:
            any_rows = True
            out += ["; ".join(f"linear - {x['vs']}: {x['extra_fp']} extra FP, germline classes {x['germline']} "
                              f"({L.fmt(x['germline_share'])}; of these HPRC v1.1 alleles {x['germline_hprc']}), "
                              f"population classes {x['population']} ({L.fmt(x['population_share'])})" for x in share), ""]
    if not any_rows:
        out += ["No sample / platform / kind has a linear arm and a graph arm finished yet.", ""]

    # flags
    out += [f"## {sec}. Non-exclusive flags and the call's read AF (in the BED)", "",
            "`germline <= 10 bp`: same allele, overlapping, or nearest germline record within 10 bp. `dbSNP / 1000G`: "
            "allele in dbSNP (not somatic-only) or in the 1000G PoN; the PoN does not use them for INDELs. "
            "`outside germline BED`: the call's interval does not overlap the germline truth's BED (dipcall / "
            "DeepVariant confident regions), where a germline allele can be missing from the truth; in brackets those "
            "of class other. Median AF = FORMAT AF of the call (tumor reads).", ""]
    sec += 1
    fh = ["call set", "FP", "germline same (PASS rec.)", "germline same and HPRC v1.1", "germline <= 10 bp", "HPRC v1.1", "gnomAD / CoLoRSdb any AF",
          "dbSNP / 1000G", "HP>=7", "STR", "outside germline BED (of class other)", "median AF all",
          "median AF germline same", "median AF other"]
    rs = [s for s in summ if s["eval"] == "eval_bed"]
    out += L.md_table(fh, [[s["label"], s["fp_rows"],
                            f"{s['flags'].get('germline_same', 0)} ({s['flags'].get('germline_same_pass', 0)})",
                            L.pct(s["flags"].get("germline_same_hprc", 0), s["flags"].get("germline_same", 0)),
                            L.pct(s["flags"].get("germline_any_within10", 0), s["fp_rows"]),
                            L.pct(s["flags"].get("hprc", 0), s["fp_rows"]),
                            L.pct(s["flags"].get("gnomad_or_colorsdb", 0), s["fp_rows"]),
                            L.pct(s["flags"].get("dbsnp_or_1000g", 0), s["fp_rows"]),
                            L.pct(s["flags"].get("HP>=7", 0), s["fp_rows"]), L.pct(s["flags"].get("STR", 0), s["fp_rows"]),
                            f"{L.pct(s['flags'].get('outside_germline_bed', 0), s['fp_rows'])} "
                            f"({s['flags'].get('outside_germline_bed:other', 0)})",
                            s["median_call_af"].get("all"), s["median_call_af"].get("germline_same"),
                            s["median_call_af"].get("other")] for s in rs]) + [""]

    out += [f"## {sec}. Without the BED (eval_nobed), all call sets", ""]
    rs = [s for s in summ if s["eval"] == "eval_nobed"]
    out += L.md_table(head, [[s["label"], s["tp"], s["fp"], s["precision"], s["f1"]] + class_cells(s) for s in rs]) + [""]
    return out, comp


def main():
    t0 = time.time()
    sets = sorted(L.discover(), key=L.sort_key)
    only = os.environ.get("ONLY")  # testing: substring of the call-set id
    if only:
        sets = [cs for cs in sets if only in cs.id]
    L.log(f"{len(sets)} finished call sets")
    ctx = Ctx()
    summ = []
    for cs in sets:
        for ev in EVALS:
            t = time.time()
            rows = process(cs, ev, ctx)
            summ.append(summarize(cs, ev, rows))
            L.log(f"{cs.id} {ev}: {len(rows)} FP rows ({time.time() - t:.0f} s)")
    out, comp = tables(summ)
    (HERE / "tables.md").write_text("\n".join(out) + "\n")
    (HERE / "summary.json").write_text(json.dumps(dict(generated=time.strftime("%Y-%m-%d %H:%M"), call_sets=summ,
                                                       arm_comparison=comp), indent=1) + "\n")
    L.log(f"done in {time.time() - t0:.0f} s")


if __name__ == "__main__":
    main()
