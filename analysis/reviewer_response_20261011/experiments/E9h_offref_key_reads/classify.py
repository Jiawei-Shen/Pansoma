"""E9h step 3: classify every assembly off-reference somatic SNV key by where its A1-read consensus query (step 1)
occurs in the HG008-N v6.2 and HG008-T v3.2 assemblies (step 2, minimap2 -ax sr, one SAM per haplotype), and join the
E1.1 per-key scores (base model, the off-reference rescue, all_H).

Hit rules (a query alignment to one assembly haplotype; primary or secondary):
  earlier  the rule of pansoma_net_v2_runs/offref_negatives_20260929/v2_classify.py (adapted, the original unchanged):
           >= 90 % of the query aligned, NM <= 1 (a mismatch or a 1-bp indel), the site base inside a match block
  strict   the E9d rule (E9c_offref_filter_validation/e9d_asm.py, adapted): >= 90 % aligned, NM <= 1, no insertion or
           deletion anywhere in the alignment, no mismatch within one base of the site (a repeat-length difference is
           never a hit)
Event locus: the alignments of the event's N and T windows (>= 90 % aligned, NM <= 5) +- 100 bp in each haplotype,
plus the event's own position (+- 200 bp) in its HG008-N haplotype.

Classes of a key (per query length and hit rule):
  tumor_only  the query hits HG008-T and not HG008-N: the alt-supporting reads carry a haplotype that only the tumor
              has -> read-supported somatic
  normal      the query hits HG008-N (the patient's germline sequence), except the case below
  ambiguous   no query (the A1 rows give no consensus at the site, or < 25 bp), or the query hits HG008-N only away from
              the event's locus while it hits HG008-T at the event's locus (the 61 / 101 bases occur elsewhere in the
              normal genome: the tensor window cannot tell a somatic allele here from a germline copy elsewhere)
  neither     no hit in either assembly

Controls (build_controls.py): the events' own T (somatic) and N (germline) allele contexts through the same rules.

usage: python classify.py DATA_DIR   (reads DATA_DIR/keys.tsv, windows.tsv, controls.tsv, {N,T}_hap{1,2}[.controls].sam; writes
       results/key_classes.tsv, results/summary.json and tables.md next to this script)"""
import csv
import hashlib
import json
import re
import sys
from collections import Counter, defaultdict
from pathlib import Path

HERE = Path(__file__).resolve().parent
E11 = HERE.parent / "E1.1_ch6_counterfactual" / "results"
E11_DATA = Path("/scratch/jshen/data/reviewer_experiments_20261013/E1.1_ch6_counterfactual")
OLD = Path("/scratch/jshen/data/pansoma_net_v2_runs/offref_negatives_20260929")
E9C = Path("/scratch/jshen/data/reviewer_experiments_20261013/E9c_offref_filter_validation")
TARGETS = ("N_hap1", "N_hap2", "T_hap1", "T_hap2")
HALVES = (30, 50)
RULES = ("earlier", "strict")
CLASSES = ("tumor_only", "normal", "neither", "ambiguous")
LOCUS_PAD, OWN_PAD, WIN_NM = 100, 200, 5
CIG = re.compile(r"(\d+)([MIDNSHP=X])")
CS_OLD = re.compile(r"(:\d+|\*[a-z][a-z]|[+-][a-z]+)")                    # v2_classify.py
CS_E9D = re.compile(r"(:\d+|\*[a-z][a-z]|[+-][a-z]+|~[a-z]{2}\d+[a-z]{2})")  # e9d_asm.py


def site_ok(cs, start, site):
    """v2_classify.py: is query base `site` (0-based, in the SAM record's orientation) inside a match block?"""
    p = start
    for op in CS_OLD.findall(cs):
        if op[0] == ":":
            n = int(op[1:])
            if p <= site < p + n:
                return True
            p += n
        elif op[0] == "*":
            if p == site:
                return False
            p += 1
        elif op[0] == "+":
            if p <= site < p + len(op) - 1:
                return False
            p += len(op) - 1
    return False


def site_exact(cs, start, lo, hi):
    """e9d_asm.py: [lo, hi) inside the alignment with no mismatch / insertion inside and no deletion strictly inside."""
    p = start
    if lo < start:
        return False
    for op in CS_E9D.findall(cs):
        t = op[0]
        if t == ":":
            p += int(op[1:])
        elif t == "*":
            if lo <= p < hi:
                return False
            p += 1
        elif t == "+":
            n = len(op) - 1
            if p < hi and p + n > lo:
                return False
            p += n
        elif t == "-":
            if lo < p < hi:
                return False
        elif t == "~":
            return False
        if p >= hi and t == ":":
            return True
    return p >= hi


def read_sam(path):
    """Mapped records: (qname, reverse, contig, start0, end0, lead clip, clipped bases, NM, cs)."""
    for line in open(path):
        if line[0] == "@":
            continue
        f = line.rstrip("\n").split("\t")
        flag = int(f[1])
        if flag & 4:
            continue
        cig = CIG.findall(f[5])
        lead = int(cig[0][0]) if cig and cig[0][1] in "SH" else 0
        clips = sum(int(n) for n, o in cig if o in "SH")
        rlen = sum(int(n) for n, o in cig if o in "MDN=X")
        tags = {t.split(":", 2)[0]: t.split(":", 2)[2] for t in f[11:]}
        start = int(f[3]) - 1
        yield f[0], bool(flag & 16), f[2], start, start + rlen, lead, clips, int(tags.get("NM", 99)), tags.get("cs", "")


def query_hits(cs, rev, lead, clips, nm, qlen, site):
    """The hit rules met by one alignment of a query (qlen, site) -> (earlier, strict, has_indel)."""
    if qlen - clips < 0.9 * qlen or nm > 1 or not cs:
        return False, False, False
    indel = "+" in cs or "-" in cs
    earlier = site_ok(cs, lead, qlen - 1 - site if rev else site)
    lo, hi = max(0, site - 1), min(qlen, site + 2)
    if rev:
        lo, hi = qlen - hi, qlen - lo
    return earlier, not indel and site_exact(cs, lead, lo, hi), indel


def classify_hits(hs, at_locus):
    """(class, why, N hits, T hits, N at locus, T at locus) of one query from its passing hits under one rule."""
    nh = [x for x in hs if x[0][0] == "N"]
    th = [x for x in hs if x[0][0] == "T"]
    nloc = any(at_locus(x) for x in nh)
    tloc = any(at_locus(x) for x in th)
    if nh:
        cls, why = ("ambiguous", "paralog") if (not nloc and tloc) else ("normal", "at_locus" if nloc else "elsewhere")
    elif th:
        cls, why = "tumor_only", "at_locus" if tloc else "elsewhere"
    else:
        cls, why = "neither", ""
    return cls, why, nh, th, nloc, tloc


def homopolymer_str(w, lo, hi):
    """E9c-style repeat context of the interval [lo, hi) of w: homopolymer >= 7 touching it, else a 2-6 bp unit tandem
    repeat of >= 3 copies overlapping it, else other. Returns (context, homopolymer length)."""
    hp = 0
    i = 0
    while i < len(w):
        j = i
        while j < len(w) and w[j] == w[i]:
            j += 1
        if i <= hi and j >= lo:   # touches [lo, hi) (adjacent runs count: an insertion extends them)
            hp = max(hp, j - i)
        i = j
    if hp >= 7:
        return "homopolymer>=7", hp
    for u in range(2, 7):
        i = 0
        while i + u < len(w):
            j = i + u
            while j < len(w) and w[j] == w[j - u]:
                j += 1
            if (j - i) // u >= 3 and i <= hi and j >= lo:
                return "STR", hp
            i = max(i + 1, j - u + 1)
    return "other", hp


def e9d_classes():
    """The E9d assembly class of the gold keys per set (E9c collect.py asm_class, HG008T branch), read only."""
    def qid(seq, site):
        return "q" + hashlib.sha1(f"{seq}|{site}".encode()).hexdigest()[:16]
    d = E9C / "asm" / "HG008T"
    hits = defaultdict(lambda: defaultdict(int))
    for r in csv.DictReader(open(d / "hits.tsv"), delimiter="\t"):
        hits[r["qid"]][r["role"]] = max(hits[r["qid"]][r["role"]], int(r["n_loci"]))
    mapped = {r["qid"] for r in csv.DictReader(open(d / "queries.tsv"), delimiter="\t")}
    out = {}
    csv.field_size_limit(sys.maxsize)
    for p in sorted((E9C / "gold").glob("*/calls.tsv")):
        name = p.parent.name
        for r in csv.DictReader(open(p), delimiter="\t"):
            alt = "alt_read" if len(r.get("alt_read", "")) >= 40 else "alt_graph"
            c = "not_checked"
            if r.get(alt) and r.get(f"{alt}_site") and qid(r[alt], r[f"{alt}_site"]) in mapped:
                ha = hits.get(qid(r[alt], r[f"{alt}_site"]), {})
                refs = [qid(r[c2], r[f"{c2}_site"]) for c2 in ("ref_graph", "ref_read") if r.get(c2) and r.get(f"{c2}_site")]
                ref_any = any(hits.get(q, {}).get("normal", 0) + hits.get(q, {}).get("tumor", 0) > 0 for q in refs)
                if ha.get("normal", 0) > 0:
                    c = "germline"
                elif ha.get("tumor", 0) > 0:
                    c = "tumor_only"
                else:
                    c = "unsupported" if ref_any else "unverifiable"
            out[(name, r["candidate_id"])] = c
    return out


def old_classes():
    """offref_negatives_20260929 classes (chr1 run and the many-chromosome run) and query sequences, by candidate id."""
    cls, seqs = {}, {}
    for sfx in ("", "_many"):
        q = {r["qid"]: r for r in csv.DictReader(open(OLD / f"queries{sfx}.tsv"), delimiter="\t")}
        name = None
        for line in open(OLD / f"queries{sfx}.fa"):
            if line[0] == ">":
                name = line[1:].strip()
            else:
                seqs[q[name]["candidate_id"]] = line.strip()
        for r in csv.DictReader(open(OLD / f"classes{sfx}.tsv"), delimiter="\t"):
            cls[r["candidate_id"]] = (r["group"], r["class"])
    return cls, seqs


def main():
    data = Path(sys.argv[1])
    keys = list(csv.DictReader(open(data / "keys.tsv"), delimiter="\t"))
    windows = {r["qid"]: r for r in csv.DictReader(open(data / "windows.tsv"), delimiter="\t")}
    qinfo = {}
    for k in keys:
        for h in HALVES:
            if k[f"q{h}_id"]:
                qinfo[k[f"q{h}_id"]] = (int(k[f"q{h}_len"]), int(k[f"q{h}_site"]))
    seqs = {}
    name = None
    for line in open(data / "queries.fa"):
        if line[0] == ">":
            name = line[1:].strip()
        else:
            seqs[name] = line.strip()
    # passing hits per query, rule and target; event loci from the windows
    hits = defaultdict(list)          # (qid, rule) -> [(target, contig, start, end, has_indel, nm)]
    loci = defaultdict(list)          # (event, target) -> [(contig, lo, hi)]
    win_best = defaultdict(lambda: 99)  # (window qid, assembly N/T) -> min NM of a >= 90 %-aligned alignment
    win_full = defaultdict(lambda: 99)  # the same, full length (no clipping), no indel
    for t in TARGETS:
        for q, rev, contig, s, e, lead, clips, nm, cs in read_sam(data / f"{t}.sam"):
            if q in windows:
                L = int(windows[q]["len"])
                if L - clips >= 0.9 * L:
                    win_best[(q, t[0])] = min(win_best[(q, t[0])], nm)
                    if nm <= WIN_NM:
                        loci[(windows[q]["event"], t)].append((contig, s - LOCUS_PAD, e + LOCUS_PAD))
                if clips == 0 and "+" not in cs and "-" not in cs:
                    win_full[(q, t[0])] = min(win_full[(q, t[0])], nm)
                continue
            earlier, strict, indel = query_hits(cs, rev, lead, clips, nm, *qinfo[q])
            if earlier:
                hits[(q, "earlier")].append((t, contig, s, e, indel, nm))
            if strict:
                hits[(q, "strict")].append((t, contig, s, e, indel, nm))
    # control queries (build_controls.py): the events' own T / N allele contexts
    controls = {r["qid"]: r for r in csv.DictReader(open(data / "controls.tsv"), delimiter="\t")}
    chits = defaultdict(list)
    for t in TARGETS:
        for q, rev, contig, s, e, lead, clips, nm, cs in read_sam(data / f"{t}.controls.sam"):
            c = controls[q]
            earlier, strict, indel = query_hits(cs, rev, lead, clips, nm, int(c["qlen"]), int(c["site"]))
            if earlier:
                chits[(q, "earlier")].append((t, contig, s, e, indel, nm))
            if strict:
                chits[(q, "strict")].append((t, contig, s, e, indel, nm))
    for w in windows.values():
        if w["which"] == "N":
            t = "N_" + w["contig"].rsplit("_", 1)[1]
            p = int(w["pos"])
            loci[(w["event"], t)].append((w["contig"], p - OWN_PAD, p + OWN_PAD))

    def at_locus(event, h):
        t, contig, s, e = h[:4]
        return any(c == contig and s < hi and e > lo for c, lo, hi in loci.get((event, t), []))

    own_hap = {w["event"]: "N_" + w["contig"].rsplit("_", 1)[1] for w in windows.values() if w["which"] == "N"}

    def locus_hap(event, nh):
        """Which HG008-N haplotype(s) carry the at-locus hits: own (the event's), other, both."""
        tg = {x[0] for x in nh if at_locus(event, x)}
        return "" if not tg else "both" if len(tg) == 2 else "own" if own_hap[event] in tg else "other"

    # control classes per (event, context, half, rule)
    cclass = {}
    for c in controls.values():
        for rule in RULES:
            cclass[(c["event"], c["which"], int(c["half"]), rule)] = classify_hits(
                chits.get((c["qid"], rule), []), lambda x, e=c["event"]: at_locus(e, x))[0]

    # E1.1 scores and thresholds
    e11 = json.loads((E11 / "assembly_keys.json").read_text())
    runs = json.loads((E11 / "summary.json").read_text())["runs"]
    thr = {run: float(runs[run]["threshold_val_recall_0_9"]) for run in e11}
    scores = {(run, r["key"]): r for run, rows in e11.items() for r in rows}
    e9d = e9d_classes()
    old_cls, old_seqs = old_classes()
    ev_rows = {w["event"]: w for w in windows.values() if w["which"] == "N"}
    nwin = {w["event"]: seqs[q] for q, w in windows.items() if w["which"] == "N"}

    out = []
    for k in keys:
        ev = ev_rows[k["event"]]
        flank, ref, alt = int(ev["flank"]), ev["ref"], ev["alt"]
        ctx, hp = homopolymer_str(nwin[k["event"]], flank, flank + len(ref))
        row = dict(set=k["set"], graph=k["graph"], run=k["run"], key=k["key"], event=k["event"], chrom=k["chrom"],
                   event_contig=ev["contig"], event_pos=ev["pos"], event_allele=f"{ref}>{alt}",
                   event_kind="SNV" if len(ref) == len(alt) == 1 else "INDEL", repeat_context=ctx, homopolymer=hp,
                   af=k["af"], n_a1=k["n_a1"], n_a1_site=k["n_a1_site"], site_alt_frac=k["site_alt_frac"])
        for h in HALVES:
            q = k[f"q{h}_id"]
            row[f"q{h}_status"] = k[f"q{h}_status"]
            row[f"q{h}_len"] = k[f"q{h}_len"] if q else ""
            if h == 30:
                row["q30_seq"] = seqs.get(q, "")
            for rule in RULES:
                tag = f"h{h}_{rule}"
                if not q:
                    row.update({f"{tag}_event_T_context": cclass[(k["event"], "T", h, rule)],
                            f"{tag}_class": "ambiguous", f"{tag}_why": f"no_query:{k[f'q{h}_status']}",
                                f"{tag}_N_hits": "", f"{tag}_T_hits": "", f"{tag}_N_at_locus": "", f"{tag}_N_locus_hap": "",
                                f"{tag}_T_at_locus": "", f"{tag}_N_indel_only": ""})
                    continue
                cls, why, nh, th, nloc, tloc = classify_hits(hits.get((q, rule), []), lambda x: at_locus(k["event"], x))
                row[f"{tag}_event_T_context"] = cclass[(k["event"], "T", h, rule)]
                row.update({f"{tag}_class": cls, f"{tag}_why": why, f"{tag}_N_hits": len(nh), f"{tag}_T_hits": len(th),
                            f"{tag}_N_at_locus": nloc, f"{tag}_N_locus_hap": locus_hap(k["event"], nh),
                            f"{tag}_T_at_locus": tloc, f"{tag}_N_indel_only": bool(nh) and all(x[4] for x in nh)})
        sc = scores.get((k["run"], k["key"])) if k["run"] else None
        if sc:
            t = thr[k["run"]]
            row.update(threshold=f"{t:.4f}", p_base=sc["p_base"], p_site10_high=sc["p_site10_high"], p_all_H=sc["p_all_H"],
                       base_t=sc["p_base"] >= t, base_argmax=sc["argmax_base"] == 1,
                       rescue_argmax=sc["argmax_site10_high"] == 1, rescue_t=sc["p_site10_high"] >= t,
                       allH_t=sc["p_all_H"] >= t, allH_argmax=sc["argmax_all_H"] == 1, e11_label=sc["label"],
                       e11_reason=sc["reason"])
        row["e9d_class"] = e9d.get((k["set"], k["key"]), "")
        oc = old_cls.get(k["key"]) if k["set"] == "HG008T_Illumina" else None
        row["old_class"] = oc[1] if oc else ""
        row["old_query_same"] = (old_seqs.get(k["key"]) == row["q30_seq"]) if oc else ""
        out.append(row)

    res = HERE / "results"
    res.mkdir(exist_ok=True)
    cols = list(dict.fromkeys(c for r in out for c in r))
    with open(res / "key_classes.tsv", "w") as f:
        w = csv.DictWriter(f, fieldnames=cols, delimiter="\t", restval="")
        w.writeheader()
        w.writerows(out)

    # ------------------------------------------------------------------ summaries
    sets = list(dict.fromkeys(r["set"] for r in out))
    S = dict(sets={}, checks={}, windows={})
    for s in sets:
        rs = [r for r in out if r["set"] == s]
        d = dict(n=len(rs), chr1_n=sum(r["chrom"] == "chr1" for r in rs), graph=rs[0]["graph"], run=rs[0]["run"],
                 query_status={h: dict(Counter(r[f"q{h}_status"] for r in rs)) for h in HALVES},
                 n_a1_median=sorted(int(r["n_a1"]) for r in rs)[len(rs) // 2],
                 event_kind=dict(Counter(r["event_kind"] for r in rs)),
                 repeat_context=dict(Counter(r["repeat_context"] for r in rs)), classes={})
        for h in HALVES:
            for rule in RULES:
                tag = f"h{h}_{rule}"
                c = Counter(r[f"{tag}_class"] for r in rs)
                cw = Counter((r[f"{tag}_class"], r[f"{tag}_why"]) for r in rs)
                e = dict(counts={x: c.get(x, 0) for x in CLASSES},
                         chr1={x: sum(r[f"{tag}_class"] == x and r["chrom"] == "chr1" for r in rs) for x in CLASSES},
                         why={f"{a}:{b}": n for (a, b), n in sorted(cw.items())},
                         by_event_kind={kd: {x: sum(r[f"{tag}_class"] == x and r["event_kind"] == kd for r in rs)
                                             for x in CLASSES} for kd in ("SNV", "INDEL")},
                         by_repeat={rc: {x: sum(r[f"{tag}_class"] == x and r["repeat_context"] == rc for r in rs)
                                         for x in CLASSES} for rc in ("homopolymer>=7", "STR", "other")})
                if rs[0]["run"]:
                    calls = {}
                    for x in CLASSES + ("all",):
                        sub = [r for r in rs if x == "all" or r[f"{tag}_class"] == x]
                        calls[x] = dict(n=len(sub), chr1_n=sum(r["chrom"] == "chr1" for r in sub),
                                        **{m: sum(bool(r[m]) for r in sub) for m in
                                           ("base_t", "base_argmax", "rescue_argmax", "rescue_t", "allH_t", "allH_argmax")},
                                        chr1_rescue_argmax=sum(bool(r["rescue_argmax"]) and r["chrom"] == "chr1" for r in sub))
                    e["calls"] = calls
                d["classes"][tag] = e
        d["earlier_normal_indel_only"] = {h: sum(r[f"h{h}_earlier_class"] == "normal" and r[f"h{h}_earlier_N_indel_only"] is True
                                                 for r in rs) for h in HALVES}
        d["transition_h30_earlier_to_strict"] = {f"{a}->{b}": n for (a, b), n in sorted(Counter(
            (r["h30_earlier_class"], r["h30_strict_class"]) for r in rs).items())}
        d["normal_at_locus_hap"] = {f"h{h}_{rule}": dict(Counter(r[f"h{h}_{rule}_N_locus_hap"] for r in rs
                                                                 if r[f"h{h}_{rule}_class"] == "normal" and r[f"h{h}_{rule}_why"] == "at_locus"))
                                    for h in HALVES for rule in RULES}
        d["class_by_event_T_context"] = {f"h{h}_{rule}": {f"{a}|{b}": n for (a, b), n in sorted(Counter(
            (r[f"h{h}_{rule}_class"], r[f"h{h}_{rule}_event_T_context"]) for r in rs).items())} for h in HALVES for rule in RULES}
        d["tumor_only_at_locus"] = {f"h{h}_{rule}": sum(r[f"h{h}_{rule}_class"] == "tumor_only" and r[f"h{h}_{rule}_T_at_locus"] is True
                                                        for r in rs) for h in HALVES for rule in RULES}
        if any(r["e9d_class"] for r in rs):
            d["vs_e9d_h50_strict"] = {f"{a}|{b}": n for (a, b), n in sorted(Counter(
                (r["e9d_class"] or "absent", r["h50_strict_class"]) for r in rs).items())}
            d["vs_e9d_h30_earlier"] = {f"{a}|{b}": n for (a, b), n in sorted(Counter(
                (r["e9d_class"] or "absent", r["h30_earlier_class"]) for r in rs).items())}
        if s == "HG008T_Illumina":
            d["vs_old"] = {f"{a}|{b}": n for (a, b), n in sorted(Counter(
                (r["old_class"] or "absent", r["h30_earlier_class"]) for r in rs).items())}
            d["old_query_same"] = dict(Counter(str(r["old_query_same"]) for r in rs))
        S["sets"][s] = d
    # E1.1 recount against E1.1's own summary (d9: summary.json; d46: report_assembly.json)
    for run in e11:
        rs = [r for r in out if r["run"] == run]
        mine = dict(base=(sum(r["base_t"] for r in rs), sum(r["base_argmax"] for r in rs)),
                    all_H=(sum(r["allH_t"] for r in rs), sum(r["allH_argmax"] for r in rs)),
                    site10_high=(sum(r["rescue_t"] for r in rs), sum(r["rescue_argmax"] for r in rs)))
        res_ = runs[run]["results"]
        if res_["base"].get("assembly"):
            ref = {m: (res_[m]["assembly"]["ge_thr"], res_[m]["assembly"]["argmax_somatic"]) for m in mine}
        else:
            rep = json.loads((E11_DATA / run / "report_assembly.json").read_text())["summary"]
            ref = {m: (rep[m]["ge_thr"], rep[m]["argmax_somatic"]) for m in mine}
        S["checks"][run] = dict(n=len(rs), recount=mine, e11=ref, equal=mine == ref)
    # controls: the events' own allele contexts through the same rules
    kind = {e: "SNV" if len(w["ref"]) == len(w["alt"]) == 1 else "INDEL" for e, w in ev_rows.items()}
    ctrl = {}
    for which in ("T", "N"):
        for h in HALVES:
            for rule in RULES:
                cs_ = [c for c in controls.values() if c["which"] == which and int(c["half"]) == h]
                res_c = {"all": Counter(), "SNV": Counter(), "INDEL": Counter()}
                for c in cs_:
                    cls = cclass[(c["event"], which, h, rule)]
                    res_c["all"][cls] += 1
                    res_c[kind[c["event"]]][cls] += 1
                ctrl[f"{which}_h{h}_{rule}"] = {g: {x: v.get(x, 0) for x in CLASSES} for g, v in res_c.items()}
    S["controls"] = dict(events=len(ev_rows), event_kind=dict(Counter(kind.values())), classes=ctrl)
    # the events' own windows in the assemblies
    evs = sorted({w["event"] for w in windows.values()})
    S["windows"] = dict(
        events=len(evs),
        N_window_in_N_exact=sum(win_full[(f"w_{e}_N", "N")] == 0 for e in evs),
        T_window_in_T_exact=sum(win_full[(f"w_{e}_T", "T")] == 0 for e in evs),
        T_window_in_T_le1=sum(win_full[(f"w_{e}_T", "T")] <= 1 for e in evs),
        T_window_in_N_exact=sum(win_full[(f"w_{e}_T", "N")] == 0 for e in evs),
        N_window_in_T_exact=sum(win_full[(f"w_{e}_N", "T")] == 0 for e in evs),
        T_locus_found=sum(bool(loci.get((e, "T_hap1")) or loci.get((e, "T_hap2"))) for e in evs))
    (res / "summary.json").write_text(json.dumps(S, indent=1) + "\n")
    write_tables(S, out)
    print(json.dumps(S["checks"], indent=1))
    for s, d in S["sets"].items():
        print(s, d["n"], {t: d["classes"][t]["counts"] for t in d["classes"]})


def write_tables(S, out):
    L = ["# E9h tables (generated by classify.py from results/key_classes.tsv; do not edit)", ""]
    sets = S["sets"]

    def cls_cell(c):
        return " / ".join(str(c[x]) for x in CLASSES)

    L += ["## 0. Read-supported somatic keys and the E1.1 calls", "",
          "Read-supported somatic = `tumor_only` (the A1-read consensus occurs in HG008-T v3.2 and not in HG008-N v6.2). "
          "Two tests: `h30 earlier` = the offref_negatives_20260929 method (61-bp query, a 1-bp indel allowed in a hit); "
          "`h50 strict` = the whole 101-column tensor and no indel allowed (the E9d rule). Calls among the tumor-only keys: "
          "base = trained model at the PASS threshold t; rescue = off-reference rescue (`site10_high`, argmax); all_H = "
          "ch6 set to H, at t. `rescue calls: tumor-only / all` = how many of all the rescue's calls on these keys are "
          "read-supported somatic. All autosomes (chr1 in brackets).", "",
          "| set | keys (chr1) | tumor-only, h30 earlier (chr1) | base / rescue / all_H | tumor-only, h50 strict (chr1) | base / rescue / all_H | rescue calls: tumor-only h30 earlier / h50 strict / all |",
          "|---|---|---|---|---|---|---|"]
    for s, d in sets.items():
        a, b = d["classes"]["h30_earlier"], d["classes"]["h50_strict"]
        if "calls" in a:
            ca, cb, call = a["calls"]["tumor_only"], b["calls"]["tumor_only"], a["calls"]["all"]
            L.append(f"| `{s}` | {d['n']} ({d['chr1_n']}) | {ca['n']} ({ca['chr1_n']}) | {ca['base_t']} / {ca['rescue_argmax']} / {ca['allH_t']} | "
                     f"{cb['n']} ({cb['chr1_n']}) | {cb['base_t']} / {cb['rescue_argmax']} / {cb['allH_t']} | "
                     f"{ca['rescue_argmax']} / {cb['rescue_argmax']} / {call['rescue_argmax']} |")
        else:
            L.append(f"| `{s}` (no E1.1 scores) | {d['n']} ({d['chr1_n']}) | {a['counts']['tumor_only']} ({a['chr1']['tumor_only']}) | - | "
                     f"{b['counts']['tumor_only']} ({b['chr1']['tumor_only']}) | - | - |")
    L.append("")

    L += ["## 1. Keys and queries", "",
          "Keys = assembly off-reference somatic SNP keys with a tensor in the set (all autosomes). Query h30 = the A1-read "
          "consensus within +-30 columns of the site (the earlier method), h50 = +-50 (the whole tensor). A1 rows = "
          "sampled rows of the A1 block (median). Event kind = the event's allele in the HG008-N assembly frame; repeat "
          "context of the event in its HG008-N window (homopolymer >= 7 bp touching it, else a 2-6 bp unit repeat of >= 3 "
          "copies, else other).", "",
          "| set | graph | keys (chr1) | A1 rows (median) | h30 query: ok / no query | h50 query: ok / no query | event kind in the assembly frame | repeat context |",
          "|---|---|---|---|---|---|---|---|"]
    for s, d in sets.items():
        qs = {h: d["query_status"][h] for h in HALVES}
        cell = {h: f"{qs[h].get('ok', 0)} / {sum(v for k, v in qs[h].items() if k != 'ok')}"
                   + (f" ({', '.join(f'{k} {v}' for k, v in sorted(qs[h].items()) if k != 'ok')})"
                      if any(k != 'ok' for k in qs[h]) else "") for h in HALVES}
        L.append(f"| `{s}` | {d['graph']} | {d['n']} ({d['chr1_n']}) | {d['n_a1_median']} | {cell[30]} | {cell[50]} | "
                 f"{', '.join(f'{k} {v}' for k, v in sorted(d['event_kind'].items()))} | "
                 f"{', '.join(f'{k} {v}' for k, v in sorted(d['repeat_context'].items()))} |")
    L += ["", "## 2. Classes per set", "",
          "Cells: tumor_only / normal / neither / ambiguous (all autosomes; chr1 in brackets). `earlier` = the "
          "offref_negatives_20260929 hit rule (NM <= 1, a 1-bp indel allowed); `strict` = the E9d rule (no indel, NM <= 1, "
          "no mismatch within one base of the site). h30 earlier = the offref_negatives_20260929 method; h50 strict = the "
          "whole tensor without indel tolerance, the version the README cites (controls in section 5).", "",
          "| set | keys | h30 earlier | h30 strict | h50 earlier | h50 strict |",
          "|---|---|---|---|---|---|"]
    for s, d in sets.items():
        cells = [f"{cls_cell(d['classes'][t]['counts'])} (chr1 {cls_cell(d['classes'][t]['chr1'])})"
                 for t in ("h30_earlier", "h30_strict", "h50_earlier", "h50_strict")]
        L.append(f"| `{s}` | {d['n']} | " + " | ".join(cells) + " |")
    L += ["", "Sub-reasons (class:why; `at_locus` = a hit at the event's locus in that assembly, `elsewhere` = only at "
          "other loci, `paralog` = ambiguous by the locus rule, `no_query:*` = ambiguous without a query):", "",
          "| set | rule | class:why counts |", "|---|---|---|"]
    for s, d in sets.items():
        for t in ("h30_earlier", "h30_strict", "h50_strict"):
            L.append(f"| `{s}` | {t} | {', '.join(f'{k} {v}' for k, v in d['classes'][t]['why'].items())} |")
    L += ["", "## 3. Calls of the E1.1 models by read class", "",
          "Per class: keys (chr1) and how many of them each rule calls: base = the trained model at the PASS threshold t "
          "(validation recall 0.9; `base argmax` in brackets); rescue = `predict`'s off-reference rescue (ch6 raised to "
          "100 at site +- 10, called by argmax: `site10_high` argmax; at t in brackets); all_H = every ch6 cell set to H, "
          "at t (argmax in brackets). E1.1 `results/assembly_keys.json`, thresholds from E1.1 `results/summary.json`.", ""]
    for t in ("h30_earlier", "h30_strict", "h50_strict"):
        L += [f"### {t}", "",
              "| set | class | keys (chr1) | base: t (argmax) | rescue: argmax (t) | all_H: t (argmax) | chr1 rescue argmax |",
              "|---|---|---|---|---|---|---|"]
        for s, d in sets.items():
            calls = d["classes"][t].get("calls")
            if not calls:
                continue
            for x in CLASSES + ("all",):
                c = calls[x]
                L.append(f"| `{s}` | {x} | {c['n']} ({c['chr1_n']}) | {c['base_t']} ({c['base_argmax']}) | "
                         f"{c['rescue_argmax']} ({c['rescue_t']}) | {c['allH_t']} ({c['allH_argmax']}) | {c['chr1_rescue_argmax']} |")
        L.append("")
    L += ["## 4. Why the two hit rules differ", "",
          "`earlier normal, indel-only` = keys that the earlier rule calls normal although every HG008-N hit has an "
          "insertion or deletion (a repeat-length difference). Transitions h30 earlier -> h30 strict. Tumor-only keys "
          "whose HG008-T hit is at the event's locus (`T at locus`).", "",
          "| set | earlier normal, indel-only (h30 / h50) | h30 earlier -> strict | tumor_only with T at locus (h30 earlier / h30 strict / h50 strict) |",
          "|---|---|---|---|"]
    for s, d in sets.items():
        tr = ", ".join(f"{k} {v}" for k, v in d["transition_h30_earlier_to_strict"].items() if k.split("->")[0] != k.split("->")[1])
        ta = d["tumor_only_at_locus"]
        L.append(f"| `{s}` | {d['earlier_normal_indel_only'][30]} / {d['earlier_normal_indel_only'][50]} | {tr or 'none'} | "
                 f"{ta['h30_earlier']} / {ta['h30_strict']} / {ta['h50_strict']} |")
    L += ["", "HG008-N haplotype of the at-locus hits of the `normal:at_locus` keys (own = the haplotype the event is on, "
          "other = the other haplotype, both):", "", "| set | h30 earlier | h30 strict | h50 strict |", "|---|---|---|---|"]
    for s, d in sets.items():
        nh = d["normal_at_locus_hap"]
        L.append(f"| `{s}` | " + " | ".join(", ".join(f"{k} {v}" for k, v in sorted(nh[t].items())) or "none"
                                           for t in ("h30_earlier", "h30_strict", "h50_strict")) + " |")
    L += ["", "Classes by the event's kind in the assembly frame and by repeat context (h30 earlier | h30 strict; "
          "tumor_only / normal / neither / ambiguous):", "",
          "| set | stratum | h30 earlier | h30 strict |", "|---|---|---|---|"]
    for s, d in sets.items():
        for kd in ("SNV", "INDEL"):
            a, b = d["classes"]["h30_earlier"]["by_event_kind"][kd], d["classes"]["h30_strict"]["by_event_kind"][kd]
            if sum(a.values()):
                L.append(f"| `{s}` | event {kd} | {cls_cell(a)} | {cls_cell(b)} |")
        for rc in ("homopolymer>=7", "STR", "other"):
            a, b = d["classes"]["h30_earlier"]["by_repeat"][rc], d["classes"]["h30_strict"]["by_repeat"][rc]
            if sum(a.values()):
                L.append(f"| `{s}` | {rc} | {cls_cell(a)} | {cls_cell(b)} |")
    L += ["", "## 5. Checks", "",
          "E1.1 recount: calls counted here from `assembly_keys.json` and the thresholds, against E1.1's own assembly "
          "summary (d9: `summary.json`, d46: `report_assembly.json`); (at t, argmax).", "",
          "| run | keys | base | all_H | site10_high | equal |", "|---|---|---|---|---|---|"]
    for run, c in S["checks"].items():
        L.append(f"| `{run}` | {c['n']} | {tuple(c['recount']['base'])} / E1.1 {tuple(c['e11']['base'])} | "
                 f"{tuple(c['recount']['all_H'])} / E1.1 {tuple(c['e11']['all_H'])} | "
                 f"{tuple(c['recount']['site10_high'])} / E1.1 {tuple(c['e11']['site10_high'])} | {c['equal']} |")
    il = sets.get("HG008T_Illumina")
    if il:
        L += ["", "Reproduction of offref_negatives_20260929 (HG008T_Illumina d9; old class | h30 earlier class here; "
              "`absent` = key not in the old tables):", "",
              f"- classes: {', '.join(f'{k} {v}' for k, v in il['vs_old'].items())}",
              f"- old query sequence identical to the h30 query here: {', '.join(f'{k} {v}' for k, v in il['old_query_same'].items())}"]
    L += ["", "E9d gold-key assembly class (E9c `gold/*/calls.tsv`, `asm/HG008T/hits.tsv`; 101-bp `alt_read`, strict rule, "
          "minimap2 -N 50) against h50 strict and h30 earlier here (E9d class | class here):", "",
          "| set | vs h50 strict | vs h30 earlier |", "|---|---|---|"]
    for s, d in sets.items():
        if "vs_e9d_h50_strict" in d:
            L.append(f"| `{s}` | {', '.join(f'{k} {v}' for k, v in d['vs_e9d_h50_strict'].items())} | "
                     f"{', '.join(f'{k} {v}' for k, v in d['vs_e9d_h30_earlier'].items())} |")
    L += ["", "Key class against the class of the key's event's own T (somatic) context under the same query length and rule "
          "(key class | event T context class). Event T context `normal`: the somatic allele's own sequence, centred on the "
          "event, also occurs in HG008-N (the key query is centred on the key's site instead, so the two can differ).", "",
          "| set | h30 earlier | h50 strict |", "|---|---|---|"]
    for s, d in sets.items():
        L.append(f"| `{s}` | " + " | ".join(", ".join(f"{k} {v}" for k, v in d["class_by_event_T_context"][t].items())
                                           for t in ("h30_earlier", "h50_strict")) + " |")
    c = S["controls"]
    L += ["", "Controls: the events' own allele contexts (`build_controls.py`; from the 201-bp windows, the allele +- 30 or "
          "+- 50 bases). The T context is somatic by construction (expected tumor_only), the N context is the patient's "
          f"germline haplotype (expected normal). {c['events']} events ({', '.join(f'{k} {v}' for k, v in sorted(c['event_kind'].items()))} "
          "in the assembly frame); cells tumor_only / normal / neither / ambiguous:", "",
          "| context | rule | all events | SNV events | INDEL events |", "|---|---|---|---|---|"]
    for key, v in c["classes"].items():
        which, h, rule = key.split("_")
        L.append(f"| {which} {h} | {rule} | {cls_cell(v['all'])} | {cls_cell(v['SNV'])} | {cls_cell(v['INDEL'])} |")
    w = S["windows"]
    L += ["", "The events' own 201-bp windows (all events with a key here; full-length alignment without indel):", "",
          f"- events {w['events']}; N window exact in HG008-N {w['N_window_in_N_exact']}; T window exact in HG008-T "
          f"{w['T_window_in_T_exact']} (NM <= 1: {w['T_window_in_T_le1']}); T window exact in HG008-N {w['T_window_in_N_exact']}; "
          f"N window exact in HG008-T {w['N_window_in_T_exact']}; event locus found in HG008-T {w['T_locus_found']}", ""]
    (HERE / "tables.md").write_text("\n".join(L) + "\n")


if __name__ == "__main__":
    main()
