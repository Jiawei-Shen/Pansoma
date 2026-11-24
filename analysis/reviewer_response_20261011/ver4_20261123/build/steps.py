"""Edit steps of build_ver4.py. Each step is fn(d, B): d = Docx, B = the build module (fill, header, load, log ...).
Table data come from the generated files under ver4_20261123/ (tables/, regions/, fpdecomp/, offref/, sbs96/,
pathcount/, masking/); numbers inside sentences are taken from the same files (FACTS) wherever a file holds them.
Rules kept after the independent check of 2026-11-23: no reference from a Reviewer 1 / 2 section to a Reviewer 3
table; the model input is described as what it is (population allele frequencies and the population-filter tag);
copied rows get their own ids and the keep-with-next pattern of their table."""
import copy
import hashlib
import json
import re
from pathlib import Path

from docx_helpers import cell_text, copy_row, para_text, replace_image, replace_in_paragraph, set_cell_text, w

V4 = Path(__file__).resolve().parent.parent
WP = "http://schemas.openxmlformats.org/drawingml/2006/wordprocessingDrawing"
WP14 = "http://schemas.microsoft.com/office/word/2010/wordprocessingDrawing"
W14 = "http://schemas.microsoft.com/office/word/2010/wordml"
A16 = "http://schemas.microsoft.com/office/drawing/2014/main"
MARK = "{urn:ver4-build}keep"

SF = "population site features as an input (population allele frequencies and the population-filter tag of each candidate)"
PASS_SNV = "PASS at the threshold that recalls 95% of the validation somatic variants after population filtering"
PASS_INDEL = ("PASS at the validation threshold after population filtering (recall 0.8 on Illumina and PacBio HiFi / "
              "Fiber-seq, best F1 on ONT)")


def _tables(B):
    return B.load("tables/tables.json")


def md_table(path, heading=None, index=0):
    """Rows (list of cell strings, bold markers removed) of the index-th Markdown table after `heading` in a file."""
    lines = Path(path).read_text().splitlines()
    start = 0
    if heading is not None:
        start = next(k for k, x in enumerate(lines) if x.startswith(heading)) + 1
    tables, cur = [], []
    for x in lines[start:]:
        if x.startswith("|"):
            cur.append([c.strip().replace("**", "") for c in x.strip().strip("|").split("|")])
        elif cur:
            tables.append(cur)
            cur = []
    if cur:
        tables.append(cur)
    t = tables[index]
    return [r for r in t if not set("".join(r)) <= set("-: ")]


# ------------------------------------------------------------------------------------------------- text helpers
def _norm(t):
    return t.replace("\xa0", " ")       # Word non-breaking spaces match a plain space in the patterns


def one(d, regex):
    rx = re.compile(regex)
    ps = [p for p in d.paragraphs() if rx.search(_norm(para_text(p)))]
    assert len(ps) == 1, (regex, len(ps))
    return ps[0]


def _resub_p(d, B, p, pattern, new, tag):
    text = para_text(p)
    ms = list(re.finditer(pattern, _norm(text)))
    assert len(ms) == 1, (tag, pattern, len(ms))
    m = ms[0]
    old = text[m.start():m.end()]
    if old != new:
        replace_in_paragraph(d, p, old, new, occurrence=[x.start() for x in re.finditer(re.escape(old), text)].index(m.start()))
    B.log(f"text [{tag[:45]}]: /{pattern[:50]}/ -> {new[:70]!r}")


def resub(d, B, para_regex, pattern, new):
    """In the single body paragraph matching para_regex, replace the text matched by the regex `pattern`
    (matched on the visible text, so the exact characters come from the document) by `new`."""
    _resub_p(d, B, one(d, para_regex), pattern, new, para_regex)


def resub_cell(d, B, title, row_label, col, pattern, new):
    tbl = d.table_by_title(title)
    for tr in d.rows(tbl):
        cells = tr.findall(w("tc"))
        if cell_text(cells[0]) == row_label:
            _resub_p(d, B, cells[col].find(w("p")), pattern, new, f"{title} / {row_label}")
            return
    raise KeyError((title, row_label))


def drop_paragraph(d, B, regex):
    p = one(d, regex)
    p.getparent().remove(p)
    B.log(f"paragraph removed [{regex[:50]}]")


def keep_next(d, title):
    """ver3 pattern: every data row except the last has w:keepNext in all its cells; restore it after rows were added."""
    tbl = d.table_by_title(title)
    rows = d.rows(tbl)
    if rows[1].find(".//" + w("keepNext")) is None:
        return
    for i, tr in enumerate(rows[1:], 1):
        last = i == len(rows) - 1
        for p in tr.iter(w("p")):
            ppr = p.find(w("pPr"))
            if ppr is None:
                continue
            kn = ppr.find(w("keepNext"))
            if last and kn is not None:
                ppr.remove(kn)
            elif not last and kn is None:
                el = ppr.makeelement(w("keepNext"), {})
                st = ppr.find(w("pStyle"))
                if st is not None:
                    st.addnext(el)
                else:
                    ppr.insert(0, el)


# ----------------------------------------------------------------------------------------------- housekeeping steps
def mark_existing(d, B):
    """Remember the empty text runs / empty tracked insertions that ver3 already has, so cleanup() removes only new ones."""
    n = 0
    for t in d.doc.iter(w("t")):
        if not (t.text or ""):
            t.set(MARK, "1")
            n += 1
    for ins in d.doc.iter(w("ins")):
        if ins.getparent().tag == w("p") and not "".join(x.text or "" for x in ins.iter(w("t"))) and ins.find(".//" + w("drawing")) is None:
            ins.set(MARK, "1")
            n += 1
    B.log(f"marked {n} empty elements of ver3")


def cleanup(d, B):
    removed_r = removed_i = 0
    for r in list(d.doc.iter(w("r"))):
        kids = [k for k in r if k.tag != w("rPr")]
        if kids and all(k.tag == w("t") and not (k.text or "") and k.get(MARK) is None for k in kids):
            r.getparent().remove(r)
            removed_r += 1
    for ins in list(d.doc.iter(w("ins"))):
        if ins.getparent().tag == w("p") and ins.get(MARK) is None and len(ins) == 0:
            ins.getparent().remove(ins)
            removed_i += 1
    for e in d.doc.iter():
        if e.get(MARK) is not None:
            del e.attrib[MARK]
    for t in ("1.5", "2.5", "3.15", "2.6", "3.8", "3.9", "2.8", "3.1", "3.2", "3.17"):
        keep_next(d, f"Response Table {t}")
    ids = [e.get("{%s}paraId" % W14) for e in d.doc.iter() if e.get("{%s}paraId" % W14)]
    assert len(ids) == len(set(ids)), "duplicate paraId"
    B.log(f"cleanup: {removed_r} empty runs and {removed_i} empty insertion wrappers removed; keep-with-next restored; {len(ids)} paraIds unique")


# ----------------------------------------------------------------------------------------------- tables from reports
def _same_but_rounding(old, new):
    a, b = [float(x) for x in old.replace("ᵃ", "").split(" / ")], [float(x) for x in new.split(" / ")]
    return all(abs(x - y) <= 0.0011 for x, y in zip(a, b))


def nofilter(d, B):
    T = _tables(B)
    for key, titles in (("nofilter_snv", ("Response Table 1.1", "Response Table 2.9", "Response Table 3.3")),
                        ("nofilter_indel", ("Response Table 1.2", "Response Table 2.10", "Response Table 3.4"))):
        rows = [r[1:] for r in T[key]["rows"]]
        for t in titles:
            tbl = d.table_by_title(t)                # comparator columns: same call sets as ver3, only rounding may move
            for tr, r in zip(d.rows(tbl)[1:], rows):
                cells = tr.findall(w("tc"))
                assert _same_but_rounding(cell_text(cells[4]), r[3]) and _same_but_rounding(cell_text(cells[5]), r[4]), (t, r)
            B.fill(d, t, rows, col0=1, resize=False)
    for n in ("1.1", "2.9", "3.3"):
        ref = " (Response Table 3.17)" if n == "3.3" else ""
        resub(d, B, rf"^Response Table {re.escape(n)}\. ",
              r"Pansoma before our population filter \(it has none of its own\): .*?one run\)\. ",
              "Pansoma: the six-tumor models trained and tested on each reference (chr1 of the three tumors and chr20 of "
              "HCC1395T not in training; one run each), PASS calls before our population filter. The HPRC v1.1 and v2.1 "
              f"models take {SF}{ref}, so their values before and after the filter are close; the GRCh38-only models do not. ")
    resub(d, B, r"^Response Table 3\.3\. ", r"After filtering: Response Table 3\.1 \(single-tumor models\)\.",
          "After filtering: Response Table 3.1.")
    for n in ("1.2", "2.10", "3.4"):
        rule = "PASS as in Response Table 3.17" if n == "3.4" else PASS_INDEL
        resub(d, B, rf"^Response Table {re.escape(n)}\. ",
              r"Pansoma: HG008T and COLO829T, single-tumor models with the final INDEL recipe.*?three seeds\. ?",
              "Pansoma: the six-tumor INDEL models trained and tested on each reference (chr1 of every tumor and chr20 "
              f"of HCC1395T not in training); GRCh38-only and HPRC v1.1, one run; HPRC v2.1, mean of three seeds; {rule}. "
              "The HPRC models take population site features as an input (gnomAD and CoLoRSdb allele frequencies and the "
              "population-filter tag of each candidate), so their values without the filter are not free of population "
              "information; the GRCh38-only models do not. ")
    resub(d, B, r"^Response Table 3\.4\. ", r"After filtering: R3\.1 \(point 3, single-tumor models\)\.", "")


def absorbed(d, B):
    rows = [r[1:] for r in _tables(B)["absorbed_v21"]["rows"]]
    for t in ("Response Table 1.7", "Response Table 2.3", "Response Table 3.6"):
        B.fill(d, t, rows, col0=1, resize=False)


def nearby(d, B):
    T = _tables(B)
    snv = T["t27_snv_nearby"]["rows"]
    assert len(snv) == 9
    B.fill(d, "Response Table 2.7", [[None] + r[1:] for r in snv], resize=False)
    ind = T["t28_indel_nearby"]["rows"]
    assert len(ind) == 12
    tbl = d.table_by_title("Response Table 2.8")
    if len(d.rows(tbl)) == 10:
        for src in (9, 6, 3):                       # bottom-up: the HPRC v2.1 row becomes two rows
            copy_row(d, tbl, src)
    assert len(d.rows(tbl)) == 13
    B.fill(d, "Response Table 2.8", [[None] + r[1:] for r in ind], resize=False)
    resub(d, B, r"^Response Table 2\.7\. ",
          r"Single-tumor models trained on chr2–22 of HG008T on each reference, one training seed per model; PASS at the validation recall-0\.9 threshold\.",
          f"Six-tumor models trained on each reference (chr1 of HG008T not in training), one training seed per model; {PASS_SNV}. "
          "The HPRC models take population site features as an input and the GRCh38-only models do not.")
    resub(d, B, r"^Response Table 2\.8\. ", r"calculated before rounding\.",
          "calculated before rounding. Six-tumor INDEL models after population filtering; GRCh38-only and HPRC v1.1, one "
          "training seed; HPRC v2.1, the best of three seeds (highest chr1 F1 on that platform) and the mean of the three "
          f"seeds; {PASS_INDEL}. The HPRC models take population site features as an input and the GRCh38-only models do not.")
    F = T["facts"]["t27_delta"]
    resub(d, B, r"^Benefit near clustered germline variants\.",
          r"On HG008T, the HPRC graphs give the largest SNV F1 gains where at least two germline variants are nearby \(\+0\.107 to \+0\.204\), mainly through higher precision \(Response Table 2\.7\)\. The benefit therefore extends",
          f"On HG008T, the HPRC graphs give larger SNV F1 gains where germline variants are nearby (+{F['one_or_more'].replace('–', ' to +')} "
          f"with one or at least two nearby germline variants, against +{F['none'].replace('–', ' to +')} with none), through "
          "higher precision and higher or equal recall (Response Table 2.7). On this tumor the benefit therefore extends")


def three_refs(d, B):
    T = _tables(B)
    rows = [list(r) for r in T["t31_three_refs"]["rows"]]
    paired = B.load("regions/whole_chr1_paired.json")
    key = {"HPRC v1.1": "dF1_v1.1_minus_GRCh38-only", "HPRC v2.1": "dF1_v2.1_minus_GRCh38-only"}
    names = {("HG008T", "PacBio HiFi"): "HG008T_PacBio", ("COLO829T", "Fiber-seq"): "COLO829T_fiberseq", ("HCC1395T", "PacBio HiFi"): "HCC1395T_PacBio"}

    def um(s):
        return s.replace("-", "−")

    n_ci = 0
    for r in rows:
        if r[2] in key:
            ds = names.get((r[0], r[1]), f"{r[0]}_{r[1]}")
            c = paired[f"SNV|{ds}|after"][key[r[2]]]
            assert um("%+.3f" % c["d"]) == r[6], (r, c["d"])                 # same call sets, same point value
            r[6] = "%s [%s, %s]" % (r[6], um("%.3f" % c["lo"]), um("%.3f" % c["hi"]))
            n_ci += 1
    assert n_ci == 18
    B.header(d, "Response Table 3.1", 6, "ΔF1, HPRC − GRCh38-only [95% interval]")
    B.fill(d, "Response Table 3.1", rows)
    resub(d, B, r"^Response Table 3\.1\. ", r"Per-seed chr1 SNV results of the three-reference comparison\.",
          "chr1 SNV results of the three-reference comparison (six-tumor models).")
    resub(d, B, r"^Precision, recall and F1 in benchmark regions after population filtering \(our four resource VCFs\)",
          r"at validation-selected checkpoints and PASS thresholds \(validation recall 0\.9 before filtering\)\. Each metric gives seeds.*$",
          "at validation-selected checkpoints and PASS thresholds (Response Table 3.17); one model (seed 0) per row. The "
          "final column gives F1 of the row minus F1 of the GRCh38-only model of the same dataset, with the 95% interval "
          "of a paired bootstrap (2,000 resamples of 1-Mb blocks of chr1) for the HPRC models with site features.")
    resub(d, B, r"^Models and training\. Models are those of the three-reference comparison", r"Models are those of the three-reference comparison.*$",
          "Six-tumor models (HCC1395T, HCC1954T, H1437T, H2009T, HG008T and COLO829T) trained separately on each "
          "reference, with chr1 of HG008T, COLO829T and HCC1395T (and chr20 of HCC1395T) held out. The HPRC models take "
          "population site features as an input (population allele frequencies and the population-filter tag of each "
          "candidate); the rows marked “no site features” are HPRC models trained without them, as the GRCh38-only models "
          "are. The HPRC v2.1 models are those used for the comparisons with other callers.")
    drop_paragraph(d, B, r"^Exploratory, not used in any comparison: COLO829T Illumina on HPRC v1\.1")
    resub(d, B, r"^Reproducibility\. Three SNV models retrained", r"Three SNV models retrained with the same seed",
          "Three models of the earlier single-tumor pipeline retrained with the same seed")
    resub(d, B, r"^Reproducibility\. Three models of the earlier", r"of 30 checkpoints on a second cluster",
          "of 30 checkpoints of the earlier single-tumor models on a second cluster")


def seeds(d, B):
    T = _tables(B)
    rows = T["t32_indel_seeds"]["rows"] + T["t32_snv_seeds"]["rows"]
    B.fill(d, "Response Table 3.2", rows)
    resub(d, B, r"^Response Table 3\.2\. ", r"Per-seed chr1 INDEL results of the final models\.",
          "Per-seed chr1 results of the final INDEL models and of the single-tumor COLO829T SNV models.")
    resub(d, B, r"^Precision, recall and F1 in benchmark regions after gnomAD/CoLoRSdb filtering",
          r"Precision, recall and F1 in benchmark regions after gnomAD/CoLoRSdb filtering, at validation-selected checkpoints and PASS thresholds \(validation recall 0\.8 after filtering\)\.",
          "Precision, recall and F1 in benchmark regions after population filtering (INDELs: gnomAD/CoLoRSdb; SNVs: our "
          "four resource VCFs), at validation-selected checkpoints and PASS thresholds (INDELs: validation recall 0.8 "
          "after filtering on Illumina and PacBio HiFi / Fiber-seq, validation best F1 after filtering on ONT; SNVs: "
          "validation recall 0.95 after filtering).")
    resub(d, B, r"^Models and training\. All models use HPRC v2\.1\.", r"All models use HPRC v2\.1\.", "The INDEL models use HPRC v2.1.")
    resub(d, B, r"^Models and training\. The INDEL models use HPRC v2\.1\.", r"AF ≥ 1e-4 elsewhere\.",
          "AF ≥ 1e-4 elsewhere. The last four rows are SNV models of the current pipeline trained on chr2–22 of COLO829T "
          "alone, on HPRC v2.1 or HPRC v1.1.")


def settings(d, B):
    T = "Response Table 3.17"
    tbl = d.table_by_title(T)
    rows = d.rows(tbl)
    names = [cell_text(r.findall(w("tc"))[0]) for r in rows]
    i = names.index("Checkpoint")
    set_cell_text(d, rows[i].findall(w("tc"))[1],
                  "The epoch with the highest validation somatic F1 after population filtering.")
    if "Calling" not in names:
        copy_row(d, tbl, i, texts=[
            "Calling",
            "PASS at a threshold set on the validation data after population filtering: for SNVs, the threshold that "
            "recalls 95% of the validation somatic variants that survive the filter; for INDELs, the threshold that "
            "recalls 80% of them on Illumina and PacBio HiFi / Fiber-seq, and the threshold with the highest validation "
            "F1 on ONT. Calls without a GRCh38 coordinate are reported in graph coordinates; no off-reference rescue is "
            "applied."])
    resub_cell(d, B, T, "Tensors", 1, r"SNV models also take population site features of the full graph as the input\.",
               "SNV and INDEL models on the HPRC graphs also take per-candidate population site features as input "
               "(population allele frequencies, the population-filter tag and the homopolymer length); the GRCh38-only "
               "models do not.")
    resub_cell(d, B, T, "Labels and splits", 1,
               r"The six-tumor models train on chr2–22 of HG008T and COLO829T, chr2–19 and chr21–22 of HCC1395, and all autosomes of the three CASTLE training lines\.",
               "The six-tumor SNV models train on chr2–22 of HG008T and COLO829T, chr2–19 and chr21–22 of HCC1395, and all "
               "autosomes of the three CASTLE training lines; the six-tumor INDEL models hold out chr1 of all six tumors and "
               "chr20 of HCC1395.")
    B.log("Response Table 3.17: Checkpoint, Tensors and Labels rows rewritten, Calling row added")


# ------------------------------------------------------------------------------------------ tables from new analyses
def _two_lines(s):
    """Letter layout of the region table: break before the arrow and before the interval."""
    return s.replace(" → ", "\n→ ").replace(" [", "\n[")


def regions(d, B):
    t = md_table(V4 / "regions/table_before_filter.md")
    rows = t[1:]
    assert len(rows) == 18 and len(rows[0]) == 9, (len(rows), len(rows[0]))
    out = []
    for r in rows:
        tumor, plat = r[0].split(" ", 1)
        out.append([f"{tumor}\n{plat}", r[1]] + [_two_lines(c) for c in r[2:]])
    for n in ("1.5", "2.5", "3.15"):
        P = rf"^Response Table {re.escape(n)}\. SNV calling by region"
        ref = " (Response Table 3.17)" if n == "3.15" else ""
        B.fill(d, f"Response Table {n}", out)
        resub(d, B, P, r"\(no filtered by the four resource VCFs\)", "(not filtered by the four resource VCFs)")
        resub(d, B, P, r"SMaHT extreme regions: PASS calls and precision on each graph",
              "Six-tumor models trained and tested on each reference (chr1 of the three tumors and chr20 of HCC1395T not in "
              f"training); the HPRC models take {SF} and the GRCh38-only models do not{ref}. "
              "SMaHT extreme regions: PASS calls and precision on each graph")
        resub(d, B, P, r"from the earlier INDEL models of the three-reference comparison(?: \(R3\.1\))?\. HP denotes homopolymers\.",
              "from the six-tumor INDEL models on each reference (one run each), before population filtering.")
    resub(d, B, r"^Response: We thank the reviewer for this comment and agree that the original manuscript lacked region-specific benchmarking",
          r"using the same sequencing reads and analysis pipeline across the GRCh38-only and HPRC graphs",
          "using the same sequencing reads and analysis pipeline across the GRCh38-only and HPRC graphs, except that the HPRC "
          "models take population site features as an input")
    resub(d, B, r"^Calling in difficult regions\.", r"in all twelve dataset–graph comparisons", "in all eighteen dataset–graph comparisons")
    resub(d, B, r"^Calling in difficult regions\.",
          r"INDELs showed a different pattern: in >=7bp homopolymers regions, Illumina gained, whereas long reads lost precision \(Response Table 1\.5\)\.",
          "For INDELs in homopolymers of at least 7 bp, F1 before filtering was higher on the HPRC graphs in 17 of the 18 comparisons (Response Table 1.5).")
    resub(d, B, r"^Fewer false calls in difficult regions\.",
          r"On COLO829T, the GRCh38-only graph's excess false calls are mainly from SMaHT extreme regions, where HPRC v2\.1 reduces false positives and improves precision across all three platforms \(Response Table 2\.5\)\.",
          "In SMaHT extreme regions, before population filtering, HPRC v2.1 reduces false-positive calls and improves precision in all three tumors on all three platforms (Response Table 2.5).")


def fpdecomp(d, B):
    snv = md_table(V4 / "fpdecomp/table_SNV_before_filter.md")[1:]
    ind = md_table(V4 / "fpdecomp/table_INDEL_before_filter.md")[1:]
    assert len(snv) == 9 and len(ind) == 9 and len(snv[0]) == 6
    B.fill(d, "Response Table 2.6", snv)
    B.fill(d, "Response Table 3.8", snv)
    B.fill(d, "Response Table 3.9", ind)
    t26 = d.table_by_title("Response Table 2.6")           # one header wording for the three tables
    head = [cell_text(c).replace(" on HPRC v2.1 model)", " on HPRC v2.1)") for c in d.rows(t26)[0].findall(w("tc"))]
    assert "model" not in " ".join(head), head
    for t in ("Response Table 2.6", "Response Table 3.8", "Response Table 3.9"):
        for c, h in enumerate(head):
            B.header(d, t, c, h)
    for n in ("2.6", "3.8"):
        resub(d, B, rf"^Response Table {re.escape(n)}\. SNV false positives",
              r"Graph-represented germline alleles account for 71–86% of GRCh38-graph FP, but only 0\.2–0\.6% retain candidates, compared with 83–88% of graph-absent germline alleles\.",
              "Six-tumor models: the GRCh38-graph models take no population site features as an input and the HPRC v2.1 "
              "models do, so the final column reflects the graph together with these features. Matched-normal germline: "
              "Dipcall of the HG008-N and COLO829BL assemblies, DeepVariant calls for HCC1395BL. Graph-represented germline "
              "alleles account for 71–87% of GRCh38-graph FP, but only 0.1–0.7% retain candidates, compared with 83–90% of "
              "graph-absent germline alleles.")
    resub(d, B, r"^Response Table 3\.9\. INDEL false positives",
          r"Graph-represented germline INDEL candidate retention is 0\.1–0\.3% on Illumina and 15–22% on long reads; non-germline calls account for 43–57% of GRCh38-graph ONT FP\.",
          "Six-tumor INDEL models (HPRC v2.1: seed 0); the GRCh38-graph models take no population site features as an input "
          "and the HPRC v2.1 models do. Graph-represented germline INDEL candidate retention is at most 0.2% on Illumina; on "
          "long reads it is 15–29% for HG008T and COLO829T and 0.8–2.6% for HCC1395T. A class with fewer than 20 alleles is "
          "given as a count.")
    resub(d, B, r"^The false positives of the GRCh38-only Pansoma model represented in the graph\.",
          r"we trained Pansoma using the same pipeline on the GRCh38-only graph \(which contains only the GRCh38 path\)\.",
          "we trained Pansoma on the GRCh38-only graph (which contains only the GRCh38 path), without the population site features.")
    resub(d, B, r"^Graph representation and precision\.", r"accounted for 71–86% of SNV false positives", "accounted for 71–87% of SNV false positives")
    resub(d, B, r"^Graph representation and precision\.", r"only 0\.2–0\.6% of these alleles retained tensors, whereas 83–88% of graph-absent",
          "only 0.1–0.7% of these alleles retained tensors, whereas 83–90% of graph-absent")
    resub(d, B, r"^Graph representation and precision\.", r"tensor retention was 0\.1–0\.3% on Illumina and 15–22% on long reads \(Response Table 3\.9\)\.",
          "tensor retention was at most 0.2% on Illumina; on long reads it was 15–29% for HG008T and COLO829T and 0.8–2.6% for HCC1395T (Response Table 3.9).")


def offref(d, B):
    t10 = [[c.replace("(100.0%)", "(100%)") for c in r] for r in md_table(V4 / "offref/table_3.10.md")[1:]]
    assert len(t10) == 9 and len(t10[0]) == 9
    assert sum(int(r[7]) for r in t10) == 27 and sum(int(r[1]) for r in t10) == 935
    assert sum(int(r[3].split()[0]) for r in t10) == 40
    B.fill(d, "Response Table 3.10", t10, resize=False)
    t11 = md_table(V4 / "offref/table_3.11.md")[1:]
    assert len(t11) == 11
    tbl = d.table_by_title("Response Table 3.11")
    rows = d.rows(tbl)
    assert len(rows) == 12
    for tr, vals in zip(rows[1:10], t11[:9]):
        for tc, v in zip(tr.findall(w("tc")), vals):
            if cell_text(tc) != v:
                set_cell_text(d, tc, v)
    for tr, vals in zip(rows[10:], t11[9:]):        # total rows: label spans two columns
        cells = tr.findall(w("tc"))
        assert len(cells) == 3 and vals[1] == "", (len(cells), vals)
        for tc, v in zip(cells, [vals[0], vals[2], vals[3]]):
            if cell_text(tc) != v:
                set_cell_text(d, tc, v)
    B.log("Response Table 3.11: 9 rows and 2 total rows rewritten")
    resub(d, B, r"^Response\. We thank the reviewer for suggesting these additional validation analyses\.",
          r"\(HPRC v1\.1 and v2\.1; HCC1395T: v2\.1 only; Response Table 3\.10\)", "(HPRC v2.1, six-tumor models; Response Table 3.10)")
    resub(d, B, r"^In the default output, 34 calls from the HPRC v2\.1 model",
          r"34 calls from the HPRC v2\.1 model in chr1 were supported by tumor-only assemblies, consistent with potential somatic origin\.",
          "40 of the 935 calls from the HPRC v2.1 models in chr1 were supported by tumor-only assemblies, consistent with "
          "potential somatic origin; 27 of them also have no ALT read in the matched normal.")
    C = r"^Response Table 3\.10\. Validation of graph-coordinate"
    resub(d, B, C, r"HCC1395T uses six-tumor models\. Of 541 calls, 34 \(6\.3%\) have tumor-only assembly support\.",
          "All nine call sets are from the six-tumor models. Of 935 calls, 40 (4.3%) have tumor-only assembly support.")
    resub(d, B, C, r"Every call was evaluated\.", "Every call was evaluated against the assemblies; 5 calls had no usable k-mer for the read test.")
    resub(d, B, C, r"neither \(matches neither assembly or cannot be aligned uniquely\)", "neither (the ALT context is found in no assembly)")
    resub(d, B, r"^Response\. We thank the reviewer for these insights and suggestions\.", r"396 of 541 v2\.1 calls \(73\.2%\)", "690 of 935 v2.1 calls (73.8%)")
    resub(d, B, r"^We also evaluated off-reference calls using matched-normal assemblies and reads",
          r"0\.5–10\.1% for HG008T/COLO829T and 0\.4–2\.2% for HCC1395T", "3.6–11.4% for HG008T/COLO829T and 1.8–3.9% for HCC1395T")
    resub(d, B, r"^Response Table 3\.11\. CHM13 representation",
          r"HCC1395T uses six-tumor models\. Overall, 396 of 541 calls \(73\.2%\) have a CHM13 allele, including 116 of 150 HCC1395T calls \(77\.3%\)\. Excluding HCC1395T, the total is 280 of 391 calls \(71\.6%\)\.",
          "All nine call sets are from the six-tumor models. Overall, 690 of 935 calls (73.8%) have a CHM13 allele, including "
          "269 of 325 HCC1395T calls (82.8%). Excluding HCC1395T, the total is 421 of 610 calls (69.0%).")
    resub(d, B, r"^ONT graph-coordinate calls\.", r"We admit that the ONT produced the most graph-coordinate PASS calls cannot be verified",
          "We admit that the graph-coordinate PASS calls, which are most numerous on ONT in total, cannot be verified")


# -------------------------------------------------------------------------------------------------------- figures
def _hex(seed, n):
    return hashlib.sha1(seed.encode()).hexdigest()[:n].upper()


def _clone_drawing_paragraph(d, drawing_index, tag):
    """Insert a copy of the paragraph that holds drawing `drawing_index` right after it, with its own ids."""
    el = d.drawings()[drawing_index][0]
    p = el
    while p.tag != w("p"):
        p = p.getparent()
    new = copy.deepcopy(p)
    if new.get("{%s}paraId" % W14) is not None:
        new.set("{%s}paraId" % W14, d.new_para_id())
    used = {int(x.get("id")) for x in d.doc.iter("{%s}docPr" % WP)}
    for dp in new.iter("{%s}docPr" % WP):
        n = max(used) + 1
        used.add(n)
        dp.set("id", str(n))
    for inl in new.iter("{%s}inline" % WP):
        for a in ("anchorId", "editId"):
            if inl.get("{%s}%s" % (WP14, a)) is not None:
                inl.set("{%s}%s" % (WP14, a), "5" + _hex(tag + a, 7))
    for c in new.iter("{%s}creationId" % A16):
        h = _hex(tag + "creation", 32)
        c.set("id", "{%s-%s-%s-%s-%s}" % (h[:8], h[8:12], h[12:16], h[16:20], h[20:32]))
    p.addnext(new)
    return new


def figures(d, B):
    draws = d.drawings()
    assert [t for _, _, t in draws][11:13] == ["media/image8.png", "media/image9.png"], [t for _, _, t in draws]
    B.log("Fig 1.2: " + str(replace_image(d, "media/image2.png", str(V4 / "masking/fig_1.2_masking.png"))))
    B.log("Fig 1.3: " + str(replace_image(d, "media/image3.png", str(V4 / "pathcount/fig_R1.3_pathcount.png"))))
    # Response Figure 3.3: SBS96 of the graph-coordinate calls (the ver3 file showed the COLO829T callers' spectra here)
    B.log("Fig 3.3: " + str(replace_image(d, "media/image9.png", str(V4 / "offref/fig_3.3_sbs96_graph_coordinate.png"))))
    # Response Figure 3.2: three tumors
    B.log("Fig 3.2 a: " + str(replace_image(d, "media/image8.png", str(V4 / "sbs96/fig_3.2_HG008T.png"))))
    _clone_drawing_paragraph(d, 11, "fig32c")        # -> drawings 12 and 13 use image8; the old 12 moves to 14
    _clone_drawing_paragraph(d, 11, "fig32b")
    assert [t for _, _, t in d.drawings()][11:15] == ["media/image8.png"] * 3 + ["media/image9.png"]
    B.log("Fig 3.2 b: " + str(replace_image(d, "media/image8.png", str(V4 / "sbs96/fig_3.2_COLO829T.png"), only=[12])))
    B.log("Fig 3.2 c: " + str(replace_image(d, "media/image8.png", str(V4 / "sbs96/fig_3.2_HCC1395T.png"), only=[13])))
    # Response Figure 1.2 (masking)
    J = json.loads((V4 / "masking/letter_facts.json").read_text())
    resub(d, B, r"^We evaluated this change by masking this context at inference", r"We evaluated this change.*$", J["sentence"])
    C12 = r"^Response Figure 1\.2 Evaluation of path-following context"
    resub(d, B, C12, r"This evaluation is done by single tumor model respectfully \(COLO829T and HG008T model\) on HPRC v2\.1 graph chromosome 1 SNVs",
          "This evaluation uses the six-tumor models on HPRC v2.1 (HG008T, COLO829T and HCC1395T), chromosome 1 SNVs")
    resub(d, B, C12, r"The main improvement is on the recall\.", J["caption_tail"])
    # Response Figure 1.3 (path-count channel)
    C13 = r"^Response Figure 1\.3 Effect of the haplotype-count channel\."
    resub(d, B, C13,
          r"The evaluation was performed separately using single tumor models trained on COLO829T and HG008T, with chromosome 1 SNVs restricted to benchmark regions\. For each tumor dataset, we trained",
          "The evaluation uses the single-tumor COLO829T SNV models of the current pipeline on Illumina and Fiber-seq, HPRC "
          "v2.1 and HPRC v1.1, with chromosome 1 SNVs restricted to benchmark regions; models without the channel exist "
          "only for these four datasets. For each dataset and graph, we trained")
    resub(d, B, C13, r"Blue circles indicate models trained with the haplotype-count channel, and orange diamonds indicate models trained without it\.",
          "Blue circles indicate the model trained with the haplotype-count channel (seed 0), and orange diamonds the model "
          "trained without it; the thin bar with ticks above each row marks the three seeds trained with the channel.")
    resub(d, B, C13, r"calculated as the value obtained with the haplotype-count channel minus that obtained without it\.",
          "calculated as the value obtained with the haplotype-count channel (seed 0) minus that obtained without it.")
    resub(d, B, r"^Haplotype-count channel\. We added a new reference channel",
          r"Removing channel 6 reduced SNV F1 before population filtering in all twelve comparisons and after filtering in eleven of twelve; the only exception was HG008T ONT on v2\.1\.",
          "Removing the haplotype-count channel reduced SNV F1 in all four comparisons that exist on the current pipeline "
          "(single-tumor COLO829T models, Illumina and Fiber-seq, on HPRC v1.1 and v2.1; seed 0): by 0.061–0.079 on Illumina, "
          "where the model without the channel made 1.9 times as many false-positive calls, and by 0.001–0.012 on "
          "Fiber-seq; 11 of the 12 models trained with the channel (three seeds per comparison) had a higher F1 than the "
          "model without it.")
    # Response Figures 1.1, 3.2, 3.3
    resub(d, B, r"^Response Figure 1\.1\. The eight channels represent", r"TITATATA>T", "TTATATATA>T")
    resub(d, B, r"^Response Figure 3\.2\. Revised Figure 4a",
          r"\(HG008-N or COLO829BL, Dipcall\), and PASS SNV calls after population filtering from Pansoma \(HPRC v2\.1, d46\)",
          "(HG008-N and COLO829BL, Dipcall; HCC1395BL, DeepVariant), and PASS SNV calls after population filtering from Pansoma (six-tumor models, HPRC v2.1)")
    resub(d, B, r"^Response Figure 3\.2\. Revised Figure 4a", r"shown separately for each sequencing platform\.",
          "shown separately for each sequencing platform. From top to bottom: HG008T, COLO829T and HCC1395T.")
    resub(d, B, r"^Response Figure 3\.3\. SBS96 spectra of Pansoma graph-coordinate",
          r"; 7 of 6,025 calls without a usable context were excluded\. For COLO829T ONT, a random subset of 3,000 of 4,569 calls was plotted\.",
          "; all 935 calls of the default output (Response Table 3.10; six-tumor models, HPRC v2.1) had a usable context. Panel "
          "titles give the number of calls and the cosine similarity to the chr1 somatic truth spectrum of the tumor.")


# ------------------------------------------------------------------------------------------------- other sentences
def sentences(d, B):
    T = _tables(B)
    F = T["facts"]
    # R1.2 / R1.4 / R1.6
    resub(d, B, r"^The HPRC graphs reduced matched-normal germline candidates by 84–91% while largely retaining",
          r"the benefit is especially evident before population filtering \(Response Tables 1\.1 and 1\.2\)\.",
          "the benefit is especially evident before population filtering (Response Tables 1.1 and 1.2), where the HPRC models also take population site features as an input.")
    S = r"^Somatic loss and germline suppression"
    resub(d, B, S, r"0\.2–0\.7% of truth SNVs and 3\.5–7\.0% of truth INDELs", "0.1–0.7% of truth SNVs and 1.1–7.0% of truth INDELs")
    resub(d, B, S, r"the HPRC graphs v2\.1 reduced matched-normal germline candidates by 84–91%", "the HPRC graphs v2.1 reduced matched-normal germline candidates by 85–91%")
    resub(d, B, r"^Response: We agree that sensitivity is essential in tumor-only calling", r"sensitivity below VAF 0\.08",
          "sensitivity below the candidate allele-fraction floor (0.05–0.08)")
    # R1.7
    m = F["snv_v21_vs_comparators"]
    hg, co = [m[f"HG008T {p}"] for p in ("Illumina", "PacBio", "ONT")], [m[f"COLO829T {p}"] for p in ("Illumina", "PacBio", "ONT")]
    assert all(x["four"] > max(x["deepsomatic"], x["clairs"]) for x in hg[:2]) and hg[2]["deepsomatic"] < hg[2]["four"] < hg[2]["clairs"]
    assert co[0]["four"] > max(co[0]["deepsomatic"], co[0]["clairs"]) and co[1]["deepsomatic"] < co[1]["four"] < co[1]["clairs"]
    assert co[2]["four"] < min(co[2]["deepsomatic"], co[2]["clairs"])
    assert all(x["six"] > max(x["deepsomatic"], x["clairs"]) for x in hg + co)
    resub(d, B, r"^We trained four-tumor models using HCC1395 and HCC1954",
          r"Four-tumor models achieve HG008T SNV F1 comparable to its own models and above both comparators \(DeepSomatic-TO and ClairS-TO\) but remain below both on COLO829T Fiber-seq and ONT\. Six-tumor models achieve higher SNV F1 than both comparators on HG008T and COLO829T;",
          "Four-tumor models achieve SNV F1 above both comparators (DeepSomatic-TO and ClairS-TO) on HG008T Illumina and "
          "PacBio HiFi (%.3f and %.3f) and on COLO829T Illumina (%.3f), but remain below ClairS-TO on COLO829T Fiber-seq (%.3f); "
          "on ONT they are %.3f below ClairS-TO on HG008T (%.3f) and below both comparators on COLO829T (%.3f). Six-tumor models "
          "achieve higher SNV F1 than both comparators on HG008T and COLO829T on all three platforms;"
          % (hg[0]["four"], hg[1]["four"], co[0]["four"], co[1]["four"], hg[2]["clairs"] - hg[2]["four"], hg[2]["four"], co[2]["four"]))
    # R2.3 / R3.7: counts of Response Table 3.16 (final classification; recount/a_absorbed.md)
    P = r"^Phased normal assemblies \(HG008N, COLO829BL\)\."
    resub(d, B, P, r"carries 110 HG008T alleles \(6\.8%\) and 27 COLO829T alleles \(6\.7%\); about 92% are absent",
          "carries 120 HG008T alleles (7.5%) and 35 COLO829T alleles (8.7%); about 91% are absent")
    resub(d, B, P, r"At 1,153 of the 1,351 HG008T INDEL loci", "At 1,148 of the 1,346 HG008T INDEL loci")
    resub(d, B, P, r"347 are INDELs; 26 are carried by the normal", "347 are INDELs; 34 are carried by the normal")
    resub(d, B, P, r"At 22 further INDELs", "At 15 further INDELs")
    resub(d, B, P, r"At these 22 sites, the tumor hifiasm assembly carries the truth allele at 7 and only normal-haplotype sequences at 11; 4 remain unresolved\. Counting the 11 normal-only cases and up to 4 unresolved cases as germline-like gives 37–41 of 347 INDELs \(10\.7–11\.8%\)\.",
          "At these 15 sites, the tumor hifiasm assembly carries the truth allele at 7 and only normal-haplotype sequences at 4; 4 remain unresolved. Counting the 4 normal-only cases and up to 4 unresolved cases as germline-like gives 38–42 of 347 INDELs (11.0–12.1%).")
    Q = r"^Response: We agree and have withdrawn the interpretation that these false negatives"
    resub(d, B, Q, r"Only 6\.8% and 6\.7%, respectively, are carried by the matched normal, while 92\.2% and 92\.5% are absent",
          "Only 7.5% and 8.7%, respectively, are carried by the matched normal, while 90.9% and 90.5% are absent")
    resub(d, B, Q, r"Of the 1,351 HG008T and 319 COLO829T INDEL loci absent from the normal, the tumor assemblies carry exactly the benchmark allele at 1,153 \(85\.3%\) and 226 \(70\.8%\)",
          "Of the 1,346 HG008T and 311 COLO829T INDEL loci absent from the normal, the tumor assemblies carry exactly the benchmark allele at 1,148 (85.3%) and 226 (72.7%)")
    # R2.4
    resub(d, B, r"^How are true low-VAF SNVs handled\?", r"an observed allele fraction of at least 0\.05 \(lowered from the submitted 0\.10\)",
          "an observed allele fraction of at least 0.05 (0.08 for ONT and for HG008T; lowered from the submitted 0.10)")
    # R2.6 point 2 (r26p2/r26p2.md: ranges over HPRC v1.1 and v2.1)
    resub(d, B, r"^Germline variants absent from the graph\.",
          r"On HG008T chr1, the model rejects 79–88% of these SNV candidates and the filter removes 74–89% of those called somatic\. Only 2\.4–3\.2% of the original germline SNV candidates remain after classification and filtering\. These residual calls, representing rare germline alleles absent from both the graph and the resources, account for 83–98% of false calls after filtering\.",
          "On HG008T chr1, the six-tumor models on HPRC v1.1 and v2.1, which take population site features as an input, reject "
          "97.8–98.9% of these SNV candidates; nearly all of those called somatic are alleles that the population resources do "
          "not hold, so the post-calling filter removes at most 15 more of them per call set. Only 1.1–2.2% of the original "
          "germline SNV candidates remain after classification and filtering. These residual calls, representing rare germline "
          "alleles absent from both the graph and the resources, account for 85–98% of false calls after filtering.")
    # R2.6 point 3
    G = r"^Performance without population filtering \(Response Tables 2\.9 and 2\.10\)\."
    bt = F["snv_raw_f1_v21_by_tumor"]
    resub(d, B, G, r"On HG008T, single-tumor SNV F1 is 0\.31–0\.33 on HPRC v2\.1 and about 0\.01 on the GRCh38-only graph at similar recall; on COLO829T it rises from 0\.07–0\.10 to 0\.61–0\.67\.",
          f"Before filtering, the six-tumor models reach SNV F1 {bt['HG008T']} on HG008T, {bt['COLO829T']} on COLO829T and "
          f"{bt['HCC1395T']} on HCC1395T on HPRC v2.1, against {F['snv_raw_f1_range']['linear_nositef']} for the models trained on the GRCh38-only graph.")
    resub(d, B, G, r"Because Pansoma’s SNV calls are compared across references before post-calling population filtering, the observed HPRC benefit cannot be attributed to that filter or to overlap between its resources and the test samples; it is consistent with the preferential reduction of germline candidates on the graphs \(Response Figure 2\.3\)\.",
          f"These values involve no post-calling filter. The HPRC models, however, take {SF} and the GRCh38-only models do not, "
          "so their values before and after the filter are close and the difference between references before filtering "
          "reflects the graph together with this input. In six-tumor models trained without it on all three references, the "
          f"HPRC graphs remain {F['t31_gain_hg_colo_no_sitef']} F1 above GRCh38-only on HG008T and COLO829T after filtering; the "
          "reduction of germline candidates on the graphs, which involves no model, is shown in Response Figure 2.3.")
    resub(d, B, G, r" The release six-tumor model reaches SNV F1 0\.44–0\.45 before filtering on HG008T, with higher precision but lower recall than the single-tumor models on the same HPRC v2\.1 reference\.", "")
    resub(d, B, G, r"INDELs show no consistent difference between references, and because their models take population allele frequencies as input, their unfiltered results retain population information\.",
          "INDELs show the same pattern before filtering (Response Table 2.10); their HPRC models also take population allele frequencies as input, so their unfiltered results retain population information.")
    # R3.1
    H = r"^Response: We thank the reviewer and agree that comparisons with other callers alone cannot isolate"
    resub(d, B, H, r"For HG008T and COLO829T, each across three sequencing platforms,", "For HG008T, COLO829T and HCC1395T, each across three sequencing platforms,")
    resub(d, B, H, r"and was constant on the GRCh38-only graph\. Models were trained separately for each reference with three seeds, using chr2–22 of the same tumor for training and holding out chr1 for testing\.",
          "and was constant on the GRCh38-only graph, and that the GRCh38-only models take no population site features as an "
          "input; HPRC models without these features were trained as a matched control. The six-tumor models were trained "
          "separately for each reference, holding out chr1 of the three benchmark tumors for testing.")
    resub(d, B, H, r"SNV results after population filtering are reported in Response Table 3\.1, with paired confidence intervals in Response Table 3\.2\.",
          "SNV results after population filtering are reported in Response Table 3.1, with paired confidence intervals.")
    hd = F["t31_hcc_delta"]
    assert max(abs(hd[k]) for k in ("Illumina|d9", "Illumina|d46", "PacBio|d9", "PacBio|d46")) < 0.03 and abs(hd["ONT|d9"]) < 0.005
    resub(d, B, H, r"After the same population filter, both HPRC graphs outperform GRCh38-only on all six datasets, increasing mean SNV F1 by 0\.028–0\.292, with higher precision and recall \(Response Table 3\.1\)\.",
          "After the same population filter, both HPRC graphs outperform GRCh38-only on HG008T and COLO829T on all three "
          "platforms, by %s SNV F1, and the gain remains without the site features (%s); on HCC1395T the three references are "
          "within 0.03 F1 on Illumina and PacBio HiFi, and on ONT HPRC v2.1 is %.3f lower (%.3f lower without the site features) "
          "(Response Table 3.1)." % (F["t31_gain_hg_colo_with_sitef"], F["t31_gain_hg_colo_no_sitef"], abs(hd["ONT|d46"]), abs(F["t31_hcc_ont_v21_no_sitef"])))
    r = F["snv_raw_f1_range"]
    resub(d, B, r"^Improvement before population filtering\.",
          r"Both HPRC graphs also outperform GRCh38-only before filtering, with HPRC v2\.1 giving the highest SNV F1 on all six datasets\. HPRC v2\.1 raises F1 from approximately 0\.01 to 0\.31–0\.33 on HG008T and from 0\.07–0\.10 to 0\.61–0\.67 on COLO829T\. The advantage therefore cannot be attributed to downstream population filtering \(Response Table 3\.3\)\.",
          f"Before filtering, both HPRC graphs are above GRCh38-only on all nine datasets: SNV F1 is {r['linear_nositef']} on "
          f"GRCh38-only against {r['d9']} on HPRC v1.1 and {r['d46']} on HPRC v2.1 (Response Table 3.3). No downstream "
          "population filter is involved in these values; because the HPRC models take population site features as an input "
          "and the GRCh38-only models do not, the matched comparison without these features, after filtering, is the one in "
          "Response Table 3.1.")
    resub(d, B, r"^Germline candidate reduction\.", r" The advantage also persists after removing the haplotype-count channel\.", "")
    # R3.2: the two table references were swapped
    C = r"^Candidate composition\. On HPRC v2\.1"
    resub(d, B, C, r"other short repetitive regions \(Response Table 3\.7\)", "other short repetitive regions (Response Table 3.6)")
    resub(d, B, C, r"found substantial overlap \(Response Table 3\.6\)", "found substantial overlap (Response Table 3.7)")
    # R3.8
    resub(d, B, r"^We expanded the training set by adding HCC1395T and HCC1954T", r"an allele fraction of at least 0\.08,",
          "an allele fraction of at least 0.05 to 0.08 (Response Table 3.17),")
    # R3.9
    si, ss = F["indel_sd_18"], F["snv_sd_4"]
    resub(d, B, r"^Training seeds\.",
          r"SNV reference models and final INDEL models were trained with three seeds and tested on held-out chr1 of HG008T and COLO829T \(Response Tables 3\.1 and 3\.2\)\. The standard deviation of F1 across seeds was 0\.001–0\.042 \(median 0\.012\) for HPRC SNV models and 0\.003–0\.027 \(median 0\.014\) for final INDEL models\. Across the three seeds, SNV and INDEL performance was broadly stable, with no large shifts observed \(Response Tables 3\.1 and 3\.2\)\.",
          "The final INDEL models (six-tumor and four-tumor) and the single-tumor COLO829T SNV models of the current pipeline "
          "(Illumina and Fiber-seq, HPRC v1.1 and v2.1) were trained with three seeds and tested on held-out chr1 "
          f"(Response Table 3.2). The standard deviation of F1 across seeds was {ss['range']} (median {ss['median']:.3f}) for "
          f"these SNV models and {si['range']} (median {si['median']:.3f}) for the INDEL models on the three benchmark tumors. "
          "Performance was broadly stable across seeds; the largest spreads were COLO829T Illumina SNV on HPRC v2.1 "
          "(F1 0.797 / 0.725 / 0.816), HCC1395T Illumina INDEL with the six-tumor model (0.617 / 0.534 / 0.660) and HCC1395T "
          "ONT INDEL with the four-tumor model (0.550 / 0.491 / 0.602) (Response Table 3.2).")
    # R3.10
    assert len(F["indel_six_above_both_hg_colo"]) == 4
    resub(d, B, r"^SNVs and INDELs are reported separately throughout\.", r"in three of the six HG008T and COLO829T tumor–platform settings",
          "in four of the six HG008T and COLO829T tumor–platform settings (mean INDEL F1 of three seeds)")
    resub(d, B, r"^A new Supplementary Table stratifies INDEL calls", r"errors in long homopolymers \(44–70% of HG008T false calls\)",
          "errors in long homopolymers (35–61% of HG008T false calls)")


STEPS = [
    ("mark_existing", mark_existing),
    ("nofilter", nofilter),
    ("absorbed", absorbed),
    ("nearby", nearby),
    ("three_refs", three_refs),
    ("seeds", seeds),
    ("settings", settings),
    ("regions", regions),
    ("fpdecomp", fpdecomp),
    ("offref", offref),
    ("figures", figures),
    ("sentences", sentences),
    ("cleanup", cleanup),
]
