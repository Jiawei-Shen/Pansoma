#!/wanglab/jshen/anaconda3/bin/python
"""Reviewer letter ver4.6 = the user's Reviewer_Letter_Response_ver4.5.docx (md5 713dbb63...) with the items of the
user's list of 2026-11-26 on the latest results. Plain (untracked) edits of word/document.xml with the helpers of
../../ver4_20261123/build; the existing tracked changes and comments are untouched.

  T2  Response Tables 1.2 = 2.10 = 3.4 (INDEL, no population filter): HPRC v1.1 / v2.1 columns on the six-tumor INDEL
      models without site features, seed 0                                    <- ../tables/tables46.json (t12)
  T3  Response Tables 1.5 = 2.5 = 3.15 (regions, before the filter): SNV columns and the INDEL column on the
      no-site-feature models, paired intervals recomputed                     <- ../regions/table_before_filter.md
  T7  Response Table 2.8 (HG008T INDEL by nearby germline variants): one row per reference, no site features,
      seed 0 (the "best seed" / "mean of 3 seeds" rows go)                    <- ../strata/strata_indel_nositef.md
  T9  Response Table 3.2 (per seed): the four single-tumor COLO829T SNV rows are replaced by the nine six-tumor SNV
      rows, seeds 0 / 1 / 2                                                   <- ../tables/tables46.json (t32_snv)
  T10 Response Table 3.9 (INDEL false positives by class): last column on the no-site-feature v2.1 INDEL models
                                                                              <- ../fpdecomp/table_INDEL_before_filter.md
  F   Response Figure 1.3: six-tumor SNV models with / without the haplotype-count channel, nine rows
                                                                              <- ../pathcount/fig_R1.3_six_tumor.png
  and the captions / sentences that describe these models or restate these numbers.
Not changed: T5 (Response Table 2.6 = 3.8; recomputed, identical, ../fpdecomp/checks.txt), T11 (Response Table 3.10).

usage: /wanglab/jshen/anaconda3/bin/python -I build_ver46.py [--in ver4.5.docx] [--out ver4.6.docx]
"""
import argparse
import hashlib
import json
import re
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
V46 = HERE.parent
sys.path.insert(0, str(V46.parent / "ver4_20261123/build"))
import build_ver4 as B  # noqa: E402
import steps as S  # noqa: E402
from docx_helpers import Docx, cell_text, check, delete_row, replace_image, w  # noqa: E402

LETTER = Path("/wanglab/jshen/scratch_jshen/Pansoma_papers/ver1/Pansoma_Resubmit")
IN_MD5 = "713dbb63da17e4eb393e02f0c7a6d964"
NOSF = "None of these models takes population site features as an input."
T = json.loads((V46 / "tables/tables46.json").read_text())


def t2_indel_nofilter(d):
    rows = T["t12"]["rows"]
    assert len(rows) == 9
    for n, prev in (("1.2", "1.1"), ("2.10", "2.9"), ("3.4", "3.3")):
        tbl = d.table_by_title(f"Response Table {n}")
        for tr, r in zip(d.rows(tbl)[1:], rows):
            cells = tr.findall(w("tc"))
            assert cell_text(cells[0]) == r[0] and cell_text(cells[1]) == r[1], (n, cell_text(cells[0]), cell_text(cells[1]), r)
        B.fill(d, f"Response Table {n}", [[None, None, r[2], r[3]] for r in rows], resize=False)
        P = rf"^Response Table {re.escape(n)}\. INDELs without any population filter"
        S.resub(d, B, P, r"GRCh38-only and HPRC v1\.1, one run; HPRC v2\.1, mean of three seeds;", "one run (seed 0) on each reference;")
        S.resub(d, B, P, r"The HPRC models take population site features as an input \(gnomAD and CoLoRSdb allele frequencies and "
                r"the population-filter tag of each candidate\), so their values without the filter are not free of population "
                r"information; the GRCh38-only models do not\.", NOSF)
    # Reviewer 2 comment 6 point 3: the INDEL models of the table no longer take population allele frequencies as input
    S.resub(d, B, r"^Performance without population filtering \(Response Tables 2\.9 and 2\.10\)\.",
            r"INDELs show the same pattern before filtering \(Response Table 2\.10\); their HPRC models take population allele "
            r"frequencies as input, so their unfiltered results retain population information\.",
            "INDELs show the same pattern before filtering (Response Table 2.10).")


def t3_regions(d):
    t = S.md_table(V46 / "regions/table_before_filter.md")
    rows = t[1:]
    assert len(rows) == 18 and len(rows[0]) == 9
    out = []
    for r in rows:
        tumor, plat = r[0].split(" ", 1)
        out.append([f"{tumor}\n{plat}", r[1]] + [S._two_lines(c) for c in r[2:]])
    for n in ("1.5", "2.5", "3.15"):
        tbl = d.table_by_title(f"Response Table {n}")
        for tr, r in zip(d.rows(tbl)[1:], out):                   # same datasets, same GRCh38-only side
            cells = tr.findall(w("tc"))
            assert cell_text(cells[0]) == r[0] and cell_text(cells[1]) == r[1], (n, cell_text(cells[0]))
            assert cell_text(cells[2]).split("\n")[0] == r[2].split("\n")[0], (n, r[0], cell_text(cells[2]), r[2])
            assert cell_text(cells[8]).split("\n")[0] == r[8].split("\n")[0], (n, r[0], cell_text(cells[8]), r[8])
        B.fill(d, f"Response Table {n}", out, resize=False)
        S.resub(d, B, rf"^Response Table {re.escape(n)}\. SNV calling by region",
                r"the HPRC models take population site features as an input \(population allele frequencies and the "
                r"population-filter tag of each candidate\) and the GRCh38-only models do not\.",
                "none of these models takes population site features as an input.")


def t7_indel_nearby(d):
    t = S.md_table(V46 / "strata/strata_indel_nositef.md")
    rows = t[1:]
    assert len(rows) == 9 and [r[1] for r in rows] == ["GRCh38-only", "HPRC v1.1", "HPRC v2.1"] * 3
    tbl = d.table_by_title("Response Table 2.8")
    assert len(d.rows(tbl)) == 13
    for idx in (12, 8, 4):                                        # the three "mean of 3 seeds" rows
        assert cell_text(d.rows(tbl)[idx].findall(w("tc"))[1]) == "HPRC v2.1, mean of 3 seeds"
        delete_row(tbl, idx)
    for tr, r in zip(d.rows(tbl)[1:], rows):                      # GRCh38-only rows must be the letter's
        if r[1] == "GRCh38-only":
            assert [cell_text(c) for c in tr.findall(w("tc"))[1:]] == r[1:], r
    B.fill(d, "Response Table 2.8", [[None] + r[1:] for r in rows], resize=False)
    S.resub(d, B, r"^Response Table 2\.8\. ",
            r"Six-tumor INDEL models after population filtering; GRCh38-only and HPRC v1\.1, one training seed; HPRC v2\.1, the "
            r"best of three seeds \(highest chr1 F1 on that platform\) and the mean of the three seeds;",
            "Six-tumor INDEL models after population filtering, one training seed (seed 0) per model;")
    S.resub(d, B, r"^Response Table 2\.8\. ",
            r"The HPRC models take population site features as an input and the GRCh38-only models do not\.", NOSF)


def t9_seeds(d):
    snv = T["t32_snv"]
    rows = snv["rows"]
    assert len(rows) == 9
    tbl = d.table_by_title("Response Table 3.2")
    old = d.rows(tbl)
    assert len(old) == 23
    assert all(cell_text(tr.findall(w("tc"))[2]).startswith("COLO829T only, SNV") for tr in old[19:])
    assert all(cell_text(tr.findall(w("tc"))[2]) in ("6 tumors", "4 tumors") for tr in old[1:19])
    # the 18 INDEL rows stay; their training-set cells now say INDEL, as the SNV rows say SNV
    indel = [[None, None, cell_text(tr.findall(w("tc"))[2]) + ", INDEL", None, None, None] for tr in old[1:19]]
    B.fill(d, "Response Table 3.2", indel + rows)
    S.resub(d, B, r"^Response Table 3\.2\. ",
            r"Per-seed chr1 results of the final INDEL models and of the single-tumor COLO829T SNV models\.",
            "Per-seed chr1 results (seeds 0 / 1 / 2, mean ± s.d.; after population filtering) of the final INDEL models "
            "(six-tumor and four-tumor) and of the six-tumor SNV models on HPRC v2.1.")
    sd = snv["sd"]
    big = max(sd, key=sd.get)
    assert big == "HCC1395T ONT"
    f1 = " / ".join(f"{x['f1']:.3f}" for x in snv["cells"][big])
    P = r"^Training seeds\. "
    S.resub(d, B, P,
            r"The final INDEL models \(six-tumor and four-tumor\) and the single-tumor COLO829T SNV models of the current pipeline "
            r"\(Illumina and Fiber-seq, HPRC v1\.1 and v2\.1\) were trained",
            "The final INDEL models (six-tumor and four-tumor) and the six-tumor SNV models on HPRC v2.1 were trained")
    S.resub(d, B, P, r"was 0\.003–0\.048 \(median 0\.018\) for these SNV models",
            f"was {snv['sd_range']} (median {snv['sd_median']:.3f}) for these SNV models")
    S.resub(d, B, P, r"the largest spreads were COLO829T Illumina SNV on HPRC v2\.1 \(F1 0\.797 / 0\.725 / 0\.816\),",
            f"the largest spreads were HCC1395T ONT SNV (F1 {f1}),")
    # Response Table 1.8 (training examples per model) says how many runs each model has
    S.resub(d, B, r"^Response Table 1\.8\. ", r"SNV models: one run \(seed 0\) each;",
            "SNV models: seed 0 (the six-tumor SNV models on HPRC v2.1 were also trained with seeds 1 and 2);")


def t10_indel_fp(d):
    t = S.md_table(V46 / "fpdecomp/table_INDEL_before_filter.md")
    rows = t[1:]
    assert len(rows) == 9
    tbl = d.table_by_title("Response Table 3.9")
    for tr, r in zip(d.rows(tbl)[1:], rows):                      # the five GRCh38-only / candidate columns are unchanged
        assert [cell_text(c) for c in tr.findall(w("tc"))[:5]] == r[:5], (r, [cell_text(c) for c in tr.findall(w("tc"))])
    B.fill(d, "Response Table 3.9", [[None] * 5 + [r[5]] for r in rows], resize=False)
    S.resub(d, B, r"^Response Table 3\.9\. ",
            r"Six-tumor INDEL models \(HPRC v2\.1: seed 0\); the GRCh38-graph models take no population site features as an input "
            r"and the HPRC v2\.1 models do\.",
            "Six-tumor INDEL models without population site features as an input on both references (seed 0).")


def fig13(d):
    v = json.loads((V46 / "pathcount/values.json").read_text())["rows"]
    assert len(v) == 9
    da = [r["delta"]["F1"] for r in v]
    db = [r["delta"]["F1_before"] for r in v]
    assert all(round(x, 3) > 0 for x in da + db)                  # F1 lower without the channel in 9 of 9, before and after
    assert all(round(r["delta"]["P"], 3) > 0 for r in v)
    rng = lambda x: f"{min(x):.3f}–{max(x):.3f}"                  # noqa: E731
    B.log("Fig 1.3: " + str(replace_image(d, "media/image3.png", str(V46 / "pathcount/fig_R1.3_six_tumor.png"))))
    C = r"^Response Figure 1\.3 Effect of the haplotype-count channel\."
    S.resub(d, B, C,
            r"The evaluation was performed separately using single tumor models trained on COLO829T and HG008T, with chromosome 1 "
            r"SNVs restricted to benchmark regions\. For each tumor dataset, we trained two models",
            "The evaluation uses the six-tumor models on HPRC v2.1 (HG008T, COLO829T and HCC1395T), with chromosome 1 SNVs "
            "restricted to benchmark regions. For each sequencing platform, we trained two models")
    S.resub(d, B, r"^Haplotype-count channel\. We added a new reference channel",
            r"we retrained the SNV models with and without this channel using the same dataset and pipeline\. Removing channel 6 "
            r"reduced SNV F1 before population filtering in all twelve comparisons and after filtering in eleven of twelve; the "
            r"only exception was HG008T ONT on v2\.1\. \(Response Figure 1\.3\)",
            "we retrained the six-tumor SNV models with and without this channel using the same dataset and pipeline. Removing "
            f"the channel reduced SNV F1 in all nine comparisons, before population filtering (by {rng(db)}) and after it "
            f"(by {rng(da)}), with lower precision in all nine (Response Figure 1.3).")


STEPS = (("T2", t2_indel_nofilter), ("T3", t3_regions), ("T7", t7_indel_nearby), ("T9", t9_seeds), ("T10", t10_indel_fp),
         ("Fig1.3", fig13))


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--in", dest="src", default=str(LETTER / "Reviewer_Letter_Response_ver4.5.docx"))
    ap.add_argument("--out", default=str(HERE / "Reviewer_Letter_Response_ver4.6.docx"))
    a = ap.parse_args()
    md5 = hashlib.md5(Path(a.src).read_bytes()).hexdigest()
    assert md5 == IN_MD5, f"{a.src} is not the ver4.5 this build was written for ({md5})"
    d = Docx(a.src)
    S.mark_existing(d, B)
    for name, fn in STEPS:
        before = len(B.LOG)
        fn(d)
        print(f"[{name}] {len(B.LOG) - before} operations")
    S.cleanup(d, B)
    d.save(a.out)
    print(check(a.out))
    Path(a.out).with_suffix(".buildlog.txt").write_text("\n".join(B.LOG) + "\n")
    print("wrote", a.out, "operations:", len(B.LOG))


if __name__ == "__main__":
    main()
