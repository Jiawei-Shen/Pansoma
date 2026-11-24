"""ver4 offref (3): graph-coordinate share of all PASS SNV calls per six-tumor call set (HPRC v2.1, chr1, default output).

Two denominators, both whole chr1 (no benchmark BED) and after the GRCh38 population filter (the filter cannot tag a
graph-coordinate call, so all of them stay):
  share_rtg      U / (TP_call + FP + U), TP_call / FP = rtg vcfeval of the placed PASS calls without a BED
                 (<vcfdir>/eval_nobed/report.json, pon.at_threshold). The definition of the ver3 sentence for
                 HG008T / COLO829T ("0.5-10.1%"; points/R3.4.md, 2026-10-21 note: U / (TP + FP + U), after PoN, no BED).
  share_records  U / (placed PASS records after the filter + U), records counted in the PoN-tagged linear VCF
                 (<vcfdir>/eval_bed/*.linear.pon.vcf.gz, FILTER == PASS; the file holds every chr1 record). The
                 definition of the ver3 numbers for HCC1395T ("0.4-2.2%"; tmp/r33_hcc1395_20261104/r34/r34_results.py,
                 pass_counts.placed_pon of E9a evaluate_complete.py). Independent of rtg and of the truth set.
U = PASS records of <vcfdir>/*.linear.unplaced.vcf.gz (asserted equal to report.json unplaced_pass of both evaluations).
share_raw_records = the same before the population filter (PASS records of *.linear.vcf.gz).

usage: python share.py            -> share.md, share.json (the nine call sets of runs.txt)
       python share.py --repro    -> prints the two definitions on the ver3 call sets (old numbers; writes nothing)
"""
import gzip
import json
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
R = Path("/scratch/jshen/data/pansoma_net_v2_runs")
PLAT = {"Illumina": "Illumina", "PacBio_HiFi": "PacBio HiFi", "Fiber-seq": "Fiber-seq", "ONT": "ONT"}
REPRO = [("HG008T", "Illumina", R / "reviewer_cpu_20261013/HG008_Illumina_SNV_d46_7fc7396__norescue/vcf_chr1"),
         ("HG008T", "PacBio_HiFi", R / "reviewer_cpu_20261013/HG008_PacBio_SNV_d46_7fc7396__norescue/vcf_chr1"),
         ("HG008T", "ONT", R / "reviewer_cpu_20261013/HG008_ONT_SNV_d46_7fc7396__norescue/vcf_chr1"),
         ("COLO829T", "Illumina", R / "reviewer_cpu_20261013/COLO829_Illumina_SNV_d46_7fc7396__norescue/vcf_chr1"),
         ("COLO829T", "Fiber-seq", R / "reviewer_cpu_20261013/COLO829_fiberseq_SNV_d46_7fc7396__norescue/vcf_chr1"),
         ("COLO829T", "ONT", R / "reviewer_cpu_20261013/COLO829_ONT_SNV_d46_7fc7396__norescue/vcf_chr1"),
         ("HCC1395T", "Illumina", R / "r33_hcc1395_20261104/mc6_Illumina_SNV_norescue/vcf_chr1"),
         ("HCC1395T", "PacBio_HiFi", R / "r33_hcc1395_20261104/mc6_PacBio_SNV_norescue/vcf_chr1"),
         ("HCC1395T", "ONT", R / "r33_hcc1395_20261104/mc6_ONT_SNV_norescue/vcf_chr1")]


def count_pass(path):
    n = 0
    with gzip.open(path, "rt") as f:
        for line in f:
            if line[0] != "#" and line.split("\t", 8)[6] == "PASS":
                n += 1
    return n


def one(vdir):
    vdir = Path(vdir)
    (unplaced,) = vdir.glob("*.linear.unplaced.vcf.gz")
    (linear,) = vdir.glob("*.linear.vcf.gz")
    nobed = json.loads((vdir / "eval_nobed/report.json").read_text())
    bed = json.loads((vdir / "eval_bed/report.json").read_text())
    u = count_pass(unplaced)
    assert u == nobed["unplaced_pass"] == bed["unplaced_pass"], (vdir, u, nobed["unplaced_pass"], bed["unplaced_pass"])
    assert nobed["bed"] is None and bed["bed"], vdir
    assert Path(nobed["calls"]).name == linear.name
    pon_vcf = vdir / "eval_bed" / Path(bed["pon_vcf"]).name
    at = nobed["pon"]["at_threshold"]
    raw = nobed["raw"]["at_threshold"]
    assert at["tp_baseline"] + at["fn"] == nobed["pon"]["baseline"]
    placed_pon, placed_raw = count_pass(pon_vcf), count_pass(linear)
    return dict(vcf_dir=str(vdir), unplaced_pass=u, tp_call=at["tp_call"], fp=at["fp"], placed_pass_records_pon=placed_pon,
                placed_pass_records_raw=placed_raw, raw_tp_call=raw["tp_call"], raw_fp=raw["fp"],
                share_rtg=u / (at["tp_call"] + at["fp"] + u), share_records=u / (placed_pon + u),
                share_raw_records=u / (placed_raw + u), truth_nobed=nobed["pon"]["baseline"])


def main():
    repro = "--repro" in sys.argv
    if repro:
        runs = REPRO
    else:
        runs = [(l.split()[2], l.split()[3], l.split()[1]) for l in open(HERE / "runs.txt")]
    rows = [dict(tumor=t, platform=PLAT[p], **one(v)) for t, p, v in runs]
    pct = lambda x: f"{100 * x:.1f}%"
    md = ["| Tumor | Platform | Graph-coordinate PASS (U) | Placed PASS, rtg TP + FP (no BED, after filter) | share U / (TP + FP + U) "
          "| Placed PASS records (after filter) | share U / (records + U) | Placed PASS records (before filter) | share before filter |",
          "|---|---|---|---|---|---|---|---|---|"]
    for r in rows:
        md.append(f"| {r['tumor']} | {r['platform']} | {r['unplaced_pass']} | {r['tp_call'] + r['fp']:,} | {pct(r['share_rtg'])} | "
                  f"{r['placed_pass_records_pon']:,} | {pct(r['share_records'])} | {r['placed_pass_records_raw']:,} | "
                  f"{pct(r['share_raw_records'])} |")
    rng = {}
    for label, sel in (("HG008T / COLO829T", ("HG008T", "COLO829T")), ("HCC1395T", ("HCC1395T",)),
                       ("all nine", ("HG008T", "COLO829T", "HCC1395T")), ("HG008T", ("HG008T",)), ("COLO829T", ("COLO829T",))):
        s = [r for r in rows if r["tumor"] in sel]
        rng[label] = {k: [min(r[k] for r in s), max(r[k] for r in s)] for k in ("share_rtg", "share_records")}
    md += ["", "Ranges (min-max over the call sets):", ""]
    for label, d in rng.items():
        md.append(f"- {label}: U / (TP + FP + U) {pct(d['share_rtg'][0])}-{pct(d['share_rtg'][1])}; "
                  f"U / (PASS records + U) {pct(d['share_records'][0])}-{pct(d['share_records'][1])}")
    text = "\n".join(md) + "\n"
    print(text)
    if not repro:
        head = ("# Graph-coordinate share of all PASS SNV calls (six-tumor current models, HPRC v2.1, chr1, default output)\n\n"
                "Whole chr1 (no benchmark BED), after the GRCh38 population filter. Definitions in share.py.\n\n")
        (HERE / "share.md").write_text(head + text)
        (HERE / "share.json").write_text(json.dumps(dict(rows=rows, ranges=rng), indent=1) + "\n")


if __name__ == "__main__":
    main()
