"""Letter Supplementary Tables R1.3a (SNV calling by region), R1.3b (long-read alignment by region) and R3.3
(graph-coordinate SNV calls per run, by validation class), from existing results only:
E16 strata (SMaHT classes, paired 1-Mb block bootstrap), E10b (GIAB v3.5 strata), E11a (G sample, long reads), and the
epoch-8 E9c / E9d / E9i tables (tmp/assessment_20261030/e9c_ep8). Writes tables.md next to this script."""
import csv, json, re
from pathlib import Path
EXP = Path("/scratch/jshen/Github/Pansoma/analysis/reviewer_response_20261011/experiments")
E9 = Path("/scratch/jshen/Github/Pansoma/tmp/assessment_20261030/e9c_ep8")
HERE = Path(__file__).resolve().parent
SETS = [("HG008T_Illumina", "HG008T Illumina"), ("HG008T_PacBio", "HG008T PacBio HiFi"), ("HG008T_ONT", "HG008T ONT"),
        ("COLO829T_Illumina", "COLO829T Illumina"), ("COLO829T_fiberseq", "COLO829T Fiber-seq"), ("COLO829T_ONT", "COLO829T ONT")]
G = {"d9": "HPRC v1.1", "d46": "HPRC v2.1"}
out = []
def tsv(p):
    return list(csv.DictReader(open(p), delimiter="\t"))
def ci(d, lo, hi):
    return f"{d:+.3f} [{lo:+.3f}, {hi:+.3f}]"
# ---- R1.3a
smaht = {(r["set"], r["arm"]): r for r in tsv(EXP / "E16_bootstrap/strata_paired.tsv") if r["stratum"] == "extreme"}
calls = {}
for l in open(EXP / "E16_bootstrap/strata.md"):
    m = re.match(r"\| (\S+) \| (linear|d9|d46) \| extreme \| ([\d,]+) \| [^|]+ \| ([\d.]+) \[", l)
    if m:
        calls[(m.group(1), m.group(2))] = (m.group(3), float(m.group(4)))
gp = {(r["set"], r["kind"], r["arm"], r["stratum"]): r for r in tsv(EXP / "E10b_giab_strata/strata_paired.tsv")}
gc = {}
for r in tsv(EXP / "E10b_giab_strata/strata_ci.tsv"):
    gc[(r["set"], r["kind"], r["arm"], r["stratum"])] = float(r["F1"])
ST = [("alldifficultregions", "all difficult"), ("alllowmapandsegdupregions", "low mappability or segdup"),
      ("AllTandemRepeatsandHomopolymers_slop5", "tandem repeats and homopolymers"), ("notinalldifficultregions", "not difficult")]
out.append("## R1.3a\n")
out.append("| chr1 | graph | SMaHT extreme: PASS calls (precision), GRCh38-only → HPRC | SMaHT extreme ΔF1 | " +
           " | ".join(f"GIAB {n} ΔF1" for _, n in ST) + " | INDEL F1, homopolymers ≥ 7 bp: GRCh38-only → HPRC |")
out.append("|---" * (4 + len(ST) + 1) + "|")
for s, lab in SETS:
    for a in ("d9", "d46"):
        r = smaht[(s, a)]; c0, c1 = calls[(s, "linear")], calls[(s, a)]
        cells = [f"{c0[0]} ({c0[1]:.3f}) → {c1[0]} ({c1[1]:.3f})", ci(float(r["dF1"]), float(r["dF1_lo"]), float(r["dF1_hi"]))]
        for k, _ in ST:
            q = gp[(s, "SNV", a, k)]; cells.append(ci(float(q["dF1"]), float(q["dF1_lo"]), float(q["dF1_hi"])))
        k = "AllHomopolymers_ge7bp_imperfectge11bp_slop5"
        cells.append(f"{gc[(s, 'INDEL', 'linear', k)]:.3f} → {gc[(s, 'INDEL', a, k)]:.3f}")
        out.append(f"| {lab} | {G[a]} | " + " | ".join(cells) + " |")
# ---- R1.3b
out.append("\n## R1.3b\n")
out.append("| long reads | region (reads) | mismatched bp per kb: GRCh38-only / v1.1 / v2.1 | same placement as GRCh38-only: v1.1 / v2.1 | MAPQ ≥ 30 lost: v1.1 / v2.1 | MAPQ ≥ 30 gained: v1.1 / v2.1 |")
out.append("|---|---|---|---|---|---|")
for s, lab in [("HG008T_PacBio", "HG008T PacBio HiFi"), ("HG008T_ONT", "HG008T ONT"), ("COLO829T_fiberseq", "COLO829T Fiber-seq"), ("COLO829T_ONT", "COLO829T ONT")]:
    d = json.load(open(EXP / f"E11a_alignment_stats/results/{s}.json"))["G"]["strata"]
    for reg in ("SMaHT easy", "SMaHT difficult", "SMaHT extreme"):
        pa, pr = d[reg]["per_arm"], d[reg]["pairs"]
        f = lambda p, k: pr[p][k]
        out.append(f"| {lab} | {reg.replace('SMaHT ', '')} ({d[reg]['reads']:,}) | " +
                   " / ".join(f"{pa[a]['mm_per_kb']:.2f}" for a in ("linear", "d9", "d46")) + " | " +
                   " / ".join(f"{f(p, 'concordant_pct'):.2f}%" for p in ("d9_vs_linear", "d46_vs_linear")) + " | " +
                   " / ".join(f"{f(p, 'mapq30_lost_pct'):.1f}%" for p in ("d9_vs_linear", "d46_vs_linear")) + " | " +
                   " / ".join(f"{f(p, 'mapq30_gained_pct'):.1f}%" for p in ("d9_vs_linear", "d46_vs_linear")) + " |")
# ---- R3.3
T = open(E9 / "EXP/tables.md").read().split("\n## ")
sec = {x.split(".")[0]: x for x in T[1:]}
rows = lambda x: [[c.strip() for c in l.strip("|").split("|")] for l in x.split("\n") if l.startswith("| rcpu")]
RUN = {"HG008T_Illumina": ("HG008_Illumina_SNV_7fc7396_ddp4_b1024", "HG008_Illumina_SNV_d46_7fc7396"),
       "HG008T_PacBio": ("HG008_PacBio_SNV_7fc7396", "HG008_PacBio_SNV_d46_7fc7396"),
       "HG008T_ONT": ("HG008_ONT_SNV_7fc7396", "HG008_ONT_SNV_d46_7fc7396"),
       "COLO829T_Illumina": ("COLO829_Illumina_SNV_7fc7396_ddp4_b1024", "COLO829_Illumina_SNV_d46_7fc7396"),
       "COLO829T_fiberseq": ("COLO829_fiberseq_SNV_7fc7396_ddp4_b1024", "COLO829_fiberseq_SNV_d46_7fc7396"),
       "COLO829T_ONT": ("COLO829_ONT_SNV_7fc7396_ddp4_b1024", "COLO829_ONT_SNV_d46_7fc7396")}
s1 = {r[0]: r for r in rows(sec["1"])}
s3 = {(r[0], r[1]): r for r in rows(sec["3"])}
s8 = {r[1]: r for r in rows(sec["8"])}
e9i = {}
for l in open(E9 / "code_e9i/tables.md"):
    c = [x.strip() for x in l.strip().strip("|").split("|")]
    if len(c) > 7 and c[1] == "unplaced_pass":
        e9i[c[0]] = (c[3], c[7])
num = lambda x: int(re.match(r"(\d+)", x).group(1))
out.append("\n## R3.3\n")
out.append("| chr1 SNV | graph | graph-coordinate PASS calls (evaluated) | population-represented | matched-normal supported | tumor-only | tumor-only by the assemblies alone | ALT in ≥ 2 normal reads | no normal ALT read, ALT in tumor reads | also PASS on both other platforms |")
out.append("|---|---|---|---|---|---|---|---|---|---|")
for s, lab in SETS:
    for a, run in zip(("d9", "d46"), RUN[s]):
        k = f"rcpu__{run}__vcf_chr1"; r = s1[k]
        n_eval = num(r[8])
        asm = s3[(k, "unplaced_pass")]; m = re.search(r"tumor_only (\d+)", asm[4]); ta = int(m.group(1)) if m else 0
        d8 = json.loads(s8[k][5]); all3 = sum(v.get("2", 0) for v in d8.values())
        nr = e9i.get(run, ("–", "–"))
        ev = f"{num(r[7]):,}" + (f" ({n_eval:,})" if n_eval != num(r[7]) else "")
        out.append(f"| {lab} | {G[a]} | {ev} | {r[10]} | {r[13]} | {r[14]} | {ta} ({100 * ta / int(asm[2]):.1f} %) | {nr[0]} | {nr[1]} | {all3} |")
(HERE / "tables.md").write_text("\n".join(out) + "\n")
print("\n".join(out))
