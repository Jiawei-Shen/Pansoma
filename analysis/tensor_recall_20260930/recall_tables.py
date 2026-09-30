"""Somatic truth VCF -> tensors: how many truth records are built into tensors, how many are missing and why.

Six tensor sets (HG008T and COLO829T, three platforms each), SNV and INDEL separately. Inputs, all read-only:
  <set>/truth/somatic.graph.tsv        truth alleles of chr1-22 (multi-allelic records split), keys, PASS, in BED
  <set>/tensors/somatic.recall.tsv     per truth allele: tensor_representative / tensor_non_representative_allele /
                                       filtered (builder reasons) / no_candidate
  the truth VCFs                       record counts, chrX/Y records (no tensors: the runs are chr1-22 only),
                                       COLO829T VAF_Ill / VAF_PB / RGN
  why no candidate                     HG008T: analysis/hg008_somatic_miss_20260926/<platform>/*_no_candidate.tsv;
                                       COLO829T: colo829t_miss/<platform>/*_no_candidate.tsv (same method)
  filtered_truth/<set>.<kind>.hits.ndjson  the builder's filtered_candidates records of the filtered truth keys
  labels/<set>.tsv                     label of each represented truth's tensor (collect_labels.py)
Writes tables.md and per_truth_<sample>.tsv.
"""
import csv, json, sys
from collections import Counter, defaultdict
from pathlib import Path

import numpy as np
import pysam

csv.field_size_limit(sys.maxsize)
pysam.set_verbosity(0)  # the COLO829T VCF has no INFO header lines
HERE = Path(__file__).resolve().parent
DATA = Path("/scratch/jshen/data/pansoma_v2_tensors")
VCF = {"HG008T": "/scratch/jshen/data/HG008_GIAB/draft_v02_benchmark/HG008-T_somatic_smvar_benchmark_v0.2_tumorvariants.vcf.gz",
       "COLO829T": str(DATA / "COLO829T_truth/COLO829T_somatic_snv_indel.vcf.gz")}
# (sample, platform label, tensor set, folder of the per-truth miss classes, truth VAF field for this platform)
SETS = [("HG008T", "PacBio HiFi", "HG008T_PacBio", HERE.parent / "hg008_somatic_miss_20260926/pacbio_v6", None),
        ("HG008T", "ONT", "HG008T_ONT", HERE.parent / "hg008_somatic_miss_20260926/ont", None),
        ("HG008T", "Illumina", "HG008T_Illumina", HERE.parent / "hg008_somatic_miss_20260926/illumina", None),
        ("COLO829T", "PacBio fiberseq", "COLO829T_fiberseq", HERE / "colo829t_miss/fiberseq", "VAF_PB"),
        ("COLO829T", "ONT", "COLO829T_ONT", HERE / "colo829t_miss/ONT", "VAF_PB"),
        ("COLO829T", "Illumina", "COLO829T_Illumina", HERE / "colo829t_miss/Illumina", "VAF_Ill")]
AUTOSOMES = {f"chr{i}" for i in range(1, 23)}
GROUPS = {"SNV": ("SNP",), "INDEL": ("INS", "DEL"), "INS": ("INS",), "DEL": ("DEL",)}
CLASS_TEXT = {
    "S1_low_mapq_or_no_reads": "S1 no MAPQ>10 reads span the site (segdup / low mappability)",
    "S2_alt_absent_in_reads": "S2 reads show REF: fewer than 3 ALT reads",
    "S3_alt_is_existing_graph_allele": "S3 ALT is an existing graph SNP allele (reads take that branch, no edit)",
    "S4a_repeat_snv_edit_on_branch": "S4a repeat: the SNV becomes an edit on another branch node",
    "S4b_repeat_absorbed_by_graph_path": "S4b repeat: the ALT haplotype is absorbed by a graph path",
    "S4c_repeat_no_clear_alt": "S4c repeat: no clear ALT in the reads",
    "I1_allele_fully_in_graph": "I1 allele fully in the graph (reads take the branch / skip edge, no edit)",
    "I2_partly_in_graph_residual_edit": "I2 allele partly in the graph (only a residual edit, other spelling)",
    "I3_grch38_path_different_edit": "I3 reads stay on GRCh38 but the edits spell a different allele",
    "I4_no_or_few_alt_reads": "I4 fewer than 3 ALT reads",
    "I5_low_mapq_or_no_reads": "I5 no MAPQ>10 reads span the site",
    "I6_longer_than_50bp": "I6 the read indel is longer than 50 bp (builder limit)",
    "I7_unclear": "I7 unclear (reads on GRCh38, no local edit)",
    "I8_complex_replacement_edit": "I8 complex replacement edit (not supported by the builder)",
    "I9_truth_edit_on_unbuilt_node": "I9 reads carry the truth edit, but on a node that was not a target",
    "I10_truth_edit_on_built_node": "I10 reads carry the truth edit on a target node (support/quality below the thresholds)",
}
STATUS_ORDER = ["tensor_representative", "tensor_non_representative_allele", "filtered:min_af", "filtered:min_variants",
                "filtered:min_af,min_variants", "no_candidate"]


def fmt(n, total=None):
    if total is None:
        return f"{n:,}"
    return f"{n:,} ({100 * n / total:.1f}%)" if total else f"{n:,}"


def table(header, rows):
    """Markdown table; a column is right-aligned when its first cell is a number."""
    numeric = [str(x).lstrip("*")[:1].isdigit() or str(x) == "-" for x in rows[0]]
    out = ["| " + " | ".join(header) + " |", "|" + "|".join("---:" if n else "---" for n in numeric) + "|"]
    out += ["| " + " | ".join(str(x) for x in r) + " |" for r in rows]
    return "\n".join(out) + "\n"


# ---- truth VCF records -------------------------------------------------------------------------------------------
vcf_records, vcf_info = {}, {}
for sample, path in VCF.items():
    counts = Counter()
    with pysam.VariantFile(path) as f:
        for v in f:
            kind = "SNV" if len(v.ref) == 1 and all(len(a) == 1 for a in v.alts) else "INDEL"
            where = "autosome" if v.chrom in AUTOSOMES else "chrX/Y"
            counts[(kind, where)] += 1
            if len(v.alts) > 1:
                counts[(kind, where, "multi")] += 1
            if sample == "COLO829T":
                vcf_info[(v.chrom, str(v.pos), v.ref, v.alts[0])] = {
                    k: v.info.get(k) for k in ("VAF_Ill", "VAF_PB", "RGN")}
    vcf_records[sample] = counts

# ---- per set: truth, status, class, filtered support, label ------------------------------------------------------
truth, per = {}, {}
for sample, platform, name, miss_dir, vaf_field in SETS:
    if sample not in truth:
        truth[sample] = {r["truth_id"]: r for r in csv.DictReader(open(DATA / name / "truth/somatic.graph.tsv"), delimiter="\t")}
    t = truth[sample]
    status = {}
    for r in csv.DictReader(open(DATA / name / "tensors/somatic.recall.tsv"), delimiter="\t"):
        status[r["truth_id"]] = r["status"] + (":" + r["detail"] if r["status"] == "filtered" else "")
    classes, have_classes = {}, all((miss_dir / f).exists() for f in ("snv_no_candidate.tsv", "indel_no_candidate.tsv"))
    if have_classes:
        for f in ("snv_no_candidate.tsv", "indel_no_candidate.tsv"):
            for r in csv.DictReader(open(miss_dir / f), delimiter="\t"):
                classes[(r["chrom"], r["vcf_pos"], r["vcf_ref"], r["vcf_alt"])] = r["final_class"]
    cls = {}
    for tid, s in status.items():
        if s == "no_candidate":
            a = t[tid]
            cls[tid] = classes.get((a["chrom"], a["vcf_pos"], a["vcf_ref"], a["vcf_alt"]), "pending" if not have_classes else "unclassified")
    support = {}  # truth_id -> (alt reads, coverage, af) of its best filtered key
    key_owner = defaultdict(list)
    for tid, s in status.items():
        if s.startswith("filtered"):
            for k in t[tid]["keys"].split(","):
                key_owner[k].append(tid)
    for kind in ("SNV", "INDEL"):
        hits = HERE / "filtered_truth" / f"{name}.{kind}.hits.ndjson"
        if not hits.exists():
            continue
        for line in open(hits):
            try:
                r = json.loads(line)
            except json.JSONDecodeError:  # only while filtered_scan.sh is still writing
                print(f"WARNING: {hits.name}: unfinished line", file=sys.stderr)
                continue
            # early AF filter: upper bounds (alt_support_upper_bound, af_upper_bound); later checks: alt_count, af
            alt = r.get("alt_support_upper_bound", r.get("alt_count"))
            cov = r.get("coverage")
            af = r.get("af_upper_bound", r.get("af"))
            for tid in key_owner.get(r["candidate_id"], ()):
                if tid not in support or (alt or 0) > (support[tid][0] or 0):
                    support[tid] = (alt, cov, af)
    labels = {}
    label_file = HERE / "labels" / f"{name}.tsv"
    if label_file.exists():
        for r in csv.DictReader(open(label_file), delimiter="\t"):
            # a truth whose keys are A1 of several tensors (several placements): label 1 if any of them is 1
            if r["truth_id"] not in labels or int(r["label"]) == 1:
                labels[r["truth_id"]] = (int(r["label"]), r["reason"])
    per[name] = dict(sample=sample, platform=platform, status=status, cls=cls, support=support, labels=labels,
                     have_classes=have_classes, vaf_field=vaf_field)


def alleles(sample, group, chrom=None, in_bed=None):
    kinds = GROUPS[group]
    return [tid for tid, a in truth[sample].items() if a["kind"] in kinds
            and (chrom is None or a["chrom"] == chrom) and (in_bed is None or (a["in_bed"] == "True") == in_bed)]


def status_row(name, ids):
    st = Counter(per[name]["status"][tid] for tid in ids)
    n = len(ids)
    built = st["tensor_representative"] + st["tensor_non_representative_allele"]
    missing = n - built
    return [fmt(st["tensor_representative"], n), fmt(st["tensor_non_representative_allele"]), fmt(built, n),
            fmt(st["filtered:min_af"]), fmt(st["filtered:min_variants"]), fmt(st["filtered:min_af,min_variants"]),
            fmt(st["no_candidate"], n), fmt(missing, n)], built, missing


md = []
md.append("# Somatic truth VCF -> tensors: recall tables (generated by recall_tables.py)\n")
md.append("Status per truth allele from `<set>/tensors/somatic.recall.tsv` (key match only; partial matches are not recall). "
          "`built` = a tensor contains the truth allele: as its representative allele A1, or as another allele of the "
          "site. `filtered` = the builder made the candidate but dropped it: `AF` = allele fraction below the build floor "
          "(SNV 0.06, INDEL 0.08), `<3 reads` = fewer than 3 supporting reads (min_variants), `both` = both. "
          "`no candidate` = no read produced this allele as an edit on a target node.\n")

# ---- 1. VCF records -> tensors ------------------------------------------------------------------------------------
md.append("## 1. Truth VCF records -> tensors (all chromosomes of the VCF)\n")
rows = []
for sample, platform, name, *_ in SETS:
    for group in ("SNV", "INDEL"):
        c = vcf_records[sample]
        rec, xy, multi = c[(group, "autosome")] + c[(group, "chrX/Y")], c[(group, "chrX/Y")], c[(group, "autosome", "multi")]
        ids = alleles(sample, group)
        cells, built, missing = status_row(name, ids)
        rows.append([sample, platform, group, fmt(rec), fmt(xy), fmt(len(ids)) + (f" ({multi} multi-allelic records)" if multi else "")]
                    + cells + [f"{100 * built / len(ids):.1f}% / {100 * built / (len(ids) + xy):.1f}%"])
md.append(table(["Sample", "Platform", "Kind", "VCF records", "chrX/Y records (no tensors)", "chr1-22 truth alleles",
                 "built: representative", "built: other allele", "built total", "filtered: AF", "filtered: <3 reads",
                 "filtered: both", "no candidate", "missing total", "built % of chr1-22 / of all"], rows))
md.append("chr1-22 truth alleles: multi-allelic records are split into one allele each (HG008T: 15 INDEL records -> 30 "
          "alleles). `% of all` counts the chrX/Y truth as missing (the tensor runs are chr1-22 only).\n")

# ---- 2. INS / DEL ---------------------------------------------------------------------------------------------------
md.append("## 2. INDEL split into INS and DEL (chr1-22)\n")
rows = []
for sample, platform, name, *_ in SETS:
    for group in ("INS", "DEL"):
        ids = alleles(sample, group)
        cells, built, missing = status_row(name, ids)
        rows.append([sample, platform, group, fmt(len(ids))] + cells)
md.append(table(["Sample", "Platform", "Kind", "truth alleles", "built: representative", "built: other allele",
                 "built total", "filtered: AF", "filtered: <3 reads", "filtered: both", "no candidate", "missing total"], rows))

# ---- 3. why no candidate -------------------------------------------------------------------------------------------
md.append("## 3. Why no candidate (read-level class of every no-candidate truth allele)\n")
md.append("Each missed allele: the reads around it were fetched from the GAM and decoded with the run's frozen decoder; the "
          "class says how the ALT reads represent the allele (method of analysis/hg008_somatic_miss_20260926). "
          "Percent = of that set's no-candidate alleles of the kind.\n")
for group in ("SNV", "INS", "DEL"):
    all_classes = sorted({c for p in per.values() for tid, c in p["cls"].items()
                          if truth[p["sample"]][tid]["kind"] in GROUPS[group]},
                         key=lambda c: (c[0], int(''.join(ch for ch in c.split("_")[0][1:] if ch.isdigit()) or 0), c))
    rows = []
    totals = {}
    counts = {}
    for sample, platform, name, *_ in SETS:
        ids = [tid for tid, c in per[name]["cls"].items() if truth[sample][tid]["kind"] in GROUPS[group]]
        totals[name] = len(ids)
        counts[name] = Counter(per[name]["cls"][tid] for tid in ids)
    for c in all_classes:
        rows.append([CLASS_TEXT.get(c, c)] + [fmt(counts[name][c], totals[name]) for _, _, name, *_ in SETS])
    rows.append(["**no candidate total**"] + [f"**{totals[name]:,}**" for _, _, name, *_ in SETS])
    md.append(f"**{group}**\n")
    md.append(table(["class"] + [f"{s} {p}" for s, p, *_ in SETS], rows))

md.append("**COLO829T: truth VAF and SMaHT region of the no-candidate classes** (median truth VAF of the platform's "
          "VAF field; Extreme / Difficult / Easy = `RGN` of the truth VCF)\n")
rows = []
for sample, platform, name, *_ in SETS:
    if sample != "COLO829T" or not per[name]["have_classes"]:
        continue
    field = per[name]["vaf_field"]
    for group in ("SNV", "INS", "DEL"):
        by = defaultdict(list)
        for tid, c in per[name]["cls"].items():
            if truth[sample][tid]["kind"] in GROUPS[group]:
                a = truth[sample][tid]
                by[c].append(vcf_info[(a["chrom"], a["vcf_pos"], a["vcf_ref"], a["vcf_alt"])])
        for c in sorted(by, key=lambda c: -len(by[c])):
            v = [float(i[field]) for i in by[c]]
            rg = Counter(i["RGN"] for i in by[c])
            rows.append([platform, group, c, fmt(len(v)), f"{np.median(v):.3f}", fmt(sum(x < 0.1 for x in v)),
                         f"{rg['Extreme']:,} / {rg['Difficult']:,} / {rg['Easy']:,}"])
if rows:
    md.append(table(["Platform", "Kind", "class", "alleles", "median truth VAF", "truth VAF < 0.10",
                     "Extreme / Difficult / Easy"], rows))

# ---- 4. filtered ---------------------------------------------------------------------------------------------------
md.append("## 4. Filtered truth alleles: support the builder saw\n")
md.append("From the builder's `filtered_candidates.ndjson` record of the truth key (best key per allele): ALT reads "
          "(`alt_support_upper_bound`, or `alt_count` when a later check dropped it), coverage and AF (`af_upper_bound` / "
          "`af`, only for the AF check). Medians.\n")
rows = []
for sample, platform, name, *_ in SETS:
    for group in ("SNV", "INDEL"):
        ids = alleles(sample, group)
        for reason in ("filtered:min_af", "filtered:min_variants", "filtered:min_af,min_variants"):
            sel = [tid for tid in ids if per[name]["status"][tid] == reason]
            if not sel:
                continue
            sup = [per[name]["support"].get(tid) for tid in sel]
            alt = [s[0] for s in sup if s and s[0] is not None]
            cov = [s[1] for s in sup if s and s[1] is not None]
            af = [s[2] for s in sup if s and s[2] is not None]
            med = lambda x, d=0: (f"{np.median(x):.{d}f}" if x else "-")  # noqa: E731
            rows.append([sample, platform, group, reason.split(":")[1].replace("min_variants", "<3 reads").replace("min_af", "AF"),
                         fmt(len(sel)), fmt(len(alt)), med(alt), med(cov), med(af, 3)])
md.append(table(["Sample", "Platform", "Kind", "reason", "alleles", "with record", "median ALT reads", "median coverage",
                 "median AF"], rows))

# ---- 5. COLO829T truth VAF -----------------------------------------------------------------------------------------
md.append("## 5. COLO829T: status by truth VAF (VCF INFO; Illumina set: VAF_Ill, fiberseq and ONT: VAF_PB)\n")
edges = [0, 0.06, 0.10, 0.20, 0.30, 1.01]
names = ["<0.06", "0.06-0.10", "0.10-0.20", "0.20-0.30", ">=0.30"]
rows = []
for sample, platform, name, *_ in SETS:
    if sample != "COLO829T":
        continue
    field = per[name]["vaf_field"]
    for group in ("SNV", "INDEL"):
        ids = alleles(sample, group)
        by = defaultdict(list)
        for tid in ids:
            a = truth[sample][tid]
            v = float(vcf_info[(a["chrom"], a["vcf_pos"], a["vcf_ref"], a["vcf_alt"])][field])
            by[names[int(np.searchsorted(edges, v, side="right")) - 1]].append(tid)
        for b in names:
            st = Counter(per[name]["status"][tid].split(":")[0] for tid in by[b])
            n = len(by[b])
            if not n:  # the COLO829T INDEL truth has no VAF below 0.2
                continue
            built = st["tensor_representative"] + st["tensor_non_representative_allele"]
            rows.append([platform, group, b, fmt(n), fmt(built, n), fmt(st["filtered"], n), fmt(st["no_candidate"], n)])
md.append(table(["Platform", "Kind", "truth VAF", "alleles", "built", "filtered", "no candidate"], rows))

md.append("COLO829T: status by SMaHT region (`RGN` of the truth VCF)\n")
rows = []
for sample, platform, name, *_ in SETS:
    if sample != "COLO829T":
        continue
    for group in ("SNV", "INDEL"):
        by = defaultdict(list)
        for tid in alleles(sample, group):
            a = truth[sample][tid]
            by[vcf_info[(a["chrom"], a["vcf_pos"], a["vcf_ref"], a["vcf_alt"])]["RGN"]].append(tid)
        for region in ("Easy", "Difficult", "Extreme"):
            st = Counter(per[name]["status"][tid].split(":")[0] for tid in by[region])
            n = len(by[region])
            built = st["tensor_representative"] + st["tensor_non_representative_allele"]
            rows.append([platform, group, region, fmt(n), fmt(built, n), fmt(st["filtered"], n), fmt(st["no_candidate"], n)])
md.append(table(["Platform", "Kind", "region", "alleles", "built", "filtered", "no candidate"], rows))

# ---- 6. labels of the built truths ---------------------------------------------------------------------------------
md.append("## 6. Label of the tensor that represents the truth allele (built: representative)\n")
md.append("A built truth is a training/test positive only if its tensor is labelled 1. HG008T SNV labels use the AF floor "
          "0.08 and COLO829T Illumina SNV 0.07 (`below_snv_min_af` = -1); no INDEL floor.\n")
rows = []
for sample, platform, name, *_ in SETS:
    for group in ("SNV", "INDEL"):
        ids = [tid for tid in alleles(sample, group) if per[name]["status"][tid] == "tensor_representative"]
        lab = Counter(per[name]["labels"].get(tid, (None, "no label row")) for tid in ids)
        one = sum(n for (l, _), n in lab.items() if l == 1)
        other = "; ".join(f"{l} {r} {n:,}" for (l, r), n in sorted(lab.items(), key=lambda x: -x[1]) if l != 1)
        rows.append([sample, platform, group, fmt(len(ids)), fmt(one, len(ids)), other or "-"])
md.append(table(["Sample", "Platform", "Kind", "represented truths", "label 1", "other labels (label reason count)"], rows))

# ---- 7. chr1 -------------------------------------------------------------------------------------------------------
md.append("## 7. chr1 (the pansoma_net_v2 test chromosome): tensor recall = ceiling of any chr1 model result\n")
rows = []
for sample, platform, name, *_ in SETS:
    for group in ("SNV", "INDEL"):
        for bed in ((None, True) if sample == "HG008T" else (None,)):
            ids = alleles(sample, group, chrom="chr1", in_bed=bed)
            cells, built, missing = status_row(name, ids)
            rows.append([sample, platform, group, "chr1 in somatic BED" if bed else "chr1", fmt(len(ids))] + cells)
md.append(table(["Sample", "Platform", "Kind", "chr1 set", "chr1 truth alleles", "built: representative",
                 "built: other allele", "built total", "filtered: AF", "filtered: <3 reads", "filtered: both",
                 "no candidate", "missing total"], rows))

# ---- 8. across platforms -------------------------------------------------------------------------------------------
md.append("## 8. The same truth allele across the three platforms of a sample (chr1-22)\n")
md.append("T = built (tensor), F = filtered, N = no candidate; order of the columns above (PacBio, ONT, Illumina).\n")
rows = []
for sample in ("HG008T", "COLO829T"):
    names_ = [name for s, _, name, *_ in SETS if s == sample]
    for group in ("SNV", "INS", "DEL"):
        ids = alleles(sample, group)
        pat = Counter()
        for tid in ids:
            pat["".join("T" if per[n]["status"][tid].startswith("tensor") else "F" if per[n]["status"][tid].startswith("filtered")
                        else "N" for n in names_)] += 1
        n = len(ids)
        rows.append([sample, group, fmt(n), fmt(pat["TTT"], n), fmt(n - sum(v for k, v in pat.items() if "T" not in k), n),
                     fmt(sum(v for k, v in pat.items() if "T" not in k), n), fmt(pat["NNN"], n),
                     ", ".join(f"{k} {v:,}" for k, v in pat.most_common(6))])
md.append(table(["Sample", "Kind", "alleles", "built on all 3", "built on >=1", "built on none", "no candidate on all 3",
                 "most common patterns"], rows))

(HERE / "tables.md").write_text("\n".join(md))

# ---- per-truth tables ----------------------------------------------------------------------------------------------
for sample in ("HG008T", "COLO829T"):
    names_ = [(p, name) for s, p, name, *_ in SETS if s == sample]
    with open(HERE / f"per_truth_{sample}.tsv", "w") as out:
        head = ["truth_id", "chrom", "vcf_pos", "vcf_ref", "vcf_alt", "kind", "passed", "in_bed"]
        if sample == "COLO829T":
            head += ["VAF_Ill", "VAF_PB", "RGN"]
        for p, name in names_:
            short = name.split("_", 1)[1]
            head += [f"{short}_status", f"{short}_class", f"{short}_label", f"{short}_filtered_alt_reads", f"{short}_filtered_af"]
        out.write("\t".join(head) + "\n")
        for tid, a in truth[sample].items():
            row = [tid, a["chrom"], a["vcf_pos"], a["vcf_ref"], a["vcf_alt"], a["kind"], a["passed"], a["in_bed"]]
            if sample == "COLO829T":
                info = vcf_info[(a["chrom"], a["vcf_pos"], a["vcf_ref"], a["vcf_alt"])]
                row += [info["VAF_Ill"], info["VAF_PB"], info["RGN"]]
            for p, name in names_:
                d = per[name]
                lab = d["labels"].get(tid)
                sup = d["support"].get(tid, (None, None, None))
                row += [d["status"][tid], d["cls"].get(tid, ""), f"{lab[0]}:{lab[1]}" if lab else "",
                        "" if sup[0] is None else sup[0], "" if sup[2] is None else round(sup[2], 4)]
            out.write("\t".join(str(x) for x in row) + "\n")
print((HERE / "tables.md").read_text())
