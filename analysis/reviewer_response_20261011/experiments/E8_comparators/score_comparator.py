"""Score a comparator caller's tumor-only VCF (ClairS-TO, DeepSomatic) on chr1 with Pansoma's protocol.

    /wanglab/jshen/anaconda3/bin/python score_comparator.py --vcf snv.vcf.gz [--vcf indel.vcf.gz] \
        --set HG008T_PacBio --kind SNV --caller clairs_to --output DIR [--threads 4]

The protocol is machine_learning/pansoma_net_v2/vcfeval.py, used as a library (not modified): vcfeval.py itself only
accepts linear_vcf outputs (it reads ##pansoma_kind / ##pansoma_chromosomes / ##pansoma_threshold), so this script
prepares the comparator's records and then calls the same functions:
  - truth: vcfeval.write_truth (PASS truth records of the kind on chr1), truth VCF and BED from the tensor set's
    tensors/<KIND>/labels.manifest.json (the same files Pansoma's chr1 test uses);
  - PoN: vcfeval.run_pon = scripts/filter_panel_of_normals.py with the four PoN VCFs of jobs/vcf.sh (SNV gnomAD /
    CoLoRSdb AF >= 1e-4, dbSNP non-somatic, 1000G; INDEL gnomAD / CoLoRSdb, AF >= 0.01 in homopolymers >= 7 bp, else
    >= 0.001), FASTA = the labels' FASTA; the truth records are tagged too (the recall bound of any caller after the PoN);
  - rtg: vcfeval.vcfeval = rtg vcfeval --squash-ploidy --sample ALT,ALT -f QUAL, chr1 as --bed-regions, -e BED for
    eval_bed and no -e for eval_nobed, --all-records for the curves; vcfeval.describe for the result at the
    caller's PASS, the best F1 on the curve, precision at recall 0.5/0.8/0.9/0.95 and the curve's ceiling.

Preparation of the calls (<output>/calls.all.vcf.gz):
  1. bcftools norm -m -any -f FASTA (split multi-allelic records, left-align; linear_vcf left-aligns Pansoma's calls
     too, and the PoN matches exact alleles), chr1 only;
  2. records whose ALT is a real allele of the kind (SNV: ALT of REF's length, INDEL: of another length; the
     vcfeval.record_kinds rule that also selects the truth); RefCall / NoCall records (the caller's "no variant")
     are left out; one record per (POS, REF, ALT), the highest QUAL kept;
  3. sites-only VCF: QUAL and FILTER of the caller ('.' QUAL -> 0), INFO TOOL_AF / TOOL_DP = the first sample's
     FORMAT VAF or AF / DP (for VAF strata later), ##contig from the FASTA index.

Four call sets, each scored inside the BED (eval_bed) and on the whole chromosome (eval_nobed):
  pass      the caller's PASS records                              (setting 1: tool PASS)
  pass_pon  pass minus records the PoN tags                         (setting 2: tool PASS + our PoN)
  all       every record (any FILTER), curve by QUAL                (curve without our PoN)
  all_pon   all minus records the PoN tags, curve by QUAL           (setting 3: all calls + our PoN, tool's quality)
'all' and 'all_pon' report the curve; their at_threshold row is the one of 'pass' and 'pass_pon'.
"""
import argparse
import collections
import json
import os
import subprocess
import sys
from pathlib import Path

REPO = Path("/scratch/jshen/Github/Pansoma")
sys.path.insert(0, str(REPO / "machine_learning"))
from pansoma_net_v2 import vcfeval as V  # noqa: E402
from pansoma_net_v2.metrics import RECALLS  # noqa: E402
from pansoma_net_v2.vcfio import read_vcf, write_vcf  # noqa: E402

TENSORS = Path("/scratch/jshen/data/pansoma_v2_tensors")
PON_DIR = Path("/scratch/jshen/data/Pansoma/panel_of_normal_VCFs")
PONS = [PON_DIR / "af-only-gnomad.hg38.vcf.gz", PON_DIR / "Homo_sapiens_assembly38.dbsnp138.vcf.gz",
        PON_DIR / "1000g_pon.hg38.vcf.gz", PON_DIR / "CoLoRSdb.GRCh38.v1.1.0.deepvariant.glnexus.vcf.gz"]
SDF = "/scratch/jshen/REF/GRCh38.primary.sdf"
RTG = "/scratch/jshen/bin/rtg-tools-3.12.1/rtg"
NO_CALL_FILTERS = {"RefCall", "NoCall"}
SETS = ("pass", "pass_pon", "all", "all_pon")


FASTA = "/scratch/jshen/data/HapMap/GCA_000001405.15_GRCh38_no_alt_analysis_set.fasta"  # the labels' FASTA
CASTLE_TRUTH = {  # sets not yet built (CASTLE Illumina): the truth files their labels will use (CASTLE/README.txt)
    line: (f"/scratch/jshen/data/CASTLE/{line}/truth/{line}_somatic_CASTLE.vcf.gz",
           f"/scratch/jshen/data/CASTLE/{line}/truth/{line}_high-confidence-regions_DeepSomatic-HCC1395-model.bed")
    for line in ("HCC1954", "H1437", "H2009")}


def labels_manifest(set_name, kind):
    path = TENSORS / set_name / "tensors" / kind / "labels.manifest.json"
    if not path.exists():
        line = set_name.split("T_")[0]
        if line in CASTLE_TRUTH:
            print(f"no {path}: CASTLE truth of {line}", flush=True)
            return (*CASTLE_TRUTH[line], FASTA)
        raise SystemExit(f"no {path}; give --truth and --bed")
    m = json.loads(path.read_text())
    return m["truth"]["somatic"]["vcf"], m["truth"]["somatic"]["bed"], m["fasta"]


def fasta_lengths(fasta):
    lengths = {}
    for line in Path(str(fasta) + ".fai").read_text().splitlines():
        f = line.split("\t")
        lengths[f[0]] = int(f[1])
    return lengths


def normalized_records(vcf, fasta, chrom):
    """bcftools norm -m -any -f FASTA of one VCF, chrom only: (FILTER header lines, records as column lists)."""
    command = ["bcftools", "norm", "-m", "-any", "-f", str(fasta), "-c", "w", "-t", chrom, "-O", "v", str(vcf)]
    result = subprocess.run(command, capture_output=True, text=True)
    if result.returncode != 0:
        raise SystemExit(f"{' '.join(command)} failed:\n{result.stderr}")
    filters, records = [], []
    for line in result.stdout.splitlines():
        if line.startswith("##FILTER="):
            filters.append(line)
        elif line and not line.startswith("#"):
            records.append(line.split("\t"))
    return filters, records, " ".join(command)


def sample_value(record, keys):
    if len(record) < 10:
        return None
    fmt, sample = record[8].split(":"), record[9].split(":")
    for key in keys:
        if key in fmt and fmt.index(key) < len(sample):
            value = sample[fmt.index(key)]
            if value not in (".", ""):
                return value.split(",")[0]
    return None


def prepare_calls(vcfs, kind, fasta, chrom, output, caller):
    lengths = fasta_lengths(fasta)
    filter_lines, best, commands = {}, {}, []
    counts = collections.Counter()
    for vcf in vcfs:
        filters, records, command = normalized_records(vcf, fasta, chrom)
        commands.append(command)
        for line in filters:
            filter_lines.setdefault(line.split(",")[0], line)
        for r in records:
            ref, alt = r[3].upper(), r[4].upper()
            counts["records_after_norm"] += 1
            if alt in ("*", ".", "<*>", "<NON_REF>") or alt.startswith("<") or alt == ref:
                counts["no_alt"] += 1
                continue
            if kind not in V.record_kinds(ref, alt):
                counts["other_kind"] += 1
                continue
            filt = r[6]
            if set(filt.split(";")) & NO_CALL_FILTERS:
                counts["refcall_or_nocall"] += 1
                continue
            qual = r[5] if r[5] not in (".", "") else "0"
            info = []
            af, dp = sample_value(r, ("VAF", "AF")), sample_value(r, ("DP",))
            if af is not None:
                info.append(f"TOOL_AF={af}")
            if dp is not None:
                info.append(f"TOOL_DP={dp}")
            row = [r[0], r[1], ".", ref, alt, qual, filt if filt else ".", ";".join(info) or "."]
            key = (int(r[1]), ref, alt)
            if key in best:
                counts["duplicates"] += 1
                if float(qual) <= float(best[key][5]):
                    continue
            best[key] = row
    header = ["##fileformat=VCFv4.2", f"##source=E8_comparators.score_comparator ({caller})",
              f"##reference=file://{fasta}", f"##contig=<ID={chrom},length={lengths[chrom]}>"]
    header += sorted(filter_lines.values())
    if '##FILTER=<ID=PASS' not in "".join(header):
        header.append('##FILTER=<ID=PASS,Description="All filters passed">')
    header += ['##INFO=<ID=TOOL_AF,Number=1,Type=Float,Description="Allele fraction of the caller (FORMAT VAF or AF)">',
               '##INFO=<ID=TOOL_DP,Number=1,Type=Integer,Description="Depth of the caller (FORMAT DP)">']
    columns = "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO"
    rows = [best[k] for k in sorted(best)]
    declared = {line.split(",")[0][len("##FILTER=<ID="):] for line in header if line.startswith("##FILTER=")}
    for row in rows:  # FILTER values the input did not declare
        for name in row[6].split(";"):
            if name not in declared and name != ".":
                header.append(f'##FILTER=<ID={name},Description="Filter of the caller (not declared in its VCF)">')
                declared.add(name)
    write_vcf(output, header, columns, rows)
    counts["records"] = len(rows)
    counts["filters"] = dict(collections.Counter(row[6] for row in rows))
    return header, dict(counts), commands


def subset(source, output, keep):
    return V.subset(source, output, keep)


def evaluate(args):
    out = Path(args.output)
    out.mkdir(parents=True, exist_ok=True)
    truth_vcf, bed, fasta = labels_manifest(args.set, args.kind)
    if args.truth:
        truth_vcf = args.truth
    if args.bed:
        bed = args.bed
    chroms = [args.chrom]
    calls_all = out / "calls.all.vcf.gz"
    header, counts, norm_commands = prepare_calls(args.vcf, args.kind, fasta, args.chrom, calls_all, args.caller)
    calls = {"all": calls_all, "pass": out / "calls.pass.vcf.gz"}
    counts["pass"] = subset(calls_all, calls["pass"], lambda f: f == {"PASS"})
    truth = out / f"truth.{args.kind}.vcf.gz"
    n_truth, n_truth_filtered = V.write_truth(truth_vcf, args.kind, set(chroms), truth, V.contig_lines(header))
    regions = out / "regions.bed"
    V.regions_bed(header, chroms, regions)
    pon = {}
    pon["calls_all"] = V.run_pon(calls_all, out / "calls.all.pon.vcf.gz", PONS, fasta)
    truth_pon = V.run_pon(truth, out / f"truth.{args.kind}.pon.vcf.gz", PONS, fasta)
    calls["all_pon"] = out / "calls.all_pon.vcf.gz"
    counts["all_pon"] = subset(out / "calls.all.pon.vcf.gz", calls["all_pon"], lambda f: "PanelOfNormals" not in f)
    calls["pass_pon"] = out / "calls.pass_pon.vcf.gz"
    counts["pass_pon"] = subset(out / "calls.all.pon.vcf.gz", calls["pass_pon"], lambda f: f == {"PASS"})
    report = dict(caller=args.caller, set=args.set, kind=args.kind, chromosomes=chroms, input_vcfs=[str(v) for v in args.vcf],
                  normalization=norm_commands, truth=str(truth_vcf), bed=str(bed), fasta=str(fasta),
                  pon_files=[str(p) for p in PONS], truth_records=n_truth, truth_not_pass=n_truth_filtered,
                  truth_pon_tags=truth_pon, pon_tags=pon, counts=counts)
    rtg_args = argparse.Namespace(rtg=RTG, rtg_mem=args.rtg_mem, sdf=SDF, threads=args.threads, bed=None)
    for mode, mode_bed in (("eval_bed", bed), ("eval_nobed", None)):
        rtg_args.bed = mode_bed
        work = out / mode
        work.mkdir(exist_ok=True)
        runs = {}
        for name in SETS:
            target = work / name
            if target.exists():
                raise SystemExit(f"{target} exists; give a new --output")
            runs[name] = V.vcfeval(rtg_args, truth, calls[name], target, regions, curve=name.startswith("all"))
        result = {}
        for name in SETS:
            at_run = runs[name.replace("all", "pass")] if name.startswith("all") else runs[name]
            d = V.describe(at_run, runs[name])
            d["command"] = runs[name]["command"]
            d["curve_roc"] = runs[name]["roc"]
            result[name] = d
        report[mode] = result
    (out / "report.json").write_text(json.dumps(report, indent=1) + "\n")
    text = table(report)
    (out / "report.txt").write_text(text + "\n")
    print(text, flush=True)
    return report


def fmt(value, digits=3):
    return "-" if value is None else f"{value:.{digits}f}"


def table(report):
    lines = [f"{report['caller']} {report['set']} {report['kind']} chr1; truth {report['truth_records']:,} records; "
             f"calls {report['counts']['records']:,} (PASS {report['counts']['pass']:,}, PASS after PoN "
             f"{report['counts']['pass_pon']:,}, all after PoN {report['counts']['all_pon']:,})"]
    t = report["truth_pon_tags"]
    if t and t["records"]:
        lines.append(f"PoN tags {t['tagged']:,} of {t['records']:,} chr1 truth records: after the PoN any caller's "
                     f"recall <= {1 - t['tagged'] / t['records']:.3f}")
    for mode in ("eval_bed", "eval_nobed"):
        lines.append(f"== chr1 {mode}")
        lines.append(f"{'set':<9} {'TP':>5} {'FP':>7} {'FN':>5} {'P':>6} {'R':>6} {'F1':>6} | {'best F1':>7} {'(QUAL':>7} "
                     f"{'P':>6} {'R':>6}) | " + " ".join(f"{'P@R' + str(r):>7}" for r in RECALLS) + f" | {'ceiling':>7}")
        for name in SETS:
            d = report[mode][name]
            a, b = d["at_threshold"], d["best"] or {}
            par = [d["precision_at_recall"][str(r)] for r in RECALLS]
            lines.append(f"{name:<9} {a['tp_baseline']:>5} {a['fp']:>7} {a['fn']:>5} {fmt(a['precision']):>6} "
                         f"{fmt(a['recall']):>6} {fmt(a['f1']):>6} | {fmt(b.get('f1')):>7} {fmt(b.get('score'), 1):>7} "
                         f"{fmt(b.get('precision')):>6} {fmt(b.get('recall')):>6} | "
                         + " ".join(f"{fmt(p['precision'] if p else None):>7}" for p in par) + f" | {fmt(d['ceiling']):>7}")
    return "\n".join(lines)


def parse_args(argv=None):
    p = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    p.add_argument("--vcf", action="append", required=True, help="caller VCF (repeat for ClairS-TO's snv + indel)")
    p.add_argument("--set", required=True, help="tensor set whose labels manifest gives truth, BED and FASTA")
    p.add_argument("--kind", required=True, choices=["SNV", "INDEL"])
    p.add_argument("--caller", required=True)
    p.add_argument("--output", required=True)
    p.add_argument("--chrom", default="chr1")
    p.add_argument("--truth", help="override the manifest's truth VCF")
    p.add_argument("--bed", help="override the manifest's BED")
    p.add_argument("--threads", type=int, default=int(os.environ.get("SLURM_CPUS_PER_TASK", 4)), help="rtg threads (default: the job's CPUs, else 4)")
    p.add_argument("--rtg-mem", default="4g")
    return p.parse_args(argv)


if __name__ == "__main__":
    evaluate(parse_args())
