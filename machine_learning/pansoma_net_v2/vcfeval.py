"""Panel-of-normals filter and rtg vcfeval of a GRCh38 VCF of linear_vcf against a somatic truth VCF.

    cd machine_learning
    python -m pansoma_net_v2.vcfeval --calls <name>.SNV.linear.vcf.gz [--pon GNOMAD DBSNP 1000G COLORSDB] \
        --truth somatic.vcf.gz [--bed benchmark.bed] --sdf GRCh38.sdf --output <dir>

1. PoN (with --pon or --pon-vcf; without either only the raw calls are evaluated, as for INDELs for now): the
   repository's scripts/filter_panel_of_normals.py, with the four PoNs in its order (gnomAD, dbSNP, 1000G,
   CoLoRSdb). All match by allele; gnomAD and CoLoRSdb only at AF >= 0.001, dbSNP only non-somatic.
   Matched records get FILTER PanelOfNormals: <dir>/<name>.pon.vcf.gz. --pon-vcf reuses the tagged calls of an earlier run
   (e.g. to evaluate the same calls with and without a BED).
2. Truth: the PASS (or '.') truth records of the kind on the predicted chromosomes, written to
   <dir>/truth.<KIND>.vcf.gz, with the calls' ##contig lines and a declaration of every INFO key where the truth
   has none (COLO829T's truth declares neither).
   - SNV: a record with an ALT of REF's length. INDEL: an ALT of another length. ALL: every record.
   - The predicted chromosomes are the calls' ##pansoma_chromosomes, or --regions.
   - With --pon the PoN filter also runs over these truth records (truth_pon_tags): the share it tags bounds
     the recall of any caller after the PoN.
3. rtg vcfeval compares alleles only: --squash-ploidy --sample ALT,ALT. The truth is phased and the calls are
   not, and COLO829T's truth is sites-only.
   - Options: -f QUAL, the predicted chromosomes as --bed-regions, and -e --bed when given.
   - It runs for the calls without the PoN ("raw") and, with a PoN, with it ("pon"), twice each, under
     <dir>/vcfeval/:
     - <set>_calls: the PASS records, i.e. the call set at the model threshold. Its None row is the result.
     - <set>_curve: every PASS or LowQual record (not LowAF, not PanelOfNormals), with --all-records. This
       gives the precision/recall curve down to graph_vcf's --min-score.
4. <dir>/report.json and report.txt. For each set: the result at the threshold, the best F1 on the curve,
   precision at recall 0.5/0.8/0.9/0.95 (at the highest threshold that reaches it), and the ceiling (the
   curve's recall at the lowest score). Below the table: the truth records the PoN tags (with a PoN), and the
   PASS calls without a GRCh38 position (linear_vcf's unplaced records), which are not evaluated, with the
   precision they would give as false calls.
"""
import argparse
import gzip
import json
import math
import subprocess
import sys
from pathlib import Path

from .vcfio import REPO, meta, read_vcf, stats_path, write_vcf

PON_SCRIPT = REPO / "scripts" / "filter_panel_of_normals.py"
PON_NAMES = ("PoN1_gnomAD", "PoN2_dbSNP", "PoN3_1000G", "PoN4_CoLoRSdb")
RECALLS = (0.5, 0.8, 0.9, 0.95)
CURVE_FILTERS = {"PASS", "LowQual"}


def record_kinds(ref, alts):
    kinds = set()
    for alt in alts.split(","):
        if alt.startswith("<") or alt in ("*", "."):
            continue
        kinds.add("SNV" if len(alt) == len(ref) else "INDEL")
    return kinds


def contig_lines(header):
    return {line[len("##contig=<ID="):].split(",")[0].rstrip(">"): line for line in header
            if line.startswith("##contig=<ID=")}


def write_truth(truth, kind, chroms, output, contigs):
    """The PASS truth records of the kind on chroms; ##contig lines from `contigs` for those the truth lacks.
    Returns (records written, non-PASS records of the kind left out)."""
    header, columns, records = read_vcf(truth)
    kept, filtered = [], 0
    for r in records:
        if r[0] in chroms and (kind == "ALL" or kind in record_kinds(r[3], r[4])):
            if r[6] in ("PASS", "."):
                kept.append(r)
            else:
                filtered += 1
    own = contig_lines(header)
    declared = {line[len("##INFO=<ID="):].split(",")[0] for line in header if line.startswith("##INFO=<ID=")}
    undeclared = {}
    for r in kept:
        for item in r[7].split(";") if r[7] != "." else ():
            key, has_value, _ = item.partition("=")
            if key not in declared:
                undeclared.setdefault(key, bool(has_value))
    extra = [f'##INFO=<ID={k},Number={"." if v else 0},Type={"String" if v else "Flag"},'
             f'Description="Not declared in the truth VCF">' for k, v in sorted(undeclared.items())]
    extra += [contigs[c] for c in sorted(chroms) if c not in own and c in contigs]
    write_vcf(output, header + extra, columns, kept)
    return len(kept), filtered


def subset(calls, output, keep):
    """Copy of calls with the records whose FILTER passes keep(set of filter names)."""
    header, columns, records = read_vcf(calls)
    kept = [r for r in records if keep(set(r[6].split(";")))]
    write_vcf(output, header, columns, kept)
    return len(kept)


def regions_bed(header, chroms, path):
    lengths = {}
    for line in header:
        if line.startswith("##contig=<ID="):
            fields = dict(item.split("=", 1) for item in line[len("##contig=<"):-1].split(","))
            lengths[fields["ID"]] = int(fields["length"])
    missing = [c for c in chroms if c not in lengths]
    if missing:
        raise SystemExit(f"chromosomes {missing} have no ##contig line in the calls")
    Path(path).write_text("".join(f"{c}\t0\t{lengths[c]}\n" for c in chroms))


def run_pon(calls, output, pons):
    subprocess.run([sys.executable, str(PON_SCRIPT), str(calls), str(output), "--pon", *map(str, pons)], check=True,
                   stdout=subprocess.DEVNULL)
    return count_pon(output)


def count_pon(output):
    tagged = dict(records=0, tagged=0, tagged_pass=0, **{name: 0 for name in PON_NAMES})
    for r in read_vcf(output)[2]:
        tagged["records"] += 1
        info = dict(item.partition("=")[::2] for item in r[7].split(";"))
        if "PANSOMA_PON" in info:
            tagged["tagged"] += 1
            tagged["tagged_pass"] += r[6] == "PanelOfNormals"  # the script replaces a lone PASS
            for name in info["PANSOMA_PON"].split(","):
                tagged[name] += 1
    return tagged


def vcfeval(args, truth, calls, out, regions, curve):
    command = [args.rtg, f"RTG_MEM={args.rtg_mem}", "vcfeval", "-b", str(truth), "-c", str(calls), "-t", args.sdf,
               "-o", str(out), "--squash-ploidy", "--sample", "ALT,ALT", "--bed-regions", str(regions),
               "-f", "QUAL", "-T", str(args.threads)]
    if args.bed:
        command += ["-e", args.bed]
    if curve:
        command.append("--all-records")
    result = subprocess.run(command, capture_output=True, text=True)
    if result.returncode != 0:
        raise SystemExit(f"rtg vcfeval failed ({' '.join(command)}):\n{result.stdout}\n{result.stderr}")
    return dict(command=" ".join(command), summary=read_summary(out), roc=read_roc(out))


def number(text):
    try:
        value = float(text)
    except ValueError:
        return None
    return None if math.isnan(value) else value


def read_summary(out):
    """Rows of summary.txt: threshold (None for the un-thresholded row), tp_baseline, tp_call, fp, fn, P, R, F1."""
    rows = []
    for line in (Path(out) / "summary.txt").read_text().splitlines():
        f = line.split()
        if len(f) == 8 and f[0] != "Threshold" and f[1].isdigit():
            rows.append(dict(threshold=number(f[0]), tp_baseline=int(f[1]), tp_call=int(f[2]), fp=int(f[3]),
                             fn=int(f[4]), precision=number(f[5]), recall=number(f[6]), f1=number(f[7])))
    return rows


def read_roc(out):
    path = Path(out) / "weighted_roc.tsv.gz"
    if not path.exists():
        return []
    rows = []
    with gzip.open(path, "rt") as f:
        for line in f:
            if line.startswith("#") or not line.strip():
                continue
            s, tpb, fp, tpc, fn, p, r, f1 = line.split("\t")
            rows.append(dict(score=float(s), p_somatic=round(1 - 10 ** (-float(s) / 10), 5), tp_baseline=int(float(tpb)),
                             fp=int(float(fp)), tp_call=int(float(tpc)), fn=int(float(fn)), precision=float(p),
                             recall=float(r), f1=float(f1)))
    return rows  # highest score first


def describe(calls_run, curve_run):
    """At the threshold (the calls' None row), best F1 and precision at fixed recalls on the curve, ceiling."""
    at = next(r for r in calls_run["summary"] if r["threshold"] is None)
    roc = curve_run["roc"]
    everything = next(r for r in curve_run["summary"] if r["threshold"] is None)
    result = dict(at_threshold=at, baseline=at["tp_baseline"] + at["fn"], ceiling=everything["recall"],
                  curve_records=everything["tp_call"] + everything["fp"])
    result["best"] = max(roc, key=lambda r: r["f1"]) if roc else None
    result["precision_at_recall"] = {str(r): next((row for row in roc if row["recall"] >= r), None) for r in RECALLS}
    return result


def fmt(value, digits=3):
    return "-" if value is None else f"{value:.{digits}f}"


def table(report):
    lines = [f"{'set':<5} {'records':>8} {'TP':>5} {'FP':>6} {'FN':>5} {'P':>6} {'R':>6} {'F1':>6} | "
             f"{'best F1':>7} {'(p':>7} {'P':>6} {'R)':>6} | " + " ".join(f"{'P@R' + str(r):>7}" for r in RECALLS)
             + f" | {'ceiling':>7}"]
    sets = [name for name in ("raw", "pon") if name in report]
    for name in sets:
        d = report[name]
        a, b = d["at_threshold"], d["best"] or {}
        par = [d["precision_at_recall"][str(r)] for r in RECALLS]
        lines.append(f"{name:<5} {a['tp_call'] + a['fp']:>8,} {a['tp_baseline']:>5} {a['fp']:>6} {a['fn']:>5} "
                     f"{fmt(a['precision']):>6} {fmt(a['recall']):>6} {fmt(a['f1']):>6} | {fmt(b.get('f1')):>7} "
                     f"{fmt(b.get('p_somatic'), 4):>7} {fmt(b.get('precision')):>6} {fmt(b.get('recall')):>6} | "
                     + " ".join(f"{fmt(p['precision'] if p else None):>7}" for p in par) + f" | {fmt(d['ceiling']):>7}")
    t, u = report.get("truth_pon_tags"), report.get("unplaced_pass")
    if t:
        lines.append(f"PoN tags {t['tagged']:,} of {t['records']:,} truth records of the chromosomes "
                     f"({', '.join(f'{n} {t[n]:,}' for n in PON_NAMES)}): after the PoN any caller's recall "
                     f"<= {1 - t['tagged'] / t['records']:.3f}" if t["records"] else "No truth records for the PoN")
    if u:
        p = " ".join(f"{name} {d['tp_call'] / (d['tp_call'] + d['fp'] + u):.3f}"
                     for name, d in ((n, report[n]["at_threshold"]) for n in sets))
        lines.append(f"{u:,} PASS calls without a GRCh38 position are not evaluated; precision at the threshold "
                     f"with them as false calls: {p}" + (" (pon: if the PoN kept them all)" if "pon" in sets else ""))
    return "\n".join(lines)


def same_calls(pon_vcf, calls):
    """The PoN-tagged VCF is the calls': the same ##pansoma_ lines and IDs."""
    header, _, records = read_vcf(pon_vcf)
    calls_header, _, calls_records = read_vcf(calls)
    lines = lambda h: sorted(line for line in h if line.startswith("##pansoma_"))  # noqa: E731
    return lines(header) == lines(calls_header) and [r[2] for r in records] == [r[2] for r in calls_records]


def evaluate(args):
    calls = Path(args.calls)
    header, _, _ = read_vcf(calls)
    if meta(header, "source") != "pansoma_net_v2.linear_vcf":
        raise SystemExit(f"{calls} is not a linear_vcf output")
    kind = args.kind or meta(header, "pansoma_kind")
    chroms = args.regions or meta(header, "pansoma_chromosomes").split(",")
    out = Path(args.output)
    work = out / "vcfeval"
    if work.exists():
        raise SystemExit(f"{work} exists; give a new --output")
    work.mkdir(parents=True)
    truth = out / f"truth.{kind}.vcf.gz"
    n_truth, n_filtered = write_truth(args.truth, kind, set(chroms), truth, contig_lines(header))
    if n_truth == 0:
        raise SystemExit(f"no PASS {kind} truth records on {','.join(chroms)} in {args.truth}")
    regions = out / "regions.bed"
    regions_bed(header, chroms, regions)
    pon_vcf = pon = pon_files = truth_pon = None
    if args.pon_vcf:
        pon_vcf = Path(args.pon_vcf)
        if not same_calls(pon_vcf, calls):
            raise SystemExit(f"{pon_vcf} is not the PoN-tagged {calls}")
        pon = count_pon(pon_vcf)
        earlier = pon_vcf.parent / "report.json"
        earlier = json.loads(earlier.read_text()) if earlier.exists() else {}
        pon_files, truth_pon = earlier.get("pon_files"), earlier.get("truth_pon_tags")
    elif args.pon:
        pon_vcf = out / f"{calls.name[:-len('.vcf.gz')]}.pon.vcf.gz"
        print(f"PoN filter -> {pon_vcf}", flush=True)
        pon = run_pon(calls, pon_vcf, args.pon)
        truth_pon = run_pon(truth, out / f"truth.{kind}.pon.vcf.gz", args.pon)
        pon_files = [str(Path(p).resolve()) for p in args.pon]
    steps = {}
    for step, vcf in (("graph_vcf", meta(header, "pansoma_graph_vcf")), ("linear_vcf", calls)):
        if vcf and stats_path(vcf).exists():
            steps[step] = json.loads(stats_path(vcf).read_text())
    report = dict(calls=str(calls.resolve()), kind=kind, chromosomes=chroms, truth=str(Path(args.truth).resolve()),
                  truth_records=n_truth, truth_not_pass=n_filtered, bed=str(Path(args.bed).resolve()) if args.bed else None,
                  pon_files=pon_files, pon_vcf=str(pon_vcf.resolve()) if pon_vcf else None, pon_tags=pon,
                  truth_pon_tags=truth_pon, threshold=float(meta(header, "pansoma_threshold")),
                  min_score=float(meta(header, "pansoma_min_score")),
                  unplaced_pass=steps.get("linear_vcf", {}).get("unplaced_pass"))
    for label, source in (("raw", calls), ("pon", pon_vcf)) if pon_vcf else (("raw", calls),):
        curve_vcf = work / f"{label}_curve_input.vcf.gz"
        subset(source, curve_vcf, lambda f: f <= CURVE_FILTERS)
        print(f"rtg vcfeval {label}", flush=True)
        runs = dict(calls=vcfeval(args, truth, source, work / f"{label}_calls", regions, curve=False),
                    curve=vcfeval(args, truth, curve_vcf, work / f"{label}_curve", regions, curve=True))
        report[label] = dict(describe(runs["calls"], runs["curve"]), commands=[r["command"] for r in runs.values()],
                             curve_roc=runs["curve"]["roc"])
    report["steps"] = steps
    (out / "report.json").write_text(json.dumps(report, indent=2) + "\n")
    text = table(report)
    (out / "report.txt").write_text(text + "\n")
    print(f"{kind} on {','.join(chroms)}: {report['raw']['baseline']:,} truth variants"
          f"{' in ' + args.bed if args.bed else ''}; threshold {report['threshold']:.5f}; "
          + (f"PoN tagged {pon['tagged']:,} records ({pon['tagged_pass']:,} PASS)" if pon else "no PoN"))
    print(text, flush=True)
    return report


def parse_args(argv=None):
    p = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    p.add_argument("--calls", required=True, help="<name>.linear.vcf.gz of linear_vcf")
    pon = p.add_mutually_exclusive_group()
    pon.add_argument("--pon", nargs=4, metavar=("GNOMAD", "DBSNP", "1000G", "COLORSDB"),
                     help="the four indexed PoN VCFs, in this order (without --pon / --pon-vcf: no PoN)")
    pon.add_argument("--pon-vcf", help="the calls already tagged by the PoN filter (an earlier run's <name>.pon.vcf.gz)")
    p.add_argument("--truth", required=True, help="somatic truth VCF (bgzipped)")
    p.add_argument("--bed", help="evaluate inside these regions only (rtg -e), e.g. the benchmark BED")
    p.add_argument("--sdf", required=True, help="rtg SDF of the FASTA used by linear_vcf")
    p.add_argument("--output", required=True, help="directory; its vcfeval/ must not exist")
    p.add_argument("--kind", choices=["SNV", "INDEL", "ALL"], help="truth records to compare (default: the calls' kind)")
    p.add_argument("--regions", nargs="+", help="chromosomes (default: those predicted)")
    p.add_argument("--rtg", default="rtg")
    p.add_argument("--rtg-mem", default="4g")
    p.add_argument("--threads", type=int, default=4)
    return p.parse_args(argv)


def main(argv=None):
    evaluate(parse_args(argv))


if __name__ == "__main__":
    main()
