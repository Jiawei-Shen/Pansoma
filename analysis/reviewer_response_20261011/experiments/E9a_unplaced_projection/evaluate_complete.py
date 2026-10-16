"""E9e: complete-output metrics of one run's vcf_chr1 after the E9a projection.

    python evaluate_complete.py --vcf-dir <run>/vcf_chr1 --projection <scratch>/<RUN>/<vcf_dir_name> [--threads 4]

1. rtg vcfeval of complete.linear.vcf.gz (the placed calls + the class-i projections) with the repository's evaluator
   (machine_learning/pansoma_net_v2/vcfeval.py, unmodified, as a CLI): the same truth, BED, SDF, the four PoN VCFs
   (scripts/filter_panel_of_normals.py, the run's rules) and the same QUAL/threshold as the run's own eval:
   <projection>/eval_bed (PoN run on all records) and <projection>/eval_nobed (reuses the PoN-tagged calls).
2. Class-ii calls (an interval only): a PASS class-ii call is a TP when a truth record of the kind that the rtg run left
   as FN lies within its interval (+-SLACK bases; each truth record claimed once, calls by decreasing p_somatic);
   otherwise FP. alt_is_grch38_base calls (the read shows GRCh38's own base) are always FP. In the BED evaluation
   only intervals overlapping the BED count. The PoN cannot tag class-ii or class-iii calls (no GRCh38 allele): they
   stay in the post-PoN call set.
3. Class-iii calls (no anchor): counted as FP (worst case) or dropped.
4. metrics.json: for bed / nobed x raw / pon:
   placed            the run's own eval (vcf_chr1/eval_*/report.json): rtg on the placed calls only
   placed_unplaced_fp  placed + every unplaced PASS call as FP (the report's "precision with them as false calls")
   projected_rtg     rtg on placed + class i; class ii and iii dropped
   projected_ii      + class ii by interval; class iii dropped
   projected_worst   + class ii by interval; class iii as FP (no position: counted in the BED evaluation too)
   and the unplaced fractions of the PASS calls (whole chr1, no BED) before and after the projection and the PoN.
"""
import argparse
import bisect
import gzip
import json
import subprocess
import sys
import time
from pathlib import Path

REPO = Path(__file__).resolve().parents[4]
PON_DIR = Path("/scratch/jshen/data/Pansoma/panel_of_normal_VCFs")
PONS = [PON_DIR / n for n in ("af-only-gnomad.hg38.vcf.gz", "Homo_sapiens_assembly38.dbsnp138.vcf.gz",
                              "1000g_pon.hg38.vcf.gz", "CoLoRSdb.GRCh38.v1.1.0.deepvariant.glnexus.vcf.gz")]
SDF = "/scratch/jshen/REF/GRCh38.primary.sdf"
RTG = "/scratch/jshen/bin/rtg-tools-3.12.1/rtg"
PY = "/wanglab/jshen/anaconda3/bin/python"
SLACK = 1


def read_tsv(path):
    with open(path) as f:
        cols = f.readline().rstrip("\n").split("\t")
        return [dict(zip(cols, line.rstrip("\n").split("\t"))) for line in f]


def read_bed(path, chroms):
    iv = {}
    opener = gzip.open if str(path).endswith(".gz") else open
    with opener(path, "rt") as f:
        for line in f:
            if line.startswith(("#", "track", "browser")) or not line.strip():
                continue
            c, s, e = line.split("\t")[:3]
            if c in chroms:
                iv.setdefault(c, []).append((int(s), int(e)))
    merged = {}
    for c, xs in iv.items():
        xs.sort()
        out = []
        for s, e in xs:
            if out and s <= out[-1][1]:
                out[-1][1] = max(out[-1][1], e)
            else:
                out.append([s, e])
        merged[c] = ([s for s, _ in out], [e for _, e in out])
    return merged


def overlaps(bed, chrom, s0, e0):
    if bed is None:
        return True
    if chrom not in bed:
        return False
    starts, ends = bed[chrom]
    e0 = max(e0, s0 + 1)
    k = bisect.bisect_right(starts, e0 - 1) - 1
    return k >= 0 and ends[k] > s0


def fn_records(vcfeval_dir):
    out = []
    with gzip.open(Path(vcfeval_dir) / "fn.vcf.gz", "rt") as f:
        for line in f:
            if line.startswith("#"):
                continue
            c, pos, _, ref, alt = line.split("\t")[:5]
            out.append((c, int(pos), ref, alt))
    return out


def match_intervals(calls, fns, kind):
    """Greedy: class-ii calls (dicts) by decreasing p_somatic claim one FN truth record inside their interval."""
    by_chrom = {}
    for c, pos, ref, alt in fns:
        # truth span, 0-based: SNV the base; INDEL the padding base + the changed bases
        s0 = pos - 1
        e0 = pos - 1 + max(len(ref), 1) + (0 if kind == "SNV" else 1)
        by_chrom.setdefault(c, []).append([s0, e0, False])
    for v in by_chrom.values():
        v.sort()
    tp = 0
    for call in sorted(calls, key=lambda c: -float(c["p_somatic"] or 0)):
        c, s0, e0 = call["grch38_chrom"], int(call["grch38_start0"]), int(call["grch38_end0"])
        lo, hi = s0 - SLACK, max(e0, s0 + 1) + SLACK
        for rec in by_chrom.get(c, []):
            if rec[0] >= hi:
                break
            if rec[1] > lo and not rec[2]:
                rec[2] = True
                call["_tp"] = True
                tp += 1
                break
    return tp


def projected_status(vcfeval_dir):
    """PASS class-i projections (INFO PROJ) among rtg's TP and FP call records."""
    out = {}
    for name in ("tp", "fp"):
        n = 0
        with gzip.open(Path(vcfeval_dir) / f"{name}.vcf.gz", "rt") as f:
            for line in f:
                if not line.startswith("#") and "PROJ=" in line.split("\t")[7]:
                    n += 1
        out[name] = n
    return out


def prf(tp, fp, fn):
    p = tp / (tp + fp) if tp + fp else None
    r = tp / (tp + fn) if tp + fn else None
    f = 2 * p * r / (p + r) if p and r else (0.0 if p is not None and r is not None else None)
    return dict(tp=tp, fp=fp, fn=fn, precision=p, recall=r, f1=f)


def at(report, label):
    a = report[label]["at_threshold"]
    return a["tp_baseline"], a["fp"], a["fn"], a["tp_call"]


def run_vcfeval(calls, truth, bed, out, threads, pon_vcf=None):
    if (out / "report.json").exists():
        return json.loads((out / "report.json").read_text())
    if out.exists():  # an interrupted earlier attempt: keep it aside (vcfeval.py needs a new output)
        out.rename(out.with_name(f"{out.name}.incomplete.{int(time.time())}"))
    cmd = [PY, "-m", "pansoma_net_v2.vcfeval", "--calls", str(calls), "--truth", truth, "--sdf", SDF, "--rtg", RTG,
           "--threads", str(threads), "--output", str(out)]
    cmd += ["--pon-vcf", str(pon_vcf)] if pon_vcf else ["--pon", *map(str, PONS)]
    if bed:
        cmd += ["--bed", bed]
    print(" ".join(cmd), flush=True)
    subprocess.run(cmd, check=True, cwd=REPO / "machine_learning")
    return json.loads((out / "report.json").read_text())


def count_pass(vcf, pon=False):
    n = 0
    with gzip.open(vcf, "rt") as f:
        for line in f:
            if not line.startswith("#") and line.split("\t")[6] == "PASS":
                n += 1
    return n


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("--vcf-dir", required=True)
    ap.add_argument("--projection", required=True)
    ap.add_argument("--threads", type=int, default=4)
    args = ap.parse_args()
    vdir, pdir = Path(args.vcf_dir), Path(args.projection)
    pj = json.loads((pdir / "projection.json").read_text())
    kind = pj["kind"]
    own = {b: json.loads((vdir / f"eval_{b}" / "report.json").read_text()) for b in ("bed", "nobed")}
    truth, bed = own["bed"]["truth"], own["bed"]["bed"]
    calls = pdir / "complete.linear.vcf.gz"
    rep = {"bed": run_vcfeval(calls, truth, bed, pdir / "eval_bed", args.threads)}
    rep["nobed"] = run_vcfeval(calls, truth, None, pdir / "eval_nobed", args.threads,
                               pon_vcf=pdir / "eval_bed" / "complete.linear.pon.vcf.gz")
    rows = read_tsv(pdir / "projected.tsv")
    passing = [r for r in rows if r["filter"] == "PASS"]
    ii = [r for r in passing if r["proj_class"] == "ii"]
    iii = [r for r in passing if r["proj_class"] == "iii"]
    chroms = set(own["bed"]["chromosomes"])
    bedmap = read_bed(bed, chroms) if bed else None
    metrics = dict(run=pj["run"], set=pj["set"], kind=kind, graph=pj["graph"], vcf_dir=str(vdir.resolve()),
                   threshold=own["bed"]["threshold"], unplaced_pass=pj["unplaced_pass"],
                   classes_pass=pj["classes_pass"], details_pass=pj["details_pass"], chm13_pass=pj["chm13_pass"],
                   chm13_by_class_pass=pj["chm13_by_class_pass"], merged_with_placed=pj["merged_with_placed"],
                   eval={})
    for b in ("bed", "nobed"):
        u_all = own[b].get("unplaced_pass") or 0
        ii_b = [dict(r) for r in ii if overlaps(bedmap if b == "bed" else None, r["grch38_chrom"],
                                                int(r["grch38_start0"]), int(r["grch38_end0"]))]
        for label in ("raw", "pon"):
            tp0, fp0, fn0, _ = at(own[b], label)
            tp1, fp1, fn1, _ = at(rep[b], label)
            fns = fn_records(pdir / f"eval_{b}" / "vcfeval" / f"{label}_calls")
            calls_ii = [dict(r) for r in ii_b if r["proj_detail"] != "alt_is_grch38_base"]
            tp_ii = match_intervals(calls_ii, fns, kind)
            fp_ii = len(ii_b) - tp_ii  # alt_is_grch38_base: the read shows GRCh38's base, never a GRCh38 variant
            metrics["eval"][f"{b}_{label}"] = dict(
                placed=prf(tp0, fp0, fn0),
                placed_unplaced_fp=prf(tp0, fp0 + u_all, fn0),
                projected_rtg=prf(tp1, fp1, fn1),
                projected_ii=prf(tp1 + tp_ii, fp1 + fp_ii, fn1 - tp_ii),
                projected_worst=prf(tp1 + tp_ii, fp1 + fp_ii + len(iii), fn1 - tp_ii),
                class_ii_calls=len(ii_b), class_ii_alt_is_grch38=len(ii_b) - len(calls_ii), class_ii_tp=tp_ii, class_iii_calls=len(iii),
                class_i_rtg=projected_status(pdir / f"eval_{b}" / "vcfeval" / f"{label}_calls"),
                recall_gain_truth=tp1 - tp0)
    # PASS counts on chr1 (no BED) for the unplaced fractions
    pon_own = Path(own["bed"]["pon_vcf"])
    pon_new = pdir / "eval_bed" / "complete.linear.pon.vcf.gz"
    placed_pass_raw = count_pass(vdir / Path(own["bed"]["calls"]).name)
    placed_pass_pon = count_pass(pon_own)
    complete_pass_raw = count_pass(calls)
    complete_pass_pon = count_pass(pon_new)
    u, u_left = pj["unplaced_pass"], len(ii) + len(iii)
    metrics["pass_counts"] = dict(
        placed_raw=placed_pass_raw, placed_pon=placed_pass_pon, complete_raw=complete_pass_raw,
        complete_pon=complete_pass_pon, unplaced_before=u, unplaced_after=u_left,
        projected_class_i_pass=pj["classes_pass"].get("i", 0),
        frac_unplaced_raw_before=u / (placed_pass_raw + u) if placed_pass_raw + u else None,
        frac_unplaced_pon_before=u / (placed_pass_pon + u) if placed_pass_pon + u else None,
        frac_unplaced_raw_after=u_left / (complete_pass_raw + u_left) if complete_pass_raw + u_left else None,
        frac_unplaced_pon_after=u_left / (complete_pass_pon + u_left) if complete_pass_pon + u_left else None)
    # class-i PASS projections tagged by the PoN
    tagged = 0
    with gzip.open(pon_new, "rt") as f:
        for line in f:
            if line.startswith("#"):
                continue
            x = line.split("\t")
            if "PROJ=" in x[7] and "PANSOMA_PON" in x[7] and x[6] == "PanelOfNormals":
                tagged += 1
    metrics["pass_counts"]["projected_class_i_pass_pon_tagged"] = tagged
    (pdir / "metrics.json").write_text(json.dumps(metrics, indent=2) + "\n")
    print(json.dumps(metrics, indent=2))


if __name__ == "__main__":
    main()
