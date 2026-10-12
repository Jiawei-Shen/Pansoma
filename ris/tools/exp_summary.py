"""Table of the ris INDEL experiments (exp_submit.sh runs): per run and set, the decision tuned on chr2 and the chr1
results (in the BED) with the default decision and with the tuned one. Standard library only (login node).

    python3 ris/tools/exp_summary.py [RUN_PREFIX ...]      default prefixes: ix
"""
import glob
import json
import os
import re
import sys

RUNS = "/storage2/fs1/epigenome/Active/jiawei.s/pansoma_net_v2_runs"


def row(report, key="pon"):
    if not os.path.exists(report):
        return None
    m = re.search(rf"^{key}\s+(.*)$", open(report).read(), re.M)
    if not m:
        return None
    v = m.group(1).split()
    return dict(calls=v[0], tp=v[1], p=float(v[4]), r=float(v[5]), f1=float(v[6]), best=float(v[8]), ceiling=float(v[-1]))


def main():
    prefixes = sys.argv[1:] or ["ix"]
    print(f"{'run':24s} {'set':18s} {'ep':>3s} | {'chr1 default':19s} | {'GENERAL thr (chr2 F1)':21s} {'chr1 GENERAL P/R/F1':19s} "
          f"{'calls':>6s} {'noBED F1':>8s} | {'per-set thr,floor,AF':20s} {'chr1 tuned P/R/F1':19s}")
    for d in sorted(glob.glob(f"{RUNS}/*")):
        name = os.path.basename(d)
        if not any(name.startswith(p) for p in prefixes) or not os.path.isdir(d):
            continue
        tl = open(f"{d}/train.log").read() if os.path.exists(f"{d}/train.log") else ""
        b = re.search(r"best epoch (\d+) \(f1 ([\d.]+)", tl)
        best = f"{b.group(1)}" if b else ("..." if tl else "-")
        sets = sorted({re.sub(r"^vcf_|_chr[12](_tuned|_general)?$", "", os.path.basename(x)) for x in glob.glob(f"{d}/vcf_*")}) or ["-"]
        for s in sets:
            t = f"{d}/tune_{s}.json"
            tune = json.load(open(t)) if os.path.exists(t) else None
            g = f"{d}/tune_{s}.general.json"
            gt = json.load(open(g)) if os.path.exists(g) else None
            gtt = f"{gt['threshold']:.3f} ({gt['tuning_f1']:.3f})" if gt else "-"
            tt = f"{tune['threshold']:.3f},{tune['pon_floor']},{tune['read_af_min']}" if tune else "-"
            dflt = row(f"{d}/vcf_{s}_chr1/eval_bed/report.txt")
            gen = row(f"{d}/vcf_{s}_chr1_general/eval_bed/report.txt")
            gnb = row(f"{d}/vcf_{s}_chr1_general/eval_nobed/report.txt")
            tun = row(f"{d}/vcf_{s}_chr1_tuned/eval_bed/report.txt")
            fmt = lambda x: f"{x['p']:.3f}/{x['r']:.3f}/{x['f1']:.3f}" if x else "-"
            nob = f"{gnb['f1']:.3f}" if gnb else "-"
            print(f"{name:24s} {s:18s} {best:>3s} | {fmt(dflt):19s} | {gtt:21s} {fmt(gen):19s} "
                  f"{(gen or {}).get('calls', '-'):>6s} {nob:>8s} | {tt:20s} {fmt(tun):19s}")


if __name__ == "__main__":
    main()
