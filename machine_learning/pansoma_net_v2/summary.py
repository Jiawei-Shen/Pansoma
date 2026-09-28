"""Per-epoch table of training runs, their current step, and their test results.

    cd machine_learning
    python -m pansoma_net_v2.summary HG008_Illumina_SNV_base HG008_Illumina_INDEL_base ...   # names under RUNS
    python -m pansoma_net_v2.summary /path/to/run_dir

RUNS is $PANSOMA_RUNS or /scratch/jshen/data/pansoma_net_v2_runs. Columns: training loss, speed and minutes
of the epoch, GPU peak (allocated GiB); validation per tensor: somatic AP, F1 / precision / recall at the
best-F1 threshold t (s.*), t itself; against the truth VCF at the same t (t.*, metrics.truth_report: every
truth allele once, truth alleles without a tensor are misses, found truth of the other kind are true calls)
and its ceiling; germline argmax F1; "*" marks the best epoch so far by the run's --select. combine scores a
sample's SNV and INDEL models together.
"""
import argparse
import glob
import json
import os
from pathlib import Path

RUNS = Path(os.environ.get("PANSOMA_RUNS", "/scratch/jshen/data/pansoma_net_v2_runs"))


def run_dir(name):
    p = Path(name)
    return p if p.is_dir() else RUNS / name


def current_step(d):
    """The last 'epoch N step S/T' line of the run's newest Slurm log (jobs/v2_<name>-<job>.out), if any."""
    logs = sorted(glob.glob(str(d.parent / "jobs" / f"v2_{d.name}-*.out")), key=os.path.getmtime)
    if not logs:
        return None
    last = None
    with open(logs[-1], errors="replace") as f:
        for line in f:
            if line.startswith("epoch ") and " step " in line:
                last = line.strip()
    return last


def show(d):
    print(f"== {d.name}  ({d})")
    rows = [json.loads(line) for line in open(d / "metrics.jsonl")] if (d / "metrics.jsonl").exists() else []
    if rows:
        print(f"{'epoch':>5} {'loss':>7} {'t/s':>6} {'min':>5} {'GPU':>5} | {'s.AP':>6} {'s.F1':>6} {'s.P':>6} {'s.R':>6} "
              f"{'t':>6} | {'t.AP':>6} {'t.F1':>6} {'t.P':>6} {'t.R':>6} {'ceil':>6} | {'g.F1':>6}")
        best, select = -1.0, _select(d)
        for r in rows:
            v = r["val"]
            s = v["thresholded"]["somatic"]
            tr = v.get("truth")
            score = (tr["best"]["f1"] if select == "truth_f1" and tr else v["somatic_ap"] if select == "ap" else s["f1"])
            mark = "*" if score > best else " "
            best = max(best, score)
            gpu = (r.get("gpu_peak_gib") or {}).get("allocated")
            truth = (f"{tr['ap']:>6.3f} {tr['at_threshold']['f1']:>6.3f} {tr['at_threshold']['precision']:>6.3f} "
                     f"{tr['at_threshold']['recall']:>6.3f} {tr['ceiling']:>6.3f}") if tr else f"{'-':>6} " * 4 + f"{'-':>6}"
            print(f"{r['epoch']:>5} {r['train_loss']:>7.4f} {r['tensors_per_second']:>6.0f} {r['train_seconds'] / 60:>5.1f} "
                  f"{gpu if gpu is not None else '-':>5} | {v['somatic_ap']:>6.3f} {s['f1']:>6.3f} {s['precision']:>6.3f} "
                  f"{s['recall']:>6.3f} {v['threshold']:>6.3f} | {truth} | {v['argmax']['germline']['f1']:>6.3f} {mark}")
    else:
        print("  no finished epoch yet")
    step = current_step(d)
    if step:
        print(f"  now: {step}")
    for m in sorted((d / "test_chr1").glob("*.metrics.json")) if (d / "test_chr1").exists() else []:
        t = json.loads(m.read_text())
        s = t["thresholded"]["somatic"]
        print(f"  test {m.name.replace('.metrics.json', '')}: {t['tensors']:,} tensors ({t['off_reference_in_test']:,} "
              f"off-reference), t {t['threshold']:.3f}")
        print(f"    per tensor:    somatic AP {t['somatic_ap']:.3f} F1 {s['f1']:.3f} P {s['precision']:.3f} "
              f"R {s['recall']:.3f} (support {s['support']:,}) | germline F1 {t['argmax']['germline']['f1']:.3f}")
        tr = t.get("truth")
        if tr:
            at = tr["at_threshold"]
            print(f"    vs truth VCF:  somatic AP {tr['ap']:.3f} F1 {at['f1']:.3f} P {at['precision']:.3f} "
                  f"R {at['recall']:.3f} ({at['tp']:,} of {tr['truth_alleles']:,} truth found, {at['fp']:,} false calls; "
                  f"ceiling {tr['ceiling']:.3f})")
            if "other_kind_tensors" in tr:
                print(f"                   + {at['other_tp']:,} truth of the other kind found (true calls, each once; "
                      f"{tr['other_kind_tensors']:,} tensors of {tr['other_kind_truth_with_tensor']:,} such truth)")

def _select(d):
    """The run's --select, from the args.json that train writes at its start (f1 for older runs)."""
    try:
        return json.loads((d / "args.json").read_text()).get("select", "f1")
    except (OSError, ValueError):
        return "f1"


def main(argv=None):
    p = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    p.add_argument("runs", nargs="*", help="run names under RUNS or run directories (default: every run in RUNS)")
    args = p.parse_args(argv)
    names = args.runs or sorted(x.name for x in RUNS.iterdir() if (x / "train.log").exists() and "attempt" not in x.name)
    for name in names:
        show(run_dir(name))
        print()


if __name__ == "__main__":
    main()
