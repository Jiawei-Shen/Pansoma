"""Per-epoch table of training runs, their current step, and their test results.

    cd machine_learning
    python -m pansoma_net_v2.summary HG008_Illumina_SNV_base HG008_Illumina_INDEL_base ...   # names under RUNS
    python -m pansoma_net_v2.summary /path/to/run_dir

RUNS is $PANSOMA_RUNS or /scratch/jshen/data/pansoma_net_v2_runs. Columns: training loss, speed and minutes
of the epoch, GPU peak (allocated GiB); validation somatic AP, F1 / precision / recall at the best-F1
threshold t, t itself, and the germline argmax F1; "*" marks the best epoch so far.
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
              f"{'t':>6} | {'g.F1':>6}")
        best = -1.0
        for r in rows:
            v = r["val"]
            s = v["thresholded"]["somatic"]
            mark = "*" if s["f1"] > best else " "
            best = max(best, s["f1"])
            gpu = (r.get("gpu_peak_gib") or {}).get("allocated")
            print(f"{r['epoch']:>5} {r['train_loss']:>7.4f} {r['tensors_per_second']:>6.0f} {r['train_seconds'] / 60:>5.1f} "
                  f"{gpu if gpu is not None else '-':>5} | {v['somatic_ap']:>6.3f} {s['f1']:>6.3f} {s['precision']:>6.3f} "
                  f"{s['recall']:>6.3f} {v['threshold']:>6.3f} | {v['argmax']['germline']['f1']:>6.3f} {mark}")
    else:
        print("  no finished epoch yet")
    step = current_step(d)
    if step:
        print(f"  now: {step}")
    for m in sorted((d / "test_chr1").glob("*.metrics.json")) if (d / "test_chr1").exists() else []:
        t = json.loads(m.read_text())
        s = t["thresholded"]["somatic"]
        print(f"  test {m.name.replace('.metrics.json', '')}: {t['tensors']:,} tensors ({t['off_reference_in_test']:,} "
              f"off-reference) | somatic AP {t['somatic_ap']:.3f} F1 {s['f1']:.3f} P {s['precision']:.3f} "
              f"R {s['recall']:.3f} (support {s['support']:,}, t {t['threshold']:.3f}) | germline F1 "
              f"{t['argmax']['germline']['f1']:.3f}")


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
