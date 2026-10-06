"""Sum the per-task build summaries (task_*.log JSON line: timing, tensors, filtered candidates) and /usr/bin/time
figures (task_*.resources.txt) of one run root: where the tensor build time goes, per graph."""
import json, re, sys
from pathlib import Path
from collections import Counter
def summarise(root):
    root = Path(root); t = Counter(); n = 0; extra = Counter()
    for log in sorted((root / "logs").glob("task_*.log")):
        res = log.with_name(log.stem + ".resources.txt")
        lines = [l for l in log.read_text().splitlines() if l.startswith('{"tensors"')]
        if not lines or not res.exists():
            continue
        d = json.loads(lines[-1]); n += 1
        for k, v in d.get("timing", {}).items():
            t[k] += v
        extra["tensors"] += d["tensors"]; extra["filtered"] += d.get("filtered_candidates", 0)
        for m in re.finditer(r"Batch \d+: (\d+) target nodes, (\d+) alignments, (\d+) context nodes, (\d+) candidates", log.read_text()):
            extra["targets"] += int(m[1]); extra["alignments"] += int(m[2]); extra["context"] += int(m[3]); extra["candidates"] += int(m[4])
        r = res.read_text()
        g = lambda pat: float(re.search(pat, r).group(1))
        extra["user"] += g(r"User time \(seconds\): ([\d.]+)"); extra["sys"] += g(r"System time \(seconds\): ([\d.]+)")
        w = re.search(r"Elapsed \(wall clock\) time \(h:mm:ss or m:ss\): ([\d:.]+)", r).group(1).split(":")
        extra["wall"] += sum(float(x) * 60 ** i for i, x in enumerate(reversed(w)))
        extra["fs_in_blocks"] += g(r"File system inputs: (\d+)")
    return n, t, extra
for root in sys.argv[1:]:
    n, t, e = summarise(root)
    name = Path(root).parent.name
    h = lambda s: s / 3600
    keys = sorted(t, key=lambda k: -t[k])
    print(json.dumps(dict(set=name, tasks=n, wall_h=round(h(e["wall"]), 1), user_h=round(h(e["user"]), 1), sys_h=round(h(e["sys"]), 1),
                          cpu_pct=round(100 * (e["user"] + e["sys"]) / e["wall"], 1) if e["wall"] else None,
                          timing_h={k: round(h(t[k]), 1) for k in keys},
                          tensors=e["tensors"], targets=e["targets"], alignments=e["alignments"], context_nodes=e["context"],
                          candidates=e["candidates"], filtered=e["filtered"], fs_in_GB=round(e["fs_in_blocks"] * 512 / 1e9, 1))))
