"""One-time GAM preparation: the sorted GAM and its GAI that discover, build and orchestrate read.

    python -m indexed_gam_pipeline_v4.tools.gam_prep sort --gam IN.gam --output IN.sorted.gam [--threads 8] [--tmp-dir DIR] [--vg VG]
    python -m indexed_gam_pipeline_v4.tools.gam_prep check --gam IN.sorted.gam [--index GAI] [--input IN.gam] [--threads 8] [--vg VG]

sort runs `vg gamsort -i` (vg: --vg, else $PANSOMA_VG, else vg on PATH; its temporary files go to --tmp-dir, else
$TMPDIR), checks the result and publishes <output> and <output>.gai atomically. check opens the GAI with the
pipeline's reader (the 'GAI!' magic, format number 1, offsets inside the GAM) and reads the first records to see
that they are GAM records in node order; with --input it also compares `vg stats -a` of the sorted GAM with the
unsorted input (the numbers of alignments, primary, secondary, aligned and perfect records and of matched bases,
which sorting keeps). Nothing here is imported at run time.
"""
import argparse
import json
import os
from pathlib import Path
import subprocess
import time

import pysam

from ..common import stamp
from ..gam_reader import IndexedGam, decode, group
from .graph_prep import vg_command, vg_version

KEPT_COUNTS = ("Total alignments", "Total primary", "Total secondary", "Total aligned", "Total perfect", "Matches")


def first_records(gam, limit):
    """(records read, whether their smallest node IDs never decrease) over the first `limit` GAM records."""
    seen, ordered, last = 0, True, 0
    with pysam.BGZFile(str(gam), "rb") as stream:
        while seen < limit:
            messages = group(stream)
            if messages is None:
                break
            for raw in messages:
                nodes = [m.position.node_id for m in decode(raw).path.mapping if m.position.node_id > 0]
                seen += 1
                if nodes:
                    ordered &= min(nodes) >= last
                    last = max(last, min(nodes))
    return seen, ordered


def stats_process(gam, vg, threads):
    """`vg stats -a` of a GAM, started (the two GAMs of a check run side by side)."""
    return subprocess.Popen([vg, "stats", "-p", str(threads), "-a", str(gam)], stdout=subprocess.PIPE, text=True)


def parse_counts(text):
    counts = {}
    for line in text.splitlines():
        key, _, value = line.partition(":")
        if key.strip() in KEPT_COUNTS:
            counts[key.strip()] = int(value.split()[0])
    return counts


def check(gam, index=None, unsorted=None, vg=None, threads=1, records=10000):
    """Raises unless the GAI opens with IndexedGam and the first records are GAM records in node order (and, with
    `unsorted`, the kept `vg stats -a` counts equal)."""
    started = time.perf_counter()
    reader = IndexedGam(gam, index)
    seen, ordered = first_records(gam, records)
    report = dict(gam=stamp(gam), index=stamp(reader.index), gai_version=reader.version, bins=len(reader.bins),
                  first_records=seen, first_records_in_node_order=ordered)
    if not reader.bins or not seen or not ordered:
        raise ValueError(f"Not a sorted, indexed GAM: {report}")
    if unsorted:
        vg = vg_command(vg)
        processes = [stats_process(path, vg, threads) for path in (gam, unsorted)]
        texts = [p.communicate()[0] for p in processes]
        if any(p.returncode for p in processes):
            raise RuntimeError("vg stats -a failed")
        sorted_counts, input_counts = (parse_counts(t) for t in texts)
        report.update(input=stamp(unsorted), counts=sorted_counts, input_counts=input_counts, vg=vg_version(vg))
        if sorted_counts != input_counts or set(sorted_counts) != set(KEPT_COUNTS):
            raise ValueError(f"Sorting changed the alignment counts: {sorted_counts} vs {input_counts}")
    report["seconds"] = time.perf_counter() - started
    return report


def sort(gam, output, vg=None, threads=1, tmp_dir=None):
    """vg gamsort -i into <output>.tmp and <output>.gai.tmp, check, then publish (the GAI first, the GAM last)."""
    vg, output, started = vg_command(vg), Path(output), time.perf_counter()
    index = Path(str(output) + ".gai")
    for path in (output, index):
        if path.exists():
            raise ValueError(f"Output exists: {path}")
    work, work_index = output.with_name(output.name + ".tmp"), index.with_name(index.name + ".tmp")
    env = dict(os.environ)
    if tmp_dir:
        Path(tmp_dir).mkdir(parents=True, exist_ok=True)
        env["TMPDIR"] = str(tmp_dir)
    with work.open("wb") as stream:
        subprocess.run([vg, "gamsort", "-t", str(threads), "-p", "-i", str(work_index), str(gam)], stdout=stream,
                       check=True, env=env)
    report = check(work, work_index)
    work_index.rename(index)
    work.rename(output)
    return dict(input=stamp(gam), output=stamp(output), index=stamp(index), vg=vg_version(vg), check=report,
                seconds=time.perf_counter() - started)


def main(argv=None):
    parser = argparse.ArgumentParser(prog="python -m indexed_gam_pipeline_v4.tools.gam_prep",
                                     description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    commands = parser.add_subparsers(dest="command", required=True)
    p = commands.add_parser("sort", help="vg gamsort -i, check, publish <output> and <output>.gai")
    p.add_argument("--gam", required=True, help="unsorted GAM (giraffe output)")
    p.add_argument("--output", required=True, help="new sorted GAM; its GAI is <output>.gai")
    p.add_argument("--threads", type=int, default=1)
    p.add_argument("--tmp-dir", help="vg gamsort's temporary files (default: $TMPDIR)")
    p.add_argument("--vg", help="vg executable (default: $PANSOMA_VG, else vg on PATH)")
    p = commands.add_parser("check", help="GAI readable by the pipeline, records in node order; --input: same counts")
    p.add_argument("--gam", required=True)
    p.add_argument("--index", help="default: GAM path + .gai")
    p.add_argument("--input", help="the unsorted GAM: compare vg stats -a counts (reads both GAMs)")
    p.add_argument("--threads", type=int, default=1)
    p.add_argument("--vg", help="vg executable (default: $PANSOMA_VG, else vg on PATH)")
    args = parser.parse_args(argv)
    if args.command == "sort":
        result = sort(args.gam, args.output, args.vg, args.threads, args.tmp_dir)
    else:
        result = check(args.gam, args.index, args.input, args.vg, args.threads)
    print(json.dumps(result, indent=2))


if __name__ == "__main__":
    main()
