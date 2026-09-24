#!/usr/bin/env python3
"""Render the v6 counterparts of the 100 v5 examples in ../hg008_pacbio_v5/.

Every batch holding a v5 example was rebuilt with the v6 code (--debug-rows) into
REBUILT/<task>_<batch>/{SNV,INDEL}. Each v5 example is matched to the v6 site that now
holds it: the same candidate id among a site's alleles (position unchanged), else, for an
indel moved by left-normalization, the nearest site to its left on the same node, else on a
node at most 10 IDs away in either direction (the forward-strand left of a node need not have
a lower ID), with an allele of the same type and length (the normalized allele; an
insertion may be rotated), else the nearest site of the same kind on that node. Sites of
the supplement run (indels normalized off the target nodes) are searched too. index.tsv records the match, the v5
PNG, the site's alleles and the strand of the A1 rows (reverse fraction; one-strand = all
of >= 3 A1 rows on one strand, read from the strand channel). Summaries are read with their
--debug-rows audit cut out as text, and each example is rendered from a one-tensor copy.
Run with the matplotlib-capable interpreter.
"""
from concurrent.futures import ThreadPoolExecutor
import csv
import json
from pathlib import Path
import subprocess
import sys

import numpy as np

HERE = Path(__file__).resolve().parent
REPO = HERE.parents[2]
V5 = HERE.parent / "hg008_pacbio_v5"
REBUILT = Path(sys.argv[1] if len(sys.argv) > 1 else REPO / "tmp/v2_speedup_20260923/v6c/main")
SUPPLEMENT = [Path(p) for p in sys.argv[2:]]  # supplement-run outputs (each with SNV/ and INDEL/)
TIMING = Path("/scratch/jshen/data/HG008_GIAB/pansoma_v2_tensors/Liss_lab_PacBio_Revio_20240125/v5_tensors/shared")


def batch_of(task, node):
    for line in open(TIMING / f"task_{task:04d}/batch_timing.ndjson"):
        b = json.loads(line)
        if b["first_node"] <= node <= b["last_node"]:
            return b["batch"]
    raise KeyError((task, node))


def match(example, sites):
    cid = example["candidate_id"]
    node, start, kind = int(cid.split(":")[0]), int(cid.split(":")[1]), cid.split(":")[2]
    for m in sites:
        if any(a["candidate_id"] == cid for a in m["alleles"]):
            return m, "same position"
    ref, alt = cid.split(":", 3)[3].split(">")
    length = len(alt) if kind == "INS" else len(ref)
    moved = [m for m in sites if (m["node_id"] == node and m["start"] <= start or 0 < abs(m["node_id"] - node) <= 10)
             and any(a["event_type"] == kind and a["event_length"] == length for a in m["alleles"])]
    if moved:  # the nearest such site: on this node to the left, else on a neighbouring (short) node
        m = min(moved, key=lambda m: (abs(m["node_id"] - node), -m["start"]))
        return m, "left-normalized" + (" to another node" if m["node_id"] != node else "")
    near = [m for m in sites if m["node_id"] == node]
    return (min(near, key=lambda m: abs(m["start"] - start)), "nearest site") if near else (None, "none")


def light(line):
    """A summary record without its --debug-rows audit (cut out as text: records are huge)."""
    i = line.find('"rows": [')
    if i >= 0:
        j = line.index('"sample_unit": ', i)  # the audit fields sit between "rows" and the site fields
        line = line[:i] + line[j:]
    return json.loads(line)


def a1_strand(folder, m):
    """(rows, reverse fraction, one-strand) of the A1 block, read from the strand channel."""
    x = np.load(folder / f"shard_{m['shard_index']:05d}_data.npy", mmap_mode="r")[m["index_within_shard"]]
    rows = [r for g in m["row_groups"] if g["allele"] == "A1" for r in range(g["start_row"], g["end_row"])]
    reverse = [int(np.asarray(x[7, r]).max()) == 2 for r in rows]
    fraction = sum(reverse) / len(reverse) if reverse else None
    return len(reverse), fraction, bool(len(reverse) >= 3 and fraction in (0.0, 1.0))


def stage(folder, m, directory):
    """One-tensor shard + its record + the manifest, so the renderer never parses multi-GB summaries."""
    directory.mkdir(parents=True, exist_ok=True)
    x = np.load(folder / f"shard_{m['shard_index']:05d}_data.npy", mmap_mode="r")[m["index_within_shard"]]
    np.save(directory / "shard_00000_data.npy", np.asarray(x)[None])
    (directory / "manifest.json").write_text((folder / "manifest.json").read_text())
    (directory / "variant_summary.ndjson").write_text(json.dumps(dict(m, shard_index=0, index_within_shard=0)) + "\n")
    return directory / "shard_00000_data.npy"


def render(job):
    shard, png = job
    subprocess.run([sys.executable, str(REPO / "scripts/visualize_tensor.py"), str(shard), "-i", "0", "-o", str(png)],
                   check=True, stdout=subprocess.DEVNULL)


def main():
    examples = list(csv.DictReader(open(V5 / "index.tsv"), delimiter="\t"))
    header = ["index", "png", "v5_png", "match", "site_id", "alleles", "site_coverage", "site_counts",
              "a1_rows", "a1_reverse_fraction", "a1_one_strand", "v5_candidate_id", "v5_af", "v5_alt_reverse_fraction"]
    folders, wanted = {}, []
    for ex in examples:
        task, node = int(ex["task"]), int(ex["candidate_id"].split(":")[0])
        kind = "SNV" if ex["candidate_id"].split(":")[2] == "SNP" else "INDEL"
        folder = REBUILT / f"{task}_{batch_of(task, node)}" / kind
        folders.setdefault(folder, None)
        wanted.append((ex, folder))
    supplement = {}
    for root in SUPPLEMENT:
        for kind in ("SNV", "INDEL"):
            with open(root / kind / "variant_summary.ndjson") as stream:
                supplement.setdefault(kind, []).extend((root / kind, light(line)) for line in stream)
    for folder in folders:  # stream each summary once, dropping the debug rows record by record
        with open(folder / "variant_summary.ndjson") as stream:
            # a site of the batch, or of the supplement run (indels normalized off the targets)
            folders[folder] = [(folder, light(line)) for line in stream] + supplement.get(folder.name, [])
    rows, jobs = ["\t".join(header)], []
    stage_root = REBUILT.parent / "v6_render_stage"
    for ex, folder in wanted:
        owner = {id(m): f for f, m in folders[folder]}
        m, how = match(ex, [m for _, m in folders[folder]])
        folder = owner[id(m)] if m is not None else folder
        if m is not None and folder.parent in SUPPLEMENT:
            how += " (supplement run)"
        if m is None:
            rows.append("\t".join([ex["index"], "", ex["png"], how] + [""] * (len(header) - 4)))
            continue
        name = ex["png"].replace(".png", "") + "_v6.png"
        jobs.append((stage(folder, m, stage_root / name.replace(".png", "")), HERE / name))
        a1, fraction, one = a1_strand(folder, m)
        alleles = " ".join(f"{a['label']}:{a['event_type']}:{a['ref'] or '-'}>{a['alt'] or '-'}:af={a['af']:.3f}"
                           for a in m["alleles"])
        rows.append("\t".join(map(str, [ex["index"], name, ex["png"], how, m["site_id"], alleles, m["site_coverage"],
                                        json.dumps(m["site_counts"]), a1,
                                        "" if fraction is None else round(fraction, 3), one, ex["candidate_id"],
                                        ex["af"], ex["alt_reverse_fraction"]])))
    with ThreadPoolExecutor(8) as pool:
        list(pool.map(render, jobs))
    (HERE / "index.tsv").write_text("\n".join(rows) + "\n")
    print(f"rendered {len(jobs)} of {len(rows) - 1} examples into {HERE}")


if __name__ == "__main__":
    main()
