"""E9a/E9b graph preparation: the edges of a graph as a CSR over oriented node sides (from the L lines of its GFA,
the GFA of the GBZ whose segment names are the GBZ node IDs) and the CHM13 path coordinates (vg paths -A on the GBZ) in
the reference-path format of indexed_gam_pipeline_v4 (tensor_postprocessing.reference_path.ReferencePath reads it).

    python prep_graph.py --gfa G.gfa --l-range START END --gbz G.gbz --vg VG \
        --grch38-path <graph_index>/<graph>.grch38_path --output DIR

The L lines of these GFAs are one contiguous block (d9: the end of the file; d46: between the S and the W lines), so only
the byte range --l-range is read (BeeGFS reads the whole 375 GB d46 GFA at 15-45 MB/s). The checks below prove the
block complete for the reference paths.

Outputs in DIR (new files only; the GFA, GBZ and graph index are read only):
    edges/succ_offsets.npy  int64 [2 * (max_node + 1) + 1]; side s = 2 * node + reverse
    edges/succ_targets.npy  int32 sides: the successors of side s are succ_targets[off[s]:off[s+1]]
                            (an L line a o1 b o2 gives a(o1) -> b(o2) and b(!o2) -> a(!o1))
    edges/meta.json         counts, source stamp, checks: every step of the GRCh38 and CHM13 paths is an edge
    chm13_path/             lengths, chrom, start0, reverse, visits, path_nodes, path_starts, path_reverse, meta.json
                            (format gfa-reference-path, reference_sample CHM13; one contig entry per CHM13 path fragment)
    chm13.gaf               vg paths -A -S CHM13 output (deleted at the end unless --keep)
Node lengths come from the GRCh38 reference-path directory of the same GFA (its lengths.npy has every segment).
"""
import argparse
import json
import os
import re
import subprocess
import sys
import time
from pathlib import Path

import numpy as np

REPO = Path(__file__).resolve().parents[4]
sys.path.insert(0, str(REPO))
from indexed_gam_pipeline_v4.tools.graph_prep import parse_walk  # noqa: E402
from indexed_gam_pipeline_v4.tensor_postprocessing.reference_path import FORMAT, ReferencePath  # noqa: E402

NAME = re.compile(r"^([^#]+)#([^#]+)#([^\[]+)(?:\[(\d+)(?:-(\d+))?\])?$")


def stamp(path):
    s = os.stat(path)
    return dict(path=str(Path(path).resolve()), size=s.st_size, mtime_ns=s.st_mtime_ns)


def build_edges(gfa, start, end, out, max_node):
    out.mkdir(parents=True, exist_ok=True)
    t = time.time()
    txt = out / "edges.txt"
    awk = 'BEGIN{FS="\\t"} $1=="L"{print $2, ($3=="-")?1:0, $4, ($5=="-")?1:0}'
    cmd = (f"tail -c +{start + 1} {gfa} | head -c {end - start} | LC_ALL=C grep '^L\t' | "
           f"LC_ALL=C awk '{awk}' > {txt}")
    subprocess.run(["bash", "-o", "pipefail", "-c", cmd], check=False)  # head closes the pipe: tail gets SIGPIPE
    e = np.fromfile(txt, sep=" ", dtype=np.int64).reshape(-1, 4)
    txt.unlink()
    print(f"{len(e):,} L lines in bytes {start:,}-{end:,} parsed in {time.time() - t:.0f} s", flush=True)
    if len(e) == 0 or e[:, [0, 2]].max() > max_node or e[:, [0, 2]].min() < 1:
        raise SystemExit(f"edges: {len(e)} lines, node range wrong")
    a, o1, b, o2 = e[:, 0], e[:, 1], e[:, 2], e[:, 3]
    n_l = len(e)
    del e
    src = np.concatenate([2 * a + o1, 2 * b + (1 - o2)])
    dst = np.concatenate([2 * b + o2, 2 * a + (1 - o1)])
    del a, o1, b, o2
    n_sides = 2 * (max_node + 1)
    key = np.unique(src * n_sides + dst)
    del src, dst
    src, dst = key // n_sides, key % n_sides
    del key
    offsets = np.zeros(n_sides + 1, dtype=np.int64)
    np.cumsum(np.bincount(src, minlength=n_sides), out=offsets[1:])
    if n_sides >= 2**31 - 1:
        raise SystemExit("too many sides for int32 targets")
    np.save(out / "succ_offsets.npy", offsets)
    np.save(out / "succ_targets.npy", dst.astype(np.int32))
    meta = dict(format="E9a-side-csr", source=stamp(gfa), l_range=[start, end], max_node=int(max_node),
                sides=int(n_sides), l_lines=n_l, directed_side_edges=int(len(dst)), seconds=time.time() - t, checks={})
    (out / "meta.json").write_text(json.dumps(meta, indent=2) + "\n")
    print(f"edges: {meta}", flush=True)


def has_edges(offsets, targets, src_sides, dst_sides):
    """For each (src, dst) pair: is dst among the successors of src (vectorised by sorting)."""
    lo, hi = offsets[src_sides], offsets[src_sides + 1]
    ok = np.zeros(len(src_sides), dtype=bool)
    # most sides have 1-3 successors: check the first 8 slots directly, then the rest one by one
    for k in range(8):
        idx = lo + k
        valid = idx < hi
        hit = np.zeros(len(src_sides), dtype=bool)
        hit[valid] = targets[idx[valid]] == dst_sides[valid]
        ok |= hit
    for i in np.flatnonzero(~ok & (hi - lo > 8)):
        ok[i] = dst_sides[i] in targets[lo[i]:hi[i]]
    return ok


def check_path_edges(out, path_dir, label):
    """Every consecutive step of every walk of a reference-path directory is an edge of the CSR."""
    off = np.load(out / "succ_offsets.npy")
    tgt = np.load(out / "succ_targets.npy")
    rp = ReferencePath(path_dir)
    missing = total = 0
    for c in rp.meta["contigs"]:
        sl = slice(c["offset"], c["offset"] + c["nodes"])
        nodes = np.asarray(rp.path_nodes[sl], dtype=np.int64)
        rev = np.asarray(rp.path_reverse[sl], dtype=np.int64)
        if len(nodes) < 2:
            continue
        s = 2 * nodes[:-1] + rev[:-1]
        d = 2 * nodes[1:] + rev[1:]
        ok = has_edges(off, tgt, s, d)
        total += len(ok)
        missing += int((~ok).sum())
    meta = json.loads((out / "meta.json").read_text())
    meta["checks"][label] = dict(steps=total, missing_edges=missing, passed=missing == 0)
    if "sides_without_successor" not in meta["checks"]:
        lengths = np.asarray(rp.lengths)
        nodes = np.flatnonzero(lengths > 0)
        deg = np.diff(off)
        meta["checks"]["nodes"] = int(len(nodes))
        meta["checks"]["sides_without_successor"] = int((deg[2 * nodes] == 0).sum() + (deg[2 * nodes + 1] == 0).sum())
    (out / "meta.json").write_text(json.dumps(meta, indent=2) + "\n")
    print(f"check {label}: {total:,} steps, {missing:,} without an edge", flush=True)
    return missing


def build_chm13(gaf, out, lengths, gbz, sample):
    t = time.time()
    out.mkdir(parents=True, exist_ok=True)
    size = len(lengths)
    reference = []
    with open(gaf, "rb") as f:
        for line in f:
            if line.startswith(b"@"):  # GAF header lines (@HD, @RG ...)
                continue
            fields = line.rstrip(b"\n").split(b"\t")
            if fields[7] != b"0" or fields[8] != fields[6] or fields[1] != fields[6]:
                raise SystemExit(f"GAF record {fields[0][:60]!r} is not a whole path (lengths {fields[1]}, {fields[6:9]})")
            m = NAME.match(fields[0].decode())
            if not m or m.group(1) != sample:
                raise SystemExit(f"unexpected GAF name {fields[0][:80]!r}")
            ids, rev = parse_walk(fields[5])
            start = int(m.group(4) or 0)
            qlen = int(fields[1])
            rec = dict(sample=sample, hap=m.group(2), contig=m.group(3), start=start, end=start + qlen,
                       nodes=int(ids.size), min=int(ids.min()), max=int(ids.max()))
            reference.append((rec, ids, rev))
    if not reference:
        raise SystemExit(f"no {sample} paths in {gaf}")
    chrom = np.full(size, -1, dtype=np.int32)
    start0 = np.zeros(size, dtype=np.int64)
    reverse_of = np.zeros(size, dtype=bool)
    visits = np.zeros(size, dtype=np.int64)
    contigs, offset, order = [], 0, []
    for k, (rec, ids, rev) in enumerate(reference):
        if ids.max() >= size or not lengths[ids].all():
            raise SystemExit(f"path {rec['contig']}:{rec['start']} uses nodes without segments")
        span = lengths[ids].astype(np.int64)
        starts = rec["start"] + np.cumsum(span) - span
        if int(span.sum()) != rec["end"] - rec["start"]:
            raise SystemExit(f"path {rec['contig']}:{rec['start']}: node lengths {span.sum()} != GAF length")
        unique, first, counts = np.unique(ids, return_index=True, return_counts=True)
        visits[unique] += counts
        new_mask = chrom[unique] == -1
        chrom[unique[new_mask]] = k
        start0[unique[new_mask]] = starts[first[new_mask]]
        reverse_of[unique[new_mask]] = rev[first[new_mask]]
        contigs.append(dict(name=rec["contig"], hap=rec["hap"], start=rec["start"], end=rec["end"], nodes=int(ids.size),
                            min_node=rec["min"], max_node=rec["max"], offset=offset))
        offset += ids.size
        order.append((ids, starts, rev))
    arrays = dict(lengths=lengths.astype(np.int32), chrom=chrom, start0=start0, reverse=reverse_of,
                  visits=np.minimum(visits, np.iinfo(np.uint16).max).astype(np.uint16),
                  path_nodes=np.concatenate([o[0] for o in order]), path_starts=np.concatenate([o[1] for o in order]),
                  path_reverse=np.concatenate([o[2] for o in order]))
    for name, values in arrays.items():
        np.save(out / f"{name}.npy", values)
    meta = dict(format=FORMAT, source=stamp(gbz), reference_sample=sample, max_node=size - 1,
                segments=int(np.count_nonzero(lengths)), walks=len(reference), samples=1, contigs=contigs,
                reference_nodes=int(np.count_nonzero(visits)), ambiguous_reference_nodes=int(np.count_nonzero(visits > 1)),
                scan_seconds=time.time() - t, checks={},
                note="E9a prep_graph.py: vg paths -A -S CHM13 on the GBZ; one contig entry per path fragment")
    (out / "meta.json").write_text(json.dumps(meta, indent=2) + "\n")
    print(f"{sample}: {len(reference)} path fragments, {meta['reference_nodes']:,} nodes, "
          f"{meta['ambiguous_reference_nodes']:,} visited more than once, {time.time() - t:.0f} s", flush=True)


def main():
    p = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    p.add_argument("--gfa", required=True)
    p.add_argument("--l-range", type=int, nargs=2, required=True, metavar=("START", "END"))
    p.add_argument("--gbz", required=True)
    p.add_argument("--vg", required=True)
    p.add_argument("--grch38-path", required=True, help="the graph's GRCh38 reference-path directory (lengths.npy)")
    p.add_argument("--output", required=True)
    p.add_argument("--sample", default="CHM13")
    p.add_argument("--threads", type=int, default=2)
    p.add_argument("--keep", action="store_true", help="keep chm13.gaf")
    args = p.parse_args()
    out = Path(args.output)
    out.mkdir(parents=True, exist_ok=True)
    ref_meta = json.loads((Path(args.grch38_path) / "meta.json").read_text())
    if Path(ref_meta["source"]["path"]).resolve() != Path(args.gfa).resolve():
        raise SystemExit(f"{args.grch38_path} was scanned from {ref_meta['source']['path']}, not {args.gfa}")
    lengths = np.load(Path(args.grch38_path) / "lengths.npy")
    if not (out / "edges" / "meta.json").exists():
        build_edges(args.gfa, args.l_range[0], args.l_range[1], out / "edges", len(lengths) - 1)
    check_path_edges(out / "edges", args.grch38_path, "grch38_path_steps")
    gaf = out / "chm13.gaf"
    if not (out / "chm13_path" / "meta.json").exists():
        if not gaf.exists():
            t = time.time()
            with open(str(gaf) + ".tmp", "wb") as f:
                subprocess.run([args.vg, "paths", "-x", args.gbz, "-S", args.sample, "-A", "-t", str(args.threads)],
                               stdout=f, check=True)
            os.rename(str(gaf) + ".tmp", gaf)
            print(f"vg paths -A -S {args.sample}: {gaf.stat().st_size / 1e6:.0f} MB in {time.time() - t:.0f} s", flush=True)
        build_chm13(gaf, out / "chm13_path", lengths, args.gbz, args.sample)
    check_path_edges(out / "edges", out / "chm13_path", "chm13_path_steps")
    if not args.keep and gaf.exists():
        gaf.unlink()
    print("done", flush=True)


if __name__ == "__main__":
    main()
