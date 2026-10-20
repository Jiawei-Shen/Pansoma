"""E9b: a small GFA with what graph_prep ref-path-scan reads for the CHM13 reference path of one graph, made without
reading the graph's whole GFA (d46: 375 GB; BeeGFS gave ~10 MB/s on 2026-10-19).

    python chm13_mini_gfa.py --lengths <graph>.grch38_path/lengths.npy --gaf chm13.gaf --output <graph>.chm13_only.gfa
    python -m indexed_gam_pipeline_v4.tools.graph_prep ref-path-scan --gfa <graph>.chm13_only.gfa --output DIR \
        --reference-sample CHM13

- S lines: one per segment, `S <id> * LN:i:<length>`, from lengths.npy of the graph's GRCh38 reference-path directory
  (ref-path-scan of the same graph's GFA, so the segment IDs and lengths are those of the GFA).
- W lines: one per record of `vg paths -A -S CHM13` on the GBZ (E9a prep_graph.py left it as <scratch>/graph/<g>/chm13.gaf):
  query name CHM13#<hap>#<contig>[<start>] (or [<start>-<end>], or no range = start 0), path = the walk (>id<id...),
  which must be aligned whole and forward. W = sample, hap, contig, start, start + path length, walk.
ref-path-scan then checks that the node lengths add up to every W interval, and ref-path-check compares 100,000 random
CHM13 nodes with the GBZ graph index (lengths) and the CHM13 FASTA (sequence).
"""
import argparse
import re

import numpy as np

NAME = re.compile(r"^([^#]+)#([^#]+)#(.+?)(?:\[(\d+)(?:-(\d+))?\])?$")


def main():
    p = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    p.add_argument("--lengths", required=True)
    p.add_argument("--gaf", required=True)
    p.add_argument("--output", required=True)
    a = p.parse_args()
    lengths = np.load(a.lengths)
    ids = np.flatnonzero(lengths)
    walks = 0
    with open(a.output, "x") as out:
        out.write("H\tVN:Z:1.1\n")
        for chunk in np.array_split(ids, max(1, ids.size // 2_000_000)):
            out.write("".join(f"S\t{i}\t*\tLN:i:{n}\n" for i, n in zip(chunk.tolist(), lengths[chunk].tolist())))
        with open(a.gaf) as gaf:
            for line in gaf:
                if line.startswith("@"):
                    continue
                q = line.rstrip("\n").split("\t")
                m = NAME.match(q[0])
                if not m:
                    raise SystemExit(f"unexpected path name {q[0]}")
                sample, hap, contig, start, end = m.groups()
                start, length = int(start or 0), int(q[1])
                if (q[2], q[3], q[4]) != ("0", q[1], "+") or (end is not None and int(end) != start + length):
                    raise SystemExit(f"{q[0]}: not one whole forward path ({q[1:5]})")
                out.write(f"W\t{sample}\t{hap}\t{contig}\t{start}\t{start + length}\t{q[5]}\n")
                walks += 1
    print(f"{a.output}: {ids.size:,} S lines, {walks} W lines")


if __name__ == "__main__":
    main()
