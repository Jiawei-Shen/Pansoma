"""Pick the ten R1.m2 INDEL examples from candidates.tsv (select_examples.py).

Kept: >= 10 ALT (A1) and >= 10 REF rows, the window over >= 3 GRCh38-path nodes, and at least one strong allele-linked
branch column: outside the site, >= 80% of ALT rows and >= 80% of REF rows each share one channel-6 code, the two codes
differ by >= 10, and neither code is the other allele's dominant code within 3 columns (the ALT and REF reads pass
through different nodes there, not the same nodes shifted by the INDEL). Ranked by: allele itself on >= 2 nodes,
number of strong columns, GRCh38-path nodes under the window, ALT rows. One example per truth INDEL; first the best of
each tensor set, then a second from the sets with more examples, until ten.
usage: python choose.py candidates.tsv picks.tsv"""
import csv
import json
import sys

import numpy as np

from select_examples import D, dominant

PLATFORM = {"PacBio": "PacBio HiFi", "fiberseq": "Fiber-seq", "ONT": "ONT", "Illumina": "Illumina"}


def strong(detail):
    out = []
    for item in filter(None, detail.split(",")):
        col, codes = item.split(":")
        a, r = map(int, codes.split("/"))
        if abs(a - r) >= 10:
            out.append((int(col), a, r))
    return out


def row_groups(rows):
    """{(set, candidate_id): {allele: (start_row, end_row)}} from the chr1 summaries."""
    out = {}
    for s in sorted({r["set"] for r in rows}):
        want = {f'{{"candidate_id": "{r["candidate_id"]}",' for r in rows if r["set"] == s}
        with open(D / s / "tensors" / "INDEL" / "chr1_variant_summary.ndjson") as f:
            for line in f:
                if line[:line.index(",") + 1] in want:
                    sm = json.loads(line)
                    out[s, sm["candidate_id"]] = {b["allele"]: (b["start_row"], b["end_row"]) for b in sm["row_groups"]}
    return out


def unshifted(r, groups, cols):
    x = np.load(D / r["set"] / "tensors" / "INDEL" / r["shard_file"], mmap_mode="r")[int(r["index"])]
    ch6 = np.asarray(x[6]).astype(int)
    w = ch6.shape[1]
    dom = {k: [dominant(ch6[a:b, c][ch6[a:b, c] > 0])[0] for c in range(w)]
           for k, (a, b) in groups.items() if k in ("A1", "REF")}
    return [(c, a, b) for c, a, b in cols
            if a not in dom["REF"][max(c - 3, 0):c + 4] and b not in dom["A1"][max(c - 3, 0):c + 4]]


def label(s):
    tumor, platform = s.split("_")[:2]
    return f"{tumor} {PLATFORM[platform]}, HPRC {'v2.1 (d46)' if s.endswith('_d46') else 'v1.1 (d9)'}"


def main():
    rows = [r for r in csv.DictReader(open(sys.argv[1]), delimiter="\t")
            if int(r["alt_rows"]) >= 10 and int(r["ref_rows"]) >= 10 and int(r["ref_nodes"]) >= 3
            and strong(r["linked_detail"])]
    groups = row_groups(rows)
    keep = []
    for r in rows:
        st = unshifted(r, groups[r["set"], r["candidate_id"]], strong(r["linked_detail"]))
        if st:
            keep.append((int(r["allele_nodes"]) > 1, len(st), int(r["ref_nodes"]), int(r["alt_rows"]), r, st))
    keep.sort(key=lambda k: k[:4], reverse=True)
    picks, sets, truths = [], set(), set()
    for per_set in (1, 2):  # first one per set, then a second from the sets with more examples
        for *_, r, st in keep:
            if len(picks) == 10:
                break
            if sum(p[0]["set"] == r["set"] for p in picks) >= per_set or r["truth"] in truths:
                continue
            picks.append((r, st)); sets.add(r["set"]); truths.add(r["truth"])
    picks.sort(key=lambda p: (not p[0]["set"].startswith("HG008T"), p[0]["set"].endswith("_d46"), p[0]["set"],
                              p[0]["truth"]))
    with open(sys.argv[2], "w") as out:
        cols = ["fig", "set", "label", "candidate_id", "shard_file", "index", "truth", "event", "len", "af", "alt_rows",
                "ref_rows", "node_len", "start", "allele_nodes", "ref_nodes", "strong_cols"]
        print("\t".join(cols), file=out)
        for i, (r, st) in enumerate(picks, 1):
            print("\t".join(map(str, [f"R1.m2_INDEL_{i:02d}", r["set"], label(r["set"]), r["candidate_id"],
                                      r["shard_file"], r["index"], r["truth"], r["event"], r["len"], r["af"],
                                      r["alt_rows"], r["ref_rows"], r["node_len"], r["start"], r["allele_nodes"],
                                      r["ref_nodes"], ",".join(f"{c}:{a}/{b}" for c, a, b in st)])), file=out)
    print(f"{len(rows)} with strong columns, {len(keep)} unshifted, {len(picks)} picked from {len(sets)} sets",
          file=sys.stderr)


if __name__ == "__main__":
    main()
