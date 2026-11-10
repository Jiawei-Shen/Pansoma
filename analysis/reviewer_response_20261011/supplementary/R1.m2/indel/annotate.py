"""Caption facts for the chosen R1.m2 INDEL examples (picks.tsv from choose.py): per example, the truth record, the
allele's nodes, the GRCh38-path nodes under the window, and at each branch column the channel-6 code of the ALT (A1)
and REF rows decoded to haplotype counts. A code is matched to nodes near the window (GRCh38-path nodes of the window
and node ids within +-3 of them) whose distinct path count encodes to that code, so an exact count is given when one
node matches. Read-only. usage: python annotate.py picks.tsv > annotations.md"""
import csv
import json
import sqlite3
import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[5] / "scripts"))
import visualize_tensor as v  # noqa: E402
from select_examples import D, chr1_path, dominant  # noqa: E402


def decode(code, h):
    return [c for c in range(1, h + 1) if v.haplotype_count_code(c, h) == code]


def main():
    picks = list(csv.DictReader(open(sys.argv[1]), delimiter="\t"))
    for p in picks:
        kd = D / p["set"] / "tensors" / "INDEL"
        man = json.load(open(kd / "manifest.json"))
        h = man["parameters"]["haplotypes"]
        con = sqlite3.connect(f"file:{man['graph_index']['path']}?mode=ro", uri=True)
        sm = None
        with open(kd / "chr1_variant_summary.ndjson") as f:
            key = f'{{"candidate_id": "{p["candidate_id"]}"'
            for line in f:
                if line.startswith(key + ","):
                    sm = json.loads(line)
                    break
        x = np.asarray(np.load(kd / p["shard_file"], mmap_mode="r")[int(p["index"])])
        g = {b["allele"]: (b["start_row"], b["end_row"]) for b in sm["row_groups"]}
        a0, a1 = g["A1"]; r0, r1 = g["REF"]
        starts, nodes = chr1_path(h)
        lo = sm["grch38"]["pos0"] - 50
        hi = sm["grch38"]["pos0"] + max(len(sm["grch38"]["ref"]), 1) + 50
        i0 = max(int(np.searchsorted(starts, lo, side="right")) - 1, 0)
        i1 = int(np.searchsorted(starts, hi, side="left"))
        ref_nodes = [int(n) for n in nodes[i0:i1]]
        near = sorted({n + d for n in ref_nodes for d in range(-3, 4)} | {sm["node_id"]} |
                      {s[0] for s in sm["path"]})
        info = {r[0]: (r[1], r[2]) for r in con.execute(
            f"select node_id, seq, distinct_path_count from nodes where node_id in ({','.join('?' * len(near))})", near)}
        allele_nodes = [sm["node_id"]] + [s[0] for s in sm["path"]]
        print(f"## {p['fig']}  {p['set']}  {p['candidate_id']}  (shard {p['shard_file']} index {p['index']})\n")
        print(f"- truth {p['truth']}; {sm['event_type']} {sm['ref'] or '-'}>{sm['alt'] or '-'}; AF {sm['af']:.3f}; "
              f"rows A1 {a1 - a0}, REF {r1 - r0}, other {sm['site_counts'].get('OTHER', 0)}; H = {h}")
        print(f"- candidate node {sm['node_id']} ({len(info[sm['node_id']][0])} bp, {info[sm['node_id']][1]} paths), "
              f"site offset {sm['start']}-{sm['end']}; allele nodes {allele_nodes}; candidate columns "
              f"{sm['candidate_columns']}")
        print(f"- GRCh38-path nodes under the window: {len(ref_nodes)}: " + ", ".join(
            f"{n} ({len(info[n][0])} bp, {info[n][1]})" for n in ref_nodes))
        ch6 = x[6].astype(int)
        site = set(range(sm["candidate_columns"][0], sm["candidate_columns"][1] + 1))
        print("- branch columns (any row on a node with fewer than H paths), runs of equal (A1, REF) codes:")
        bases = "-ACGTN-"

        def name(code):
            if code is None:
                return "-"
            m = sorted(n for n, (s, k) in info.items() if v.haplotype_count_code(k, h) == code)
            d = decode(code, h)
            cnt = d[0] if len(d) == 1 else f"{d[0]}-{d[-1]}"
            return f"code {code} = {cnt} paths" + (f" (nodes {m})" if m else "")

        runs = []
        for c in range(ch6.shape[1]):
            col = ch6[:, c]
            if not ((col > 0) & (col < 100)).any():
                continue
            va, fa = dominant(col[a0:a1][col[a0:a1] > 0])
            vr, fr = dominant(col[r0:r1][col[r0:r1] > 0])
            ba = dominant(x[5, a0:a1, c][x[5, a0:a1, c] > 0])[0]
            br = dominant(x[5, r0:r1, c][x[5, r0:r1, c] > 0])[0]
            key = (va, vr, c in site)
            if runs and runs[-1][0] == key and runs[-1][2] == c - 1:
                runs[-1][2] = c; runs[-1][3].append((fa, fr)); runs[-1][4] += bases[ba or 0]; runs[-1][5] += bases[br or 0]
            else:
                runs.append([key, c, c, [(fa, fr)], bases[ba or 0], bases[br or 0]])
        for (va, vr, st), c0, c1, fr, sa, sr in runs:
            mark = " SITE" if st else (" DIFF" if va != vr else "")
            print(f"  - cols {c0}-{c1}{mark}: A1 {name(va)} [{min(f[0] for f in fr):.0%}+, {sa}] | "
                  f"REF {name(vr)} [{min(f[1] for f in fr):.0%}+, {sr}]")
        print()


if __name__ == "__main__":
    main()
