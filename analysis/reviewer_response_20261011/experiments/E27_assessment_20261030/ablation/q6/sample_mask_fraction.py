"""Share of covered cells outside the candidate's own node visit (node_mask.own_node_mask) on chr1 of merged v4 sets:
an evenly strided sample of all chr1 tensors (every shard) and an evenly strided sample of the label-1 (somatic)
tensors. Also the E1.7 rule (start < 50 or node length - end < 50) on the same tensors, as a cross-check.
Read-only (mmap). usage: python sample_mask_fraction.py N SET/KIND [SET/KIND ...] > tsv"""
import json
import sqlite3
import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent))
from node_mask import own_node_mask  # noqa: E402

D = Path("/scratch/jshen/data/pansoma_v2_tensors")


def main():
    n_sample = int(sys.argv[1])
    print("set\tkind\tgroup\ttensors\twith_masked_cells\te17_rule_beyond_node\tagree\tcovered_cells\tmasked_cells\tmasked_share")
    for sk in sys.argv[2:]:
        s, kind = sk.split("/")
        kd = D / s / "tensors" / kind
        m = json.load(open(kd / "manifest.json"))
        gi, total = m["graph_index"]["path"], m["chromosomes"]["chr1"]["tensors"]
        labels = {sh["file"]: np.load(kd / sh["file"].replace("_data.npy", "_labels.npy"))
                  for sh in m["chromosomes"]["chr1"]["shards"]}
        want = set(np.unique(np.linspace(0, total - 1, n_sample).astype(int)).tolist())
        strided, somatic = [], []
        with open(kd / "chr1_variant_summary.ndjson") as f:
            for k, line in enumerate(f):
                r = json.loads(line)
                rec = (r["shard_file"], r["index_within_shard"], r["node_id"], r["start"], r["end"],
                       r["site_layout"]["insertion_slots"])
                if k in want:
                    strided.append(rec)
                if labels[rec[0]][rec[1]] == 1:
                    somatic.append(rec)
        if len(somatic) > n_sample:
            somatic = [somatic[i] for i in np.unique(np.linspace(0, len(somatic) - 1, n_sample).astype(int))]
        con = sqlite3.connect(f"file:{gi}?mode=ro", uri=True)
        nodes = sorted({r[2] for r in strided + somatic})
        lens = dict(con.execute(f"select node_id, length(seq) from nodes where node_id in ({','.join('?' * len(nodes))})",
                                nodes).fetchall())
        maps = {}
        for group, recs in (("all_strided", strided), ("somatic_strided", somatic)):
            t = tm = geo_n = agree = cov = msk = 0
            for f, i, node, start, end, slots in recs:
                if f not in maps:
                    maps[f] = np.load(kd / f, mmap_mode="r")
                x = np.asarray(maps[f][i])
                L = lens[node]
                own = own_node_mask(x, start, L, slots)
                covered = x[4] > 0
                mm = int((covered & ~own).sum())
                geo = start < 50 or L - end < 50
                t += 1; tm += mm > 0; geo_n += geo; agree += geo == (mm > 0); cov += int(covered.sum()); msk += mm
            print(f"{s}\t{kind}\t{group}\t{t}\t{tm}\t{geo_n}\t{agree}\t{cov}\t{msk}\t{msk / max(cov, 1):.4f}", flush=True)


if __name__ == "__main__":
    main()
