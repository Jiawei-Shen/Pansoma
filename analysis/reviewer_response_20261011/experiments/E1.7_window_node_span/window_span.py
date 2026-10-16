"""E1.7: how often does a candidate's 101-column tensor window reach beyond the candidate's own node (R1.1, R1.m2)?

For every chr1 tensor of a set (tensors/<KIND>/chr1_variant_summary.ndjson: node_id, start, end of the site on the
node) the window covers 50 columns on each side of the site. It stays inside the node only if the node has >= 50
bases on both sides (start >= 50 and node length - end >= 50). Node lengths come from the set's graph index sqlite
(manifest.json graph_index.path). Somatic = label 1 in chr1_labels.ndjson (by candidate_id). Read-only; standard
library + sqlite3. usage: python window_span.py SET [SET ...] > tsv"""
import json
import sqlite3
import sys
from pathlib import Path

D = Path("/scratch/jshen/data/pansoma_v2_tensors")


def main():
    print("set\tkind\ttensors\twindow_beyond_node\tshare\tsomatic\tsomatic_beyond\tsomatic_share\tmedian_node_len")
    for s in sys.argv[1:]:
        for kind in ("SNV", "INDEL"):
            kd = D / s / "tensors" / kind
            if not (kd / "chr1_variant_summary.ndjson").exists():
                continue
            gi = json.load(open(kd / "manifest.json"))["graph_index"]["path"]
            con = sqlite3.connect(f"file:{gi}?mode=ro", uri=True)
            labels = {}
            with open(kd / "chr1_labels.ndjson") as f:
                for line in f:
                    r = json.loads(line)
                    labels[r["candidate_id"]] = r.get("label")
            recs = []
            with open(kd / "chr1_variant_summary.ndjson") as f:
                for line in f:
                    r = json.loads(line)
                    recs.append((r["candidate_id"], r["node_id"], r["start"], r["end"]))
            lens = {}
            ids = sorted({n for _, n, _, _ in recs})
            for i in range(0, len(ids), 50000):
                chunk = ids[i:i + 50000]
                q = f"select node_id, length(seq) from nodes where node_id in ({','.join('?' * len(chunk))})"
                lens.update(con.execute(q, chunk).fetchall())
            n = beyond = som = som_b = 0
            node_lens = []
            for cid, node, st, en in recs:
                L = lens[node]
                out = st < 50 or L - en < 50
                n += 1; beyond += out; node_lens.append(L)
                if labels.get(cid) == 1:
                    som += 1; som_b += out
            node_lens.sort()
            print(f"{s}\t{kind}\t{n}\t{beyond}\t{beyond / n:.4f}\t{som}\t{som_b}\t{(som_b / som if som else float('nan')):.4f}\t{node_lens[len(node_lens) // 2]}")
            sys.stdout.flush()


if __name__ == "__main__":
    main()
