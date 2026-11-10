"""R1.m2 INDEL examples whose allele itself spans several graph nodes (a deletion through >= 2 nodes).

From candidates.tsv (select_examples.py): tensors with allele_nodes >= 2 (the summary's allele path), >= 5 ALT and
>= 5 REF rows, and fewer than 10% of the REF rows carrying the ALT reads' sequence (choose.hidden_alt). One tensor per
truth INDEL: the most allele nodes, then the larger of min(ALT, REF) rows. Writes span_picks.tsv (the picks.tsv
columns, figures R1.m2_INDEL_span_NN) and prints each allele's nodes (length, haplotype paths) to span_nodes.md.
usage: python span_examples.py candidates.tsv span_picks.tsv span_nodes.md"""
import csv
import json
import sqlite3
import sys

from choose import hidden_alt, label, row_groups
from select_examples import D


def main():
    rows = [r for r in csv.DictReader(open(sys.argv[1]), delimiter="\t")
            if int(r["allele_nodes"]) >= 2 and int(r["alt_rows"]) >= 5 and int(r["ref_rows"]) >= 5]
    groups = row_groups(rows)
    best = {}
    for r in rows:
        h = hidden_alt(r, groups[r["set"], r["candidate_id"]])
        if h >= 0.1:
            print(f"dropped {r['set']} {r['candidate_id']}: {h:.0%} of REF rows carry the ALT sequence", file=sys.stderr)
            continue
        key = (int(r["allele_nodes"]), min(int(r["alt_rows"]), int(r["ref_rows"])))
        if r["truth"] not in best or key > best[r["truth"]][0]:
            best[r["truth"]] = (key, r)
    picks = sorted((r for _, r in best.values()), key=lambda r: (-int(r["allele_nodes"]), -int(r["len"]), r["truth"]))
    cols = ["fig", "set", "label", "candidate_id", "shard_file", "index", "truth", "event", "len", "af", "alt_rows",
            "ref_rows", "node_len", "start", "allele_nodes", "ref_nodes", "strong_cols"]
    with open(sys.argv[2], "w") as out, open(sys.argv[3], "w") as md:
        print("\t".join(cols), file=out)
        for i, r in enumerate(picks, 1):
            fig = f"R1.m2_INDEL_span_{i:02d}"
            print("\t".join(map(str, [fig, r["set"], label(r["set"]), r["candidate_id"], r["shard_file"], r["index"],
                                      r["truth"], r["event"], r["len"], r["af"], r["alt_rows"], r["ref_rows"],
                                      r["node_len"], r["start"], r["allele_nodes"], r["ref_nodes"], ""])), file=out)
            kd = D / r["set"] / "tensors" / "INDEL"
            man = json.load(open(kd / "manifest.json"))
            con = sqlite3.connect(f"file:{man['graph_index']['path']}?mode=ro", uri=True)
            key = f'{{"candidate_id": "{r["candidate_id"]}",'
            with open(kd / "chr1_variant_summary.ndjson") as f:
                sm = next(json.loads(line) for line in f if line.startswith(key))
            parts = [(sm["node_id"], sm["start"], None)] + [tuple(p) for p in sm["path"]]
            desc = []
            for node, start, end in parts:
                seq, paths = con.execute("select seq, distinct_path_count from nodes where node_id=?", (node,)).fetchone()
                desc.append(f"{node} ({len(seq)} bp {seq if len(seq) <= 12 else seq[:12] + '...'}, {paths} paths)")
            print(f"- {fig} {r['set']} {r['candidate_id']} {r['truth']}: ALT {r['alt_rows']}, REF {r['ref_rows']}, "
                  f"AF {float(r['af']):.2f}, H {man['parameters']['haplotypes']}; allele nodes: " + ", ".join(desc)
                  + f"; summary path {sm['path']}; GRCh38-path nodes under the window {r['ref_nodes']}", file=md)
    print(f"{len(rows)} tensors with the allele on >= 2 nodes, {len(picks)} truth INDELs picked", file=sys.stderr)


if __name__ == "__main__":
    main()
