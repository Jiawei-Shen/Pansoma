"""E9j step 3: class of every COLO829T patient-frame somatic SNV event on one graph, and the candidate keys of the
off-reference ones. The logic of assembly_offref_20260929/s3_classify.py + s4 (HG008T, d9) and
E1.1_ch6_counterfactual/asm_d46_classify.py (HG008T, d46), unchanged: each event gives two reads, the event haplotype's
window (N) and the same with the somatic allele (T), 100-bp flanks (s1_windows.py), aligned by giraffe (s2_giraffe.sh);
  A  a somatic edit (a difference in T's alignment that N's does not have) on a GRCh38 node (one GRCh38 visit)
  B  a somatic edit on a node off the GRCh38 path: an off-reference somatic variant
  C  no somatic edit, but T and N walk different nodes there (the allele is a graph branch: no candidate possible)
  D  unmapped, MAPQ < 20, or the same alignment
A B event's key is its edit in node-forward coordinates, node:start:KIND:REF>ALT (the candidate id form).
usage: python s3_classify.py {d9|d46} DATA_DIR   (writes classes_<g>.tsv, offref_somatic_<g>.tsv)
"""
import csv
import json
import sys
from collections import Counter
from pathlib import Path

CODE = Path("/scratch/jshen/data/reviewer_experiments_20261013/E9j_offref_colo829t/code_8580a9f")
sys.path.insert(0, str(CODE))
from indexed_gam_pipeline_v4.graph_index import GraphIndex  # noqa: E402
from indexed_gam_pipeline_v4.tensor_postprocessing.reference_path import ReferencePath  # noqa: E402

GRAPHS = {"d9": ("/scratch/jshen/data/pansoma_v2_tensors/graph_index/hprc-v1.1-mc-grch38.d9.grch38_path",
                 "/scratch/jshen/data/pansoma_v2_tensors/graph_index/hprc-v1.1-mc-grch38.d9.graph_index.sqlite"),
          "d46": ("/scratch/jshen/data/pansoma_v2_tensors/graph_index_hprc-v2.1-d46/hprc-v2.1-mc-grch38.d46.grch38_path",
                  "/scratch/jshen/data/pansoma_v2_tensors/graph_index_hprc-v2.1-d46/hprc-v2.1-mc-grch38.d46.graph_index.sqlite")}
COMP = str.maketrans("ACGTN", "TGCAN")


def rc(s):
    return s.translate(COMP)[::-1]


def walk(a):
    """[(read_start, read_end, node, rev, node_pos_traversal, from_len, to_len, seq, is_diff)] over the alignment."""
    out, cur = [], 0
    for m in a.get("path", {}).get("mapping", []):
        p = m["position"]
        node, rev, off = int(p["node_id"]), bool(p.get("is_reverse", False)), int(p.get("offset", 0))
        for e in m.get("edit", []):
            fl, tl, seq = int(e.get("from_length", 0)), int(e.get("to_length", 0)), e.get("sequence", "")
            out.append((cur, cur + tl, node, rev, off, fl, tl, seq, fl != tl or bool(seq)))
            cur += tl
            off += fl
    return out


def nodes_at(w, x):
    return {s[2] for s in w if s[0] <= x < s[1] or (s[6] == 0 and s[0] == x)}


def main():
    g, out = sys.argv[1], Path(sys.argv[2])
    rp = ReferencePath(GRAPHS[g][0])

    def on_ref(n):
        return 0 < n < len(rp.visits) and int(rp.visits[n]) == 1

    flank = {}
    with open(out / "windows.fq") as f:
        for k, line in enumerate(f):
            if k % 4 == 0:
                name, _, fl = line[1:].split()
                flank[name] = int(fl.split("=")[1])
    events = {r["event"]: r for r in csv.DictReader(open(out / "events.tsv"), delimiter="\t")}
    aln = {}
    for line in open(out / f"windows.{g}.gam.json"):
        a = json.loads(line)
        aln[a["name"]] = a
    rows, cls = [], Counter()
    for eid, ev in events.items():
        t, n = aln.get(f"{eid}_T"), aln.get(f"{eid}_N")
        F, L = flank[f"{eid}_T"], max(1, len(ev["alt"]))
        row = dict(event=eid, chrom=ev["chrom"], grch38_pos=ev["grch38_pos"], ref=ev["ref"], alt=ev["alt"], hap=ev["hap"],
                   how=ev["how"], contig=ev["contig"], pos0=ev["pos0"], grch38_records=ev["grch38_records"], cls="",
                   edit_node="", edit_rev="", edit_offset="", edit="", background_on_ref="",
                   mapq_T=(t or {}).get("mapping_quality", 0), mapq_N=(n or {}).get("mapping_quality", 0))
        if not t or not n or not t.get("path") or not n.get("path"):
            row["cls"] = "D_unmapped"
        elif int(row["mapq_T"]) < 20 or int(row["mapq_N"]) < 20:
            row["cls"] = "D_mapq<20"
        else:
            wt, wn = walk(t), walk(n)
            bg = nodes_at(wn, F)
            row["background_on_ref"] = all(on_ref(x) for x in bg) if bg else ""
            ndiff = {(s[2], s[4], s[5], s[7]) for s in wn if s[8]}
            tonly = [s for s in wt if s[8] and (s[2], s[4], s[5], s[7]) not in ndiff]
            som = sorted(tonly, key=lambda s: abs((s[0] + s[1]) / 2 - (F + L / 2)))
            if som:
                s = som[0]
                row.update(cls="A" if on_ref(s[2]) else "B", edit_node=s[2], edit_rev=s[3], edit_offset=s[4],
                           edit=f"{s[5]}>{s[6]}:{s[7]}")
            elif nodes_at(wt, F) != bg or [s[2] for s in wt] != [s[2] for s in wn]:
                row["cls"] = "C"
            else:
                row["cls"] = "D_same_alignment"
        rows.append(row)
        cls[row["cls"]] += 1
    with open(out / f"classes_{g}.tsv", "w") as f:
        w = csv.DictWriter(f, fieldnames=list(rows[0]), delimiter="\t")
        w.writeheader()
        w.writerows(rows)
    print(g, "events", len(rows), dict(sorted(cls.items())))
    print(g, "by event-haplotype choice:", dict(sorted(Counter((r["how"], r["cls"]) for r in rows).items())))
    b = [r for r in rows if r["cls"] == "B"]
    with GraphIndex(GRAPHS[g][1]) as gi:
        seqs = {k: v["sequence"] for k, v in gi.get_nodes({int(r["edit_node"]) for r in b}).items()}
    keys = []
    for r in b:
        node, rev, off = int(r["edit_node"]), r["edit_rev"] == "True", int(r["edit_offset"])
        fl_tl, seq = r["edit"].split(":", 1)
        fl, tl = map(int, fl_tl.split(">"))
        s, L = seqs[node], len(seqs[node])
        if fl == tl == 1:
            start = L - off - 1 if rev else off
            kind, ref, alt = "SNP", s[start], (seq.translate(COMP) if rev else seq)
        elif fl == 0:
            start = L - off if rev else off
            kind, ref, alt = "INS", "", (rc(seq) if rev else seq)
        elif tl == 0:
            start = L - off - fl if rev else off
            kind, ref, alt = "DEL", s[start:start + fl], ""
        else:
            start, kind, ref, alt = (L - off - fl if rev else off), "COMPLEX", "", ""
        keys.append(dict(event=r["event"], contig=r["contig"], pos=r["pos0"], ref=r["ref"], alt=r["alt"], hap=r["hap"],
                         how=r["how"], chrom=r["chrom"], grch38_records=r["grch38_records"],
                         normal_offref_at_site=r["background_on_ref"] is False, node=node, start=start, kind=kind,
                         key=f"{node}:{start}:{kind}:{ref}>{alt}"))
    with open(out / f"offref_somatic_{g}.tsv", "w") as f:
        w = csv.DictWriter(f, fieldnames=list(keys[0]) if keys else ["event"], delimiter="\t")
        w.writeheader()
        w.writerows(keys)
    print(g, "off-reference somatic (B):", len(keys), dict(Counter(k["kind"] for k in keys)),
          "| chr1:", sum(k["chrom"] == "chr1" for k in keys),
          "| normal haplotype off-reference at the site:", sum(r["background_on_ref"] is False for r in b))


if __name__ == "__main__":
    main()
