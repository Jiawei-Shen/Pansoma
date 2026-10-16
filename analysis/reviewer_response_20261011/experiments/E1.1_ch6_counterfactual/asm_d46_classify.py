"""The HG008-N assembly somatic events placed on the HPRC v2.1 d46 graph, to find the d46 candidate keys of the
off-reference somatic SNVs (for the d46 models' assembly scores in E1.1).

The same procedure as pansoma_net_v2_runs/assembly_offref_20260929 (s3_classify.py, s4_tensors.py) did on HPRC v1.1
d9, with that directory's windows.fq and events.tsv (read only): every event gives two reads, the normal haplotype
window (N) and the same with the somatic allele (T), 100 bp flanks; giraffe aligns them to d46 (the job script
asm_d46.sh); an event is
  A  a somatic edit (a difference in T's alignment that N's does not have) on a GRCh38 node (one GRCh38 visit)
  B  a somatic edit on a node off the GRCh38 path: an off-reference somatic variant
  C  no somatic edit, but T and N walk different nodes there (the allele is a graph branch)
  D  unmapped, MAPQ < 20, or the same alignment
and a B event's key is its edit in node-forward coordinates, node:start:KIND:REF>ALT (the candidate id form).

usage: python asm_d46_classify.py OUT_DIR   (OUT_DIR holds windows.d46.gam.json; writes classes_d46.tsv and
       offref_somatic_d46.tsv there)"""
import csv
import json
import sys
from collections import Counter

sys.path.insert(0, "/scratch/jshen/data/reviewer_experiments_20261013/E1.1_ch6_counterfactual/code_41f7830")  # frozen 41f7830
from indexed_gam_pipeline_v4.graph_index import GraphIndex  # noqa: E402
from indexed_gam_pipeline_v4.tensor_postprocessing.reference_path import ReferencePath  # noqa: E402

A = "/scratch/jshen/data/pansoma_net_v2_runs/assembly_offref_20260929"
G = "/scratch/jshen/data/pansoma_v2_tensors/graph_index_hprc-v2.1-d46"
O = sys.argv[1]
rp = ReferencePath(f"{G}/hprc-v2.1-mc-grch38.d46.grch38_path")
COMP = str.maketrans("ACGTN", "TGCAN")


def rc(s):
    return s.translate(COMP)[::-1]


flank = {}
with open(f"{A}/windows.fq") as f:
    for k, line in enumerate(f):
        if k % 4 == 0:
            name, _, fl = line[1:].split()
            flank[name] = int(fl.split("=")[1])
events = {r["event"]: r for r in csv.DictReader(open(f"{A}/events.tsv"), delimiter="\t")}
aln = {}
for line in open(f"{O}/windows.d46.gam.json"):
    a = json.loads(line)
    aln[a["name"]] = a


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


def on_ref(n):
    return 0 < n < len(rp.visits) and int(rp.visits[n]) == 1


rows, cls = [], Counter()
for eid, ev in events.items():
    t, n = aln.get(f"{eid}_T"), aln.get(f"{eid}_N")
    F, L = flank[f"{eid}_T"], max(1, len(ev["alt"]))
    row = dict(event=eid, contig=ev["contig"], pos=ev["pos"], ref=ev["ref"], alt=ev["alt"],
               grch38_records=ev["grch38_records"], cls="", edit_node="", edit_rev="", edit_offset="", edit="",
               background_on_ref="", mapq_T=(t or {}).get("mapping_quality", 0), mapq_N=(n or {}).get("mapping_quality", 0))
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
with open(f"{O}/classes_d46.tsv", "w") as f:
    w = csv.DictWriter(f, fieldnames=list(rows[0]), delimiter="\t")
    w.writeheader()
    w.writerows(rows)
print("events", len(rows), dict(sorted(cls.items())))

b = [r for r in rows if r["cls"] == "B"]
with GraphIndex(f"{G}/hprc-v2.1-mc-grch38.d46.graph_index.sqlite") as g:
    seqs = {k: v["sequence"] for k, v in g.get_nodes({int(r["edit_node"]) for r in b}).items()}
out = []
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
    out.append(dict(event=r["event"], contig=r["contig"], pos=r["pos"], ref=r["ref"], alt=r["alt"],
                    grch38_records=r["grch38_records"], normal_offref_at_site=r["background_on_ref"] is False,
                    node=node, start=start, kind=kind, key=f"{node}:{start}:{kind}:{ref}>{alt}"))
with open(f"{O}/offref_somatic_d46.tsv", "w") as f:
    w = csv.DictWriter(f, fieldnames=list(out[0]), delimiter="\t")
    w.writeheader()
    w.writerows(out)
d9 = {r["event"]: r for r in csv.DictReader(open(f"{A}/classes.tsv"), delimiter="\t")}
print("off-reference somatic (B) on d46:", len(out), dict(Counter(e["kind"] for e in out)),
      "| chr1 SNP:", sum(e["kind"] == "SNP" and e["contig"].startswith("chr1_") for e in out))
print("d9 class -> d46 class of the d9 B events:",
      dict(Counter(r["cls"] for r in rows if d9.get(r["event"], {}).get("cls") == "B")))
print("d46 B events by d9 class:", dict(Counter(d9.get(r["event"], {}).get("cls") for r in b)))
