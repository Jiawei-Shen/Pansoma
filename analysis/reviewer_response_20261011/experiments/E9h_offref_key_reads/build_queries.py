"""E9h step 1: A1-read consensus queries of the assembly off-reference somatic SNV keys, and the events' assembly
windows, for the assembly check (step 2, map.sh) and the classification (step 3, classify.py).

Keys (SNP keys of an assembly off-reference somatic event that have a tensor in the set):
  d9   pansoma_net_v2_runs/assembly_offref_20260929/offref_somatic.tsv: HG008T_Illumina (no E1.1 scores; reproduction
       of offref_negatives_20260929), HG008T_PacBio, HG008T_ONT
  d46  reviewer_experiments_20261013/E1.1_ch6_counterfactual/asm_d46/offref_somatic_d46.tsv: HG008T_Illumina_d46,
       HG008T_PacBio_d46, HG008T_ONT_d46
The tensor of a key is found by its candidate_id in <set>/tensors/SNV/chr*_variant_summary.ndjson (every autosome).
For the five sets that E1.1 scored, the keys found here must equal the keys of E1.1's results/assembly_keys.json.

Query (the method of pansoma_net_v2_runs/offref_negatives_20260929/v1_queries.py, adapted, the original unchanged):
the A1 rows of the tensor; per column the most frequent read base (channel 0; 1-4 = ACGT, 5 = N, 6 = gap) over the A1
rows that cover it (channel 4 > 0), when >= 2 rows cover it; the query is the contiguous run of such columns around the
site column (anchor_column, 50) up to HALF columns on each side; gap columns add no base. Dropped (no query) when the
site column has no base consensus, the site consensus is not the key's ALT, or the query is < 25 bp. HALF = 30 is the
earlier method (61 bp at most); HALF = 50 is the whole tensor width (101 bp at most).

Windows: the event's normal-haplotype window (N) and the same with the somatic allele (T), 100 bp flanks, from
assembly_offref_20260929/windows.fq (s1_windows.py), mapped alongside to find the event's locus in every assembly
haplotype.

usage: python build_queries.py OUT_DIR   (writes OUT_DIR/queries.fa, keys.tsv, windows.tsv)"""
import csv
import json
import sys
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
E11 = HERE.parent / "E1.1_ch6_counterfactual" / "results"
A9 = Path("/scratch/jshen/data/pansoma_net_v2_runs/assembly_offref_20260929")
A46 = Path("/scratch/jshen/data/reviewer_experiments_20261013/E1.1_ch6_counterfactual/asm_d46/offref_somatic_d46.tsv")
TENSORS = Path("/scratch/jshen/data/pansoma_v2_tensors")
SETS = [  # (set, graph, E1.1 run or None)
    ("HG008T_Illumina", "d9", None),
    ("HG008T_PacBio", "d9", "HG008_PacBio_SNV_7fc7396"),
    ("HG008T_ONT", "d9", "HG008_ONT_SNV_7fc7396"),
    ("HG008T_Illumina_d46", "d46", "HG008_Illumina_SNV_d46_7fc7396"),
    ("HG008T_PacBio_d46", "d46", "HG008_PacBio_SNV_d46_7fc7396"),
    ("HG008T_ONT_d46", "d46", "HG008_ONT_SNV_d46_7fc7396"),
]
HALVES = (30, 50)
MIN_LEN = 25
BASES = {1: "A", 2: "C", 3: "G", 4: "T", 5: "N"}
PREFIX = '{"candidate_id": "'


def consensus(x, rows):
    """Per column: the most frequent read base code of the covered A1 rows (>= 2 rows), else None (v1_queries.py)."""
    base, cov = x[0, rows], x[4, rows] > 0
    cons = []
    for c in range(x.shape[2]):
        v = base[cov[:, c], c]
        if len(v) < 2:
            cons.append(None)
            continue
        vals, counts = np.unique(v, return_counts=True)
        cons.append(int(vals[counts.argmax()]))
    return cons


def query(cons, anchor, half):
    """The contiguous consensus run around the anchor column, at most `half` columns each side (v1_queries.py)."""
    lo = anchor
    while lo > max(0, anchor - half) and cons[lo - 1] is not None:
        lo -= 1
    hi = anchor
    while hi < min(len(cons) - 1, anchor + half) and cons[hi + 1] is not None:
        hi += 1
    seq, site = [], None
    for c in range(lo, hi + 1):
        if cons[c] in BASES:
            if c == anchor:
                site = len(seq)
            seq.append(BASES[cons[c]])
    return "".join(seq), site


def find_tensors(set_dir, keys):
    """candidate_id -> its variant_summary record (the fields used here), over every autosome."""
    found = {}
    for i in range(1, 23):
        p = set_dir / f"chr{i}_variant_summary.ndjson"
        with open(p) as f:
            for line in f:
                if not line.startswith(PREFIX):
                    raise SystemExit(f"unexpected summary line start in {p}: {line[:40]!r}")
                cid = line[len(PREFIX):line.index('"', len(PREFIX))]
                if cid in keys:
                    s = json.loads(line)
                    assert cid not in found, cid
                    found[cid] = dict(chrom=s["chrom"], shard=s["shard_file"], row=int(s["index_within_shard"]),
                                      groups=s["row_groups"], anchor=int(s["anchor_column"]),
                                      layout=s["site_layout"], af=float(s["af"]), alt_count=int(s["alt_count"]),
                                      coverage=int(s["coverage"]))
    return found


def main():
    out = Path(sys.argv[1])
    out.mkdir(parents=True, exist_ok=True)
    e11 = json.loads((E11 / "assembly_keys.json").read_text())
    summ = json.loads((E11 / "summary.json").read_text())["runs"]
    tables = {"d9": A9 / "offref_somatic.tsv", "d46": A46}
    events = {}
    for graph, p in tables.items():
        rows = [r for r in csv.DictReader(open(p), delimiter="\t") if r["kind"] == "SNP"]
        events[graph] = {r["key"]: r for r in rows}   # E1.1 convention (one row per key; no duplicate SNP keys)
        assert len(events[graph]) == len(rows), graph
    fa = open(out / "queries.fa.tmp", "w")
    kt = open(out / "keys.tsv.tmp", "w")
    cols = ["set", "graph", "run", "key", "alt", "event", "chrom", "shard", "row", "af", "alt_count", "coverage",
            "n_a1", "n_a1_site", "site_alt_frac"]
    for h in HALVES:
        cols += [f"q{h}_id", f"q{h}_len", f"q{h}_site", f"q{h}_status"]
    kt.write("\t".join(cols) + "\n")
    used_events = set()
    for si, (name, graph, run) in enumerate(SETS):
        if run is not None:
            assert summ[run]["set"] == name, (run, summ[run]["set"], name)
        keys = events[graph]
        found = find_tensors(TENSORS / name / "tensors" / "SNV", set(keys))
        if run is not None:
            e11_keys = {r["key"] for r in e11[run]}
            if set(found) != e11_keys:
                raise SystemExit(f"{name}: {len(found)} keys found, E1.1 {len(e11_keys)}; "
                                 f"only here {sorted(set(found) - e11_keys)[:5]}, only E1.1 {sorted(e11_keys - set(found))[:5]}")
        shards = {}
        stat = {}
        for ki, key in enumerate(sorted(found)):
            t = found[key]
            ev = keys[key]
            used_events.add(ev["event"])
            lay = t["layout"]
            assert lay["insertion_slots"] == 0 and lay["span"] == 1 and t["anchor"] == 50, (name, key, lay)
            if t["shard"] not in shards:
                shards[t["shard"]] = np.load(TENSORS / name / "tensors" / "SNV" / t["shard"], mmap_mode="r")
            x = np.asarray(shards[t["shard"]][t["row"]])
            a1 = next(g for g in t["groups"] if g["allele"] == "A1")
            rows = np.arange(a1["start_row"], a1["end_row"])
            alt = key.split(">")[1]
            site_bases = x[0, rows, t["anchor"]][x[4, rows, t["anchor"]] > 0]
            n_site = int(len(site_bases))
            alt_frac = float(np.mean([BASES.get(int(b)) == alt for b in site_bases])) if n_site else 0.0
            cons = consensus(x, rows)
            rec = dict(set=name, graph=graph, run=run or "", key=key, alt=alt, event=ev["event"], chrom=t["chrom"],
                       shard=t["shard"], row=t["row"], af=f"{t['af']:.4f}", alt_count=t["alt_count"],
                       coverage=t["coverage"], n_a1=len(rows), n_a1_site=n_site, site_alt_frac=f"{alt_frac:.3f}")
            for h in HALVES:
                seq, site = query(cons, t["anchor"], h)
                if site is None:
                    status = "no_site_consensus"
                elif seq[site] != alt:
                    status = "site_consensus_not_alt"
                elif len(seq) < MIN_LEN:
                    status = "short"
                else:
                    status = "ok"
                stat[(h, status)] = stat.get((h, status), 0) + 1
                qid = f"s{si}k{ki:03d}h{h}" if status == "ok" else ""
                if qid:
                    fa.write(f">{qid}\n{seq}\n")
                rec.update({f"q{h}_id": qid, f"q{h}_len": len(seq), f"q{h}_site": "" if site is None else site,
                            f"q{h}_status": status})
            kt.write("\t".join(str(rec[c]) for c in cols) + "\n")
        print(f"{name} ({graph}): {len(found)} keys with a tensor of {len(keys)} SNP keys; "
              + ", ".join(f"h{h} {s} {n}" for (h, s), n in sorted(stat.items())), flush=True)
    # event windows (N = normal haplotype, T = somatic allele), from windows.fq
    ev_rows = {r["event"]: r for r in csv.DictReader(open(A9 / "events.tsv"), delimiter="\t")}
    wt = open(out / "windows.tsv.tmp", "w")
    wt.write("qid\tevent\twhich\tcontig\tpos\tref\talt\tflank\tlen\n")
    n = 0
    with open(A9 / "windows.fq") as f:
        while True:
            head = f.readline()
            if not head:
                break
            seq = f.readline().strip()
            f.readline()
            f.readline()
            name, _, fl = head[1:].split()
            eid, which = name.split("_")
            if eid not in used_events:
                continue
            qid = f"w_{eid}_{which}"
            fa.write(f">{qid}\n{seq}\n")
            e = ev_rows[eid]
            wt.write(f"{qid}\t{eid}\t{which}\t{e['contig']}\t{e['pos']}\t{e['ref']}\t{e['alt']}\t{fl.split('=')[1]}\t{len(seq)}\n")
            n += 1
    assert n == 2 * len(used_events), (n, len(used_events))
    for f in (fa, kt, wt):
        f.close()
    for stem in ("queries.fa", "keys.tsv", "windows.tsv"):
        (out / f"{stem}.tmp").replace(out / stem)
    print(f"events {len(used_events)}, windows {n}; wrote {out}/queries.fa, keys.tsv, windows.tsv")


if __name__ == "__main__":
    main()
