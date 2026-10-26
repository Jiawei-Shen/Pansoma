"""E9j step 5: read support of the COLO829T off-reference somatic SNV keys, the E9h method (E9h_offref_key_reads:
build_queries.py + classify.py, h50 strict rule) with the COLO829T assemblies:
  normal  COLO829BL donor-specific assembly, verkko 2.1 hapX / hapY (targets N_hapX, N_hapY)
  tumor   COLO829T hifiasm 0.25.0 hap1 / hap2 (targets T_hap1, T_hap2; qfu: Fiber-seq HiFi + ONT D003 ultra-long reads,
          so not independent of the Fiber-seq and ONT tensors)
  queries   the A1-read consensus of each key's tensor (per column the most frequent base of >= 2 covering A1 rows,
            the contiguous run around the site column, up to 50 columns each side; dropped when the site consensus is
            not the key's ALT or the query is < 25 bp) + the events' N / T windows (s1) -> queries.fa, keys.tsv,
            windows.tsv
  classify  minimap2 -ax sr (s5_map.sh) hits under the strict rule (>= 90 % aligned, NM <= 1, no indel, the site base
            and its neighbours matched); a key is tumor_only (hits the tumor assembly, not the normal: read-supported
            somatic), normal (hits the normal at the event's locus, or anywhere when the tumor hit is not at the locus),
            ambiguous (no query, or a normal hit only away from the locus while the tumor hit is at it) or neither;
            joined with s4's scores -> results/key_classes.tsv, tables.md, summary.json (next to this script)
usage: python s5_reads.py {queries|classify} DATA_DIR
"""
import csv
import json
import re
import sys
from collections import Counter, defaultdict
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
TENSORS = Path("/scratch/jshen/data/pansoma_v2_tensors")
SETS = [  # (set, graph, run under pansoma_net_v2_runs/reviewer_cpu_20261013)
    ("COLO829T_Illumina", "d9", "COLO829_Illumina_SNV_d9_normbq_sweep_e06"),
    ("COLO829T_fiberseq", "d9", "COLO829_fiberseq_SNV_7fc7396_ddp4_b1024"),
    ("COLO829T_ONT", "d9", "COLO829_ONT_SNV_7fc7396_ddp4_b1024"),
    ("COLO829T_Illumina_d46", "d46", "COLO829_Illumina_SNV_d46_7fc7396"),
    ("COLO829T_fiberseq_d46", "d46", "COLO829_fiberseq_SNV_d46_7fc7396"),
    ("COLO829T_ONT_d46", "d46", "COLO829_ONT_SNV_d46_7fc7396"),
]
TARGETS = ("N_hapX", "N_hapY", "T_hap1", "T_hap2")
HALF, MIN_LEN = 50, 25
BASES = {1: "A", 2: "C", 3: "G", 4: "T", 5: "N"}
PREFIX = '{"candidate_id": "'
LOCUS_PAD, OWN_PAD, WIN_NM = 100, 200, 5
CIG = re.compile(r"(\d+)([MIDNSHP=X])")
CS = re.compile(r"(:\d+|\*[a-z][a-z]|[+-][a-z]+|~[a-z]{2}\d+[a-z]{2})")


def consensus(x, rows):
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
    found = {}
    for i in range(1, 23):
        with open(set_dir / f"chr{i}_variant_summary.ndjson") as f:
            for line in f:
                cid = line[len(PREFIX):line.index('"', len(PREFIX))]
                if cid in keys:
                    s = json.loads(line)
                    found[cid] = dict(chrom=s["chrom"], shard=s["shard_file"], row=int(s["index_within_shard"]),
                                      groups=s["row_groups"], anchor=int(s["anchor_column"]), layout=s["site_layout"],
                                      af=float(s["af"]), alt_count=int(s["alt_count"]), coverage=int(s["coverage"]))
    return found


def queries(data):
    fa = open(data / "queries.fa.tmp", "w")
    cols = ["set", "graph", "run", "key", "alt", "event", "chrom", "af", "alt_count", "coverage", "n_a1", "q_id", "q_len",
            "q_site", "q_status"]
    kt = open(data / "keys.tsv.tmp", "w")
    kt.write("\t".join(cols) + "\n")
    used = set()
    for si, (name, graph, run) in enumerate(SETS):
        rows = [r for r in csv.DictReader(open(data / f"offref_somatic_{graph}.tsv"), delimiter="\t") if r["kind"] == "SNP"]
        keys = {}
        for r in rows:
            keys.setdefault(r["key"], r)
        found = find_tensors(TENSORS / name / "tensors" / "SNV", set(keys))
        shards, stat = {}, Counter()
        for ki, key in enumerate(sorted(found)):
            t, ev = found[key], keys[key]
            used.add(ev["event"])
            lay = t["layout"]
            assert lay["insertion_slots"] == 0 and lay["span"] == 1 and t["anchor"] == 50, (name, key, lay)
            if t["shard"] not in shards:
                shards[t["shard"]] = np.load(TENSORS / name / "tensors" / "SNV" / t["shard"], mmap_mode="r")
            x = np.asarray(shards[t["shard"]][t["row"]])
            a1 = next(g for g in t["groups"] if g["allele"] == "A1")
            rws = np.arange(a1["start_row"], a1["end_row"])
            alt = key.split(">")[1]
            seq, site = query(consensus(x, rws), t["anchor"], HALF)
            status = ("no_site_consensus" if site is None else "site_consensus_not_alt" if seq[site] != alt
                      else "short" if len(seq) < MIN_LEN else "ok")
            stat[status] += 1
            qid = f"s{si}k{ki:04d}" if status == "ok" else ""
            if qid:
                fa.write(f">{qid}\n{seq}\n")
            rec = dict(set=name, graph=graph, run=run, key=key, alt=alt, event=ev["event"], chrom=t["chrom"],
                       af=f"{t['af']:.4f}", alt_count=t["alt_count"], coverage=t["coverage"], n_a1=len(rws), q_id=qid,
                       q_len=len(seq), q_site="" if site is None else site, q_status=status)
            kt.write("\t".join(str(rec[c]) for c in cols) + "\n")
        print(f"{name}: {len(found)} of {len(keys)} SNP keys have a tensor; queries {dict(stat)}", flush=True)
    events = {r["event"]: r for r in csv.DictReader(open(data / "events.tsv"), delimiter="\t")}
    wt = open(data / "windows.tsv.tmp", "w")
    wt.write("qid\tevent\twhich\thap\tcontig\tpos0\tflank\tlen\n")
    with open(data / "windows.fq") as f:
        while True:
            head = f.readline()
            if not head:
                break
            seq = f.readline().strip()
            f.readline()
            f.readline()
            name, _, fl = head[1:].split()
            eid, which = name.split("_")
            if eid not in used:
                continue
            e = events[eid]
            fa.write(f">w_{eid}_{which}\n{seq}\n")
            wt.write(f"w_{eid}_{which}\t{eid}\t{which}\t{e['hap']}\t{e['contig']}\t{e['pos0']}\t{fl.split('=')[1]}\t{len(seq)}\n")
    for f in (fa, kt, wt):
        f.close()
    for stem in ("queries.fa", "keys.tsv", "windows.tsv"):
        (data / f"{stem}.tmp").replace(data / stem)
    print(f"events {len(used)}; wrote queries.fa, keys.tsv, windows.tsv")


def site_exact(cs, start, lo, hi):
    p = start
    if lo < start:
        return False
    for op in CS.findall(cs):
        t = op[0]
        if t == ":":
            p += int(op[1:])
        elif t == "*":
            if lo <= p < hi:
                return False
            p += 1
        elif t == "+":
            n = len(op) - 1
            if p < hi and p + n > lo:
                return False
            p += n
        elif t == "-":
            if lo < p < hi:
                return False
        elif t == "~":
            return False
        if p >= hi and t == ":":
            return True
    return p >= hi


def read_sam(path):
    for line in open(path):
        if line[0] == "@":
            continue
        f = line.rstrip("\n").split("\t")
        flag = int(f[1])
        if flag & 4:
            continue
        cig = CIG.findall(f[5])
        lead = int(cig[0][0]) if cig and cig[0][1] in "SH" else 0
        clips = sum(int(n) for n, o in cig if o in "SH")
        rlen = sum(int(n) for n, o in cig if o in "MDN=X")
        tags = {t.split(":", 2)[0]: t.split(":", 2)[2] for t in f[11:]}
        start = int(f[3]) - 1
        yield f[0], bool(flag & 16), f[2], start, start + rlen, lead, clips, int(tags.get("NM", 99)), tags.get("cs", "")


def strict(cs, rev, lead, clips, nm, qlen, site):
    if qlen - clips < 0.9 * qlen or nm > 1 or not cs or "+" in cs or "-" in cs:
        return False
    lo, hi = max(0, site - 1), min(qlen, site + 2)
    if rev:
        lo, hi = qlen - hi, qlen - lo
    return site_exact(cs, lead, lo, hi)


def classify(data):
    keys = list(csv.DictReader(open(data / "keys.tsv"), delimiter="\t"))
    windows = {r["qid"]: r for r in csv.DictReader(open(data / "windows.tsv"), delimiter="\t")}
    qinfo = {k["q_id"]: (int(k["q_len"]), int(k["q_site"])) for k in keys if k["q_id"]}
    hits, loci = defaultdict(list), defaultdict(list)
    for t in TARGETS:
        for q, rev, contig, s, e, lead, clips, nm, cs in read_sam(data / f"{t}.sam"):
            if q in windows:
                L = int(windows[q]["len"])
                if L - clips >= 0.9 * L and nm <= WIN_NM:
                    loci[(windows[q]["event"], t)].append((contig, s - LOCUS_PAD, e + LOCUS_PAD))
                continue
            if strict(cs, rev, lead, clips, nm, *qinfo[q]):
                hits[q].append((t, contig, s, e))
    for w in windows.values():
        if w["which"] == "N":
            p = int(w["pos0"])
            loci[(w["event"], "N_" + w["hap"])].append((w["contig"], p - OWN_PAD, p + OWN_PAD))

    def at_locus(event, h):
        t, contig, s, e = h
        return any(c == contig and s < hi and e > lo for c, lo, hi in loci.get((event, t), []))

    scores = {}
    for name, graph, run in SETS:
        p = data / f"scores_{name}.json"
        if p.exists():
            j = json.loads(p.read_text())
            for r in j["keys"]:
                scores[(name, r["key"])] = dict(r, threshold=j["threshold"])
    out = []
    for k in keys:
        row = dict(set=k["set"], graph=k["graph"], run=k["run"], key=k["key"], event=k["event"], chrom=k["chrom"],
                   af=k["af"], n_a1=k["n_a1"], q_status=k["q_status"])
        if not k["q_id"]:
            row.update(cls="ambiguous", why=f"no_query:{k['q_status']}")
        else:
            hs = hits.get(k["q_id"], [])
            nh = [h for h in hs if h[0][0] == "N"]
            th = [h for h in hs if h[0][0] == "T"]
            nloc = any(at_locus(k["event"], h) for h in nh)
            tloc = any(at_locus(k["event"], h) for h in th)
            if nh:
                cls, why = ("ambiguous", "paralog") if (not nloc and tloc) else ("normal", "at_locus" if nloc else "elsewhere")
            elif th:
                cls, why = "tumor_only", "at_locus" if tloc else "elsewhere"
            else:
                cls, why = "neither", ""
            row.update(cls=cls, why=why, n_hits_normal=len(nh), n_hits_tumor=len(th))
        sc = scores.get((k["set"], k["key"]))
        if sc:
            t = sc["threshold"]
            row.update(threshold=round(t, 5), label=sc["label"], reason=sc["reason"], p_base=sc["p_base"],
                       p_rescue=sc["p_rescue"], p_all_H=sc["p_all_H"], base_t=sc["p_base"] >= t,
                       rescue_argmax=sc["argmax_rescue"] == 1, all_H_t=sc["p_all_H"] >= t)
        out.append(row)
    (HERE / "results").mkdir(exist_ok=True)
    cols = []
    for r in out:
        for c in r:
            if c not in cols:
                cols.append(c)
    with open(HERE / "results" / "key_classes.tsv", "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=cols, delimiter="\t", lineterminator="\n")
        w.writeheader()
        w.writerows(out)
    summary = {}
    L = ["# E9j: COLO829T assembly off-reference somatic SNV keys (h50 strict read check)", "",
         "| set | keys (chr1) | tumor-only, read-supported (chr1) | base / rescue / all_H on them | "
         "rescue calls on all keys: tumor-only / all | classes (all keys) |", "|---|---|---|---|---|---|"]
    for name, graph, run in SETS:
        rs = [r for r in out if r["set"] == name]
        if not rs:
            L.append(f"| `{name}` | 0 | - | - | - | - |")
            continue
        n1 = sum(r["chrom"] == "chr1" for r in rs)
        to = [r for r in rs if r["cls"] == "tumor_only"]
        to1 = sum(r["chrom"] == "chr1" for r in to)
        scored = all("p_base" in r for r in rs)
        b = sum(bool(r.get("base_t")) for r in to)
        rr = sum(bool(r.get("rescue_argmax")) for r in to)
        ah = sum(bool(r.get("all_H_t")) for r in to)
        rc_all = sum(bool(r.get("rescue_argmax")) for r in rs)
        rc_to = rr
        cls = Counter(r["cls"] for r in rs)
        summary[name] = dict(graph=graph, run=run, keys=len(rs), keys_chr1=n1, tumor_only=len(to), tumor_only_chr1=to1,
                             base_t_on_tumor_only=b, rescue_on_tumor_only=rr, all_H_on_tumor_only=ah,
                             rescue_calls_all=rc_all, rescue_calls_tumor_only=rc_to, classes=dict(cls), scored=scored)
        sc = f"{b} / {rr} / {ah}" if scored else "(not scored)"
        L.append(f"| `{name}` ({graph}) | {len(rs)} ({n1}) | {len(to)} ({to1}) | {sc} | "
                 f"{rc_to} / {rc_all if scored else '-'} | {dict(cls)} |")
    L += ["", "base = the set's in-sample SNV model at its PASS threshold (validation recall 0.9); rescue = predict's "
          "off-reference rescue (site +- 10 raised to the every-haplotype code, argmax); all_H = ch6 set to H everywhere, "
          "at the threshold. All autosomes; chr1 in brackets."]
    (HERE / "tables.md").write_text("\n".join(L) + "\n")
    (HERE / "summary.json").write_text(json.dumps(summary, indent=1) + "\n")
    print("\n".join(L))


if __name__ == "__main__":
    {"queries": queries, "classify": classify}[sys.argv[1]](Path(sys.argv[2]))
