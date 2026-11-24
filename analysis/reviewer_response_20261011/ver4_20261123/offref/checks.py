"""ver4 offref: consistency checks of the inputs and of the read-evidence step. Writes checks.md.

1. Truth counts. For each of the nine call sets, TP + FN == baseline in eval_bed and eval_nobed (raw and pon); the
   three platforms of a tumor have the same truth count; the in-BED truth count equals the number of SNVs in the E18
   chr1 truth spectrum used for the cosine similarities of Figure 3.3.
2. Call counts. Graph-coordinate PASS calls: unplaced VCF == report.json == E9c run.json == rows of per_call.tsv.gz ==
   rows of the E9i read tables == rows of the E9a projection; expected 77 / 143 / 90, 110 / 79 / 111, 76 / 106 / 143.
3. E9i against the earlier runs. A call (same candidate key, same ALT and REF 21-mer) that was also counted in the ver3
   runs (HG008T / COLO829T single-tumor no-rescue runs, HCC1395T six-tumor no-rescue runs) must have the same read
   counts: same BAMs, same locus (for HG008T this also checks the new remote fetch of HG008-N-D, for HCC1395T the
   local full normal BAM against the earlier remote subset).
4. HG008-N-D fetch: every region of regions.bed was fetched (log), the subset BAM has an index.
usage: /wanglab/jshen/anaconda3/bin/python checks.py
"""
import collections
import csv
import glob
import gzip
import json
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
A = HERE.parents[1]
T = Path("/scratch/jshen/Github/Pansoma/tmp")
X4 = Path("/scratch/jshen/data/reviewer_experiments_20261013/ver4_20261123/offref")
EXPECTED_U = {"HG008T": (77, 143, 90), "COLO829T": (110, 79, 111), "HCC1395T": (76, 106, 143)}
OLD_READS = {"HG008T": [T / "e9i_hg008t_20261104/out"], "COLO829T": [T / "assessment_20261030/e9c_norescue/E9i_out"],
             "HCC1395T": [T / "r33_hcc1395_20261104/e9i/out"]}
COUNTS = ("normal_alt", "normal_ref", "normal_reads", "tumor_alt", "tumor_ref", "tumor_reads", "status")
csv.field_size_limit(sys.maxsize)
out = []


def say(s=""):
    print(s)
    out.append(s)


def count_pass(path):
    return sum(1 for l in gzip.open(path, "rt") if l[0] != "#" and l.split("\t", 8)[6] == "PASS")


runs = [l.split() for l in open(HERE / "runs.txt")]
# 1. truth
say("## 1. Truth counts (TP + FN) per call set")
say()
say("| call set | in BED: baseline, raw TP+FN, after-filter TP+FN | no BED: baseline, raw TP+FN, after-filter TP+FN | E18 truth spectrum n |")
say("|---|---|---|---|")
e18 = {}
for row in csv.DictReader(open(A / "experiments/E18_sbs96/sbs96.tsv"), delimiter="\t"):
    if row["kind"] == "truth":
        e18[row["sample"]] = sum(int(v) for k, v in row.items() if "[" in k)
truth = collections.defaultdict(set)
for name, vdir, tumor, plat in runs:
    cell = {}
    for ev in ("eval_bed", "eval_nobed"):
        r = json.load(open(f"{vdir}/{ev}/report.json"))
        b = {k: r[k]["baseline"] for k in ("raw", "pon")}
        s = {k: r[k]["at_threshold"]["tp_baseline"] + r[k]["at_threshold"]["fn"] for k in ("raw", "pon")}
        assert b["raw"] == b["pon"] == s["raw"] == s["pon"], (name, ev, b, s)
        cell[ev] = f"{b['raw']}, {s['raw']}, {s['pon']}"
        truth[tumor].add((ev, b["raw"]))
        if ev == "eval_bed":
            assert b["raw"] == e18[tumor], (name, b["raw"], e18[tumor])
    say(f"| {tumor} {plat.replace('_', ' ')} | {cell['eval_bed']} | {cell['eval_nobed']} | {e18[tumor]} |")
for t, s in truth.items():
    assert len(s) == 2, (t, s)
say()
say("All equal within a tumor (three platforms) and equal to the E18 truth spectrum counts: OK.")
say()
# 2. calls
say("## 2. Graph-coordinate PASS call counts per source")
say()
say("| call set | unplaced VCF PASS | report.json (BED / no BED) | E9c run.json | per_call rows | E9i rows | E9a projection PASS rows | expected |")
say("|---|---|---|---|---|---|---|---|")
pc = collections.Counter()
for r in csv.DictReader(gzip.open(HERE / "EXP/per_call.tsv.gz", "rt"), delimiter="\t"):
    if r["group"] == "unplaced_pass":
        pc[r["run"]] += 1
new_reads = collections.defaultdict(dict)
for f in sorted(glob.glob(str(X4 / "e9i/out/*.tsv"))):
    run = Path(f).name.split(".")[0]
    for r in csv.DictReader(open(f), delimiter="\t"):
        new_reads[run][(r["group"], r["candidate_id"])] = r
seen = collections.defaultdict(list)
for name, vdir, tumor, plat in runs:
    (u,) = Path(vdir).glob("*.linear.unplaced.vcf.gz")
    n = count_pass(u)
    rb = json.load(open(f"{vdir}/eval_bed/report.json"))["unplaced_pass"]
    rn = json.load(open(f"{vdir}/eval_nobed/report.json"))["unplaced_pass"]
    rj = json.load(open(X4 / "DATA/runs" / name / "run.json"))
    e9i = sum(1 for k in new_reads[name] if k[0] == "unplaced_pass")
    proj = sum(1 for r in csv.DictReader(open(X4 / "r34/e9a" / name / "projected.tsv"), delimiter="\t") if r["filter"] == "PASS")
    exp = EXPECTED_U[tumor][len(seen[tumor])]
    seen[tumor].append(name)
    assert n == rb == rn == rj["unplaced_pass"] == rj["unplaced_pass_processed"] == pc[name] == e9i == proj == exp, \
        (name, n, rb, rn, rj["unplaced_pass"], pc[name], e9i, proj, exp)
    say(f"| {tumor} {plat.replace('_', ' ')} | {n} | {rb} / {rn} | {rj['unplaced_pass']} | {pc[name]} | {e9i} | {proj} | {exp} |")
say()
say(f"All equal: OK. Total {sum(pc[n] for n, *_ in runs)} calls.")
say()
# 3. E9i vs earlier runs
say("## 3. E9i read counts against the ver3 runs (calls counted in both, same candidate key and 21-mers)")
say()
say("| tumor | calls in both (any group) | identical counts and status | different | of them graph-coordinate PASS calls: in both / identical |")
say("|---|---|---|---|---|")
for tumor, dirs in OLD_READS.items():
    old = {}
    for d in dirs:
        for f in sorted(glob.glob(str(d / "*.tsv"))):
            for r in csv.DictReader(open(f), delimiter="\t"):
                if r.get("alt_kmer"):
                    old.setdefault((r["candidate_id"], r["alt_kmer"], r["ref_kmer"], r["locus"]), r)
    both = same = ub = us = 0
    diff = []
    for name, vdir, t, plat in runs:
        if t != tumor:
            continue
        for (g, cid), r in new_reads[name].items():
            o = old.get((cid, r.get("alt_kmer"), r.get("ref_kmer"), r["locus"]))
            if not o:
                continue
            both += 1
            ok = all(str(o[c]) == str(r[c]) for c in COUNTS)
            same += ok
            if g == "unplaced_pass":
                ub += 1
                us += ok
            if not ok:
                diff.append((name, g, cid, [(c, o[c], r[c]) for c in COUNTS if str(o[c]) != str(r[c])]))
    say(f"| {tumor} | {both} | {same} | {both - same} | {ub} / {us} |")
    for d in diff[:10]:
        say(f"|  | different: {d} | | | |")
say()
# 4. fetch
say("## 4. HG008-N-D remote fetch")
say()
reg = sum(1 for _ in open(X4 / "e9i/HG008T/regions.bed"))
logs = sorted(glob.glob(str(T / "ver4_20261123/offref/logs/e9i_fetch.*.out")))
done = [l.strip() for f in logs for l in open(f) if l.startswith("done ")]
retries = sum(1 for f in logs for l in open(f) if l.startswith("retry"))
say(f"regions.bed: {reg} regions; fetch log: `{done[-1] if done else 'NOT DONE'}`; retries logged: {retries}; "
    f"index present: {(X4 / 'e9i/HG008T/normal_subset.bam.bai').exists()}.")
assert done and done[-1].split()[1] == str(reg)
say()
# 5. the call sets: model, PASS rule, test set, rescue off
say("## 5. Call sets (from the VCF headers and report.json)")
say()
say("| call set | tensor set | checkpoint | PASS threshold | call rule | records rescored / called by the off-reference rescue | OFFREF_P_SOMATIC records in the graph-coordinate VCF |")
say("|---|---|---|---|---|---|---|")
for name, vdir, tumor, plat in runs:
    (u,) = Path(vdir).glob("*.linear.unplaced.vcf.gz")
    hdr, off = {}, 0
    for l in gzip.open(u, "rt"):
        if l.startswith("##pansoma_"):
            k, v = l[2:].rstrip("\n").split("=", 1)
            hdr[k] = v
        elif l[0] != "#":
            off += "OFFREF_P_SOMATIC=" in l
    g = json.load(open(f"{vdir}/eval_nobed/report.json"))["steps"]["graph_vcf"]
    assert off == 0 and g["offref_rescored"] == 0 and g["offref_called"] == 0, (name, off, g["offref_rescored"])
    assert "hprc-v2.1-mc-grch38.d46" in hdr["pansoma_graph_index"] and hdr["pansoma_kind"] == "SNV" and hdr["pansoma_chromosomes"] == "chr1"
    ck = Path(hdr["pansoma_checkpoint"])
    say(f"| {tumor} {plat.replace('_', ' ')} | {Path(hdr['pansoma_tensors']).parent.name} | {ck.parent.parent.name}/{ck.parent.name}/{ck.name} | "
        f"{float(hdr['pansoma_threshold']):.4f} | {hdr['pansoma_call'].split('(', 1)[1][:-1]} | {g['offref_rescored']} / {g['offref_called']} | {off} |")
say()
say("All nine: HPRC v2.1 (d46) graph index, kind SNV, chr1, no off-reference rescue: OK.")
say()
# 6. assembly step reproduced on an old run
say("## 6. E9c + E9d reproduced on a ver3 run with this batch's code (HCC1395T Illumina six-tumor no-rescue run of 2026-11-04)")
say()
rn = "mc6__HCC1395T_Illumina__norescue__vcf_chr1"
old_calls = T / "r33_hcc1395_20261104/DATA/runs" / rn / "calls.tsv"
new_calls = X4 / "repro/DATA/runs" / rn / "calls.tsv"
same_calls = old_calls.read_bytes() == new_calls.read_bytes()
say(f"E9c calls.tsv (13 graph-coordinate PASS calls + controls + gold; {sum(1 for _ in open(new_calls)) - 1} rows): "
    f"{'byte-identical to the 2026-11-04 file' if same_calls else 'DIFFERENT'}.")
assert same_calls
rp = HERE / "repro/EXP/per_call.tsv.gz"
if rp.exists():
    load = lambda p: {(r["group"], r["candidate_id"]): (r["asm"], r["asm_flag"], r["category"], r["pop"])
                      for r in csv.DictReader(gzip.open(p, "rt"), delimiter="\t") if r["run"] == rn}
    o, n = load(T / "r33_hcc1395_20261104/EXP/per_call.tsv.gz"), load(rp)
    d = [k for k in o if o[k] != n.get(k)]
    up = collections.Counter(v[0] for k, v in n.items() if k[0] == "unplaced_pass")
    say(f"E9d (queries, minimap2 on the five HCC1395 assemblies, classify, collect): {len(o)} rows in the old per_call.tsv.gz, "
        f"{len(n)} in the reproduction, {len(o) - len(d)} with the same assembly class, flag, category and population tag, "
        f"{len(d)} different; graph-coordinate PASS classes in the reproduction {dict(up)} (ver3 letter row: 9 normal / 3 tumor-only / 1 neither).")
    assert set(o) == set(n) and not d
else:
    say("E9d reproduction: not finished.")
(HERE / "checks.md").write_text("# ver4 offref: checks (checks.py)\n\n" + "\n".join(out) + "\n")
