"""E9j step 1: COLO829T somatic SNVs in the patient frame (COLO829BL donor-specific assembly), the counterpart of
pansoma_net_v2_runs/assembly_offref_20260929/s1_windows.py (HG008T, where the truth gives the HG008-N assembly event).

The SMaHT COLO829T truth (chr1-22 SNVs; the set's labels.manifest.json) is on GRCh38 only, so each truth SNV is placed
on the COLO829BL verkko 2.1 haplotypes (hapX, hapY; the user's DSA paths) by its GRCh38 window (100-bp flanks):
  grch38   writes grch38_windows.fa (the window of each truth SNV; the SNV at offset 100)
  place    reads <hap>.grch38_windows.paf (minimap2 -cx sr --cs, job s1_map.sh) and writes, per truth and haplotype,
           the DSA position of the SNV (the primary alignment: >= 90 % of the window aligned, MAPQ >= 20, the site in a
           match or mismatch block), the haplotype base (REF check) and the N (haplotype) and T (haplotype + ALT) windows,
           100-bp flanks, DSA forward strand -> placed.tsv, hap_windows.fa (for the tumor-assembly check)
  pick     reads <tumor hap>.hap_windows.paf (the T windows on the COLO829T hifiasm hap1 / hap2, job s1_map.sh MODE=tumor)
           and chooses the event haplotype: the haplotype whose T window has a strict hit in the tumor assembly (>= 90 %
           aligned, NM <= 1, no indel, the SNV base matched; 'tumor_hit'); when both haplotypes' T windows hit, the one
           with an exact (NM 0) hit ('tumor_hit_exact'); identical N windows on both haplotypes -> hapX ('identical');
           both exact or both NM 1 -> 'ambiguous' (hapX kept); neither -> 'unresolved' (hapX kept when placed, else hapY). Writes events.tsv and windows.fq (N / T of the event haplotype; giraffe input).
usage: python s1_windows.py {grch38|place|pick} DATA_DIR
"""
import csv
import gzip
import re
import sys
from pathlib import Path

import pysam

TRUTH = "/scratch/jshen/data/pansoma_v2_tensors/COLO829T_truth/COLO829T_somatic_snv_indel.vcf.gz"
GRCH38 = "/scratch/jshen/data/HapMap/GCA_000001405.15_GRCh38_no_alt_analysis_set.fasta"
Q = "/scratch/qfu/COLO829BL_DSA"
HAPS = {"hapX": f"{Q}/SMHTCOLO829BL-X-X-M45-XX-uwsc-SMASFTH1B2J6-verkko_2.1_hapX.fa",
        "hapY": f"{Q}/SMHTCOLO829BL-X-X-M45-XX-uwsc-SMASFK4V33V8-verkko_2.1_hapY.fa"}
TUMOR = ("t_hap1", "t_hap2")
W = 100
CHROMS = [f"chr{i}" for i in range(1, 23)]
COMP = str.maketrans("ACGTN", "TGCAN")
CS = re.compile(r"(:\d+|\*[a-z][a-z]|[+-][a-z]+|~[a-z]{2}\d+[a-z]{2})")


def rc(s):
    return s.translate(COMP)[::-1]


def read_paf(path):
    """qname -> [record dict], every alignment line."""
    out = {}
    with open(path) as f:
        for line in f:
            c = line.rstrip("\n").split("\t")
            tags = {t[:2]: t[5:] for t in c[12:]}
            out.setdefault(c[0], []).append(dict(qlen=int(c[1]), qs=int(c[2]), qe=int(c[3]), strand=c[4], t=c[5],
                                                 tlen=int(c[6]), ts=int(c[7]), te=int(c[8]), mapq=int(c[11]),
                                                 nm=int(tags.get("NM", 99)), cs=tags.get("cs", ""),
                                                 primary=tags.get("tp", "P") == "P"))
    return out


def site_in_target(a, q_site):
    """Target position of query offset q_site (original query orientation) and the cs op there ('=' match, '*'
    mismatch), or (None, None) when the site is clipped or inside an indel."""
    qa = q_site if a["strand"] == "+" else a["qlen"] - 1 - q_site   # offset in the aligned (target) orientation
    q0 = a["qs"] if a["strand"] == "+" else a["qlen"] - a["qe"]
    if qa < q0:
        return None, None
    q, t = q0, a["ts"]
    for op in CS.findall(a["cs"]):
        if op[0] == ":":
            n = int(op[1:])
            if q <= qa < q + n:
                return t + (qa - q), "="
            q, t = q + n, t + n
        elif op[0] == "*":
            if q == qa:
                return t, "*"
            q, t = q + 1, t + 1
        elif op[0] == "+":
            if q <= qa < q + len(op) - 1:
                return None, None
            q += len(op) - 1
        elif op[0] == "-":
            t += len(op) - 1
        else:
            return None, None
    return None, None


def strict_hit(a, q_site):
    """E9d / E9h strict rule: >= 90 % aligned, NM <= 1, no indel, the site base matched."""
    if a["qe"] - a["qs"] < 0.9 * a["qlen"] or a["nm"] > 1 or "+" in a["cs"] or "-" in a["cs"]:
        return False
    _, op = site_in_target(a, q_site)
    return op == "="


def truth_snvs():
    with gzip.open(TRUTH, "rt") as f:
        for line in f:
            if line[0] == "#":
                continue
            c = line.split("\t", 8)
            if c[0] not in CHROMS:
                continue
            for alt in c[4].split(","):
                if len(c[3]) == 1 and len(alt) == 1 and alt in "ACGT":
                    yield c[0], int(c[1]), c[3].upper(), alt.upper()


def grch38(out):
    fa = pysam.FastaFile(GRCH38)
    n = 0
    with open(out / "grch38_windows.fa", "w") as w, open(out / "truth_snvs.tsv", "w") as t:
        t.write("truth\tchrom\tpos\tref\talt\n")
        for k, (chrom, pos, ref, alt) in enumerate(truth_snvs()):
            s0 = pos - 1
            seq = fa.fetch(chrom, s0 - W, s0 + W + 1).upper()
            assert seq[W] == ref, (chrom, pos, ref, seq[W])
            tid = f"s{k:05d}"
            w.write(f">{tid}\n{seq}\n")
            t.write(f"{tid}\t{chrom}\t{pos}\t{ref}\t{alt}\n")
            n += 1
    print(f"{n} truth SNVs (chr1-22) -> grch38_windows.fa")


def place(out):
    truths = list(csv.DictReader(open(out / "truth_snvs.tsv"), delimiter="\t"))
    fas = {h: pysam.FastaFile(p) for h, p in HAPS.items()}
    paf = {h: read_paf(out / f"{h}.grch38_windows.paf") for h in HAPS}
    rows = []
    with open(out / "hap_windows.fa", "w") as wf:
        for t in truths:
            for h in HAPS:
                prim = [a for a in paf[h].get(t["truth"], []) if a["primary"]]
                row = dict(truth=t["truth"], hap=h, status="", contig="", pos0="", strand="", mapq="", nm="",
                           hap_base="", ref_ok="", flank="")
                if not prim:
                    row["status"] = "unmapped"
                else:
                    a = prim[0]
                    row.update(contig=a["t"], strand=a["strand"], mapq=a["mapq"], nm=a["nm"])
                    pos0, op = site_in_target(a, W)
                    if a["qe"] - a["qs"] < 0.9 * a["qlen"] or a["mapq"] < 20:
                        row["status"] = "weak"
                    elif pos0 is None:
                        row["status"] = "site_in_indel"
                    else:
                        lo = max(0, pos0 - W)
                        win = fas[h].fetch(a["t"], lo, pos0 + 1 + W).upper()   # one read per window
                        left, base, right = win[:pos0 - lo], win[pos0 - lo], win[pos0 - lo + 1:]
                        ref_fwd = t["ref"] if a["strand"] == "+" else rc(t["ref"])
                        alt_fwd = t["alt"] if a["strand"] == "+" else rc(t["alt"])
                        row.update(pos0=pos0, hap_base=base, ref_ok=base == ref_fwd, flank=len(left), status="placed",
                                   alt_fwd=alt_fwd)
                        wf.write(f">{t['truth']}_{h}_N\n{left + base + right}\n>{t['truth']}_{h}_T\n{left + alt_fwd + right}\n")
                rows.append(row)
    cols = ["truth", "hap", "status", "contig", "pos0", "strand", "mapq", "nm", "hap_base", "ref_ok", "flank", "alt_fwd"]
    with open(out / "placed.tsv", "w") as f:
        w = csv.DictWriter(f, fieldnames=cols, delimiter="\t", extrasaction="ignore")
        w.writeheader()
        w.writerows(rows)
    from collections import Counter
    print("placement:", dict(Counter((r["hap"], r["status"], r["ref_ok"]) for r in rows)))


def pick(out):
    from collections import Counter
    truths = {r["truth"]: r for r in csv.DictReader(open(out / "truth_snvs.tsv"), delimiter="\t")}
    placed = {}
    for r in csv.DictReader(open(out / "placed.tsv"), delimiter="\t"):
        placed.setdefault(r["truth"], {})[r["hap"]] = r
    seqs = {}
    with open(out / "hap_windows.fa") as f:
        name = None
        for line in f:
            if line[0] == ">":
                name = line[1:].strip()
            else:
                seqs[name] = line.strip()
    tpaf = {h: read_paf(out / f"{h}.hap_windows.paf") for h in TUMOR}
    events, stat = [], Counter()
    with open(out / "windows.fq", "w") as fq:
        for tid, t in truths.items():
            ok = [h for h in ("hapX", "hapY") if placed.get(tid, {}).get(h, {}).get("status") == "placed"]
            if not ok:
                stat["not_placed"] += 1
                continue
            nm = {}   # best NM of a strict hit of the haplotype's T window in the tumor assembly (None: no hit)
            for h in ok:
                q = f"{tid}_{h}_T"
                fl = int(placed[tid][h]["flank"])
                v = [a["nm"] for th in TUMOR for a in tpaf[th].get(q, []) if strict_hit(a, fl)]
                nm[h] = min(v) if v else None
            hit = {h: nm[h] is not None for h in ok}
            exact = [h for h in ok if nm[h] == 0]
            same = len(ok) == 2 and seqs[f"{tid}_hapX_N"] == seqs[f"{tid}_hapY_N"]
            if same:
                h, how = "hapX", "identical" + ("" if hit["hapX"] else "_no_tumor_hit")
            elif sum(hit.values()) == 1:
                h, how = next(x for x in ok if hit[x]), "tumor_hit"
            elif len(exact) == 1:   # both within NM <= 1 (a germline difference): the exact one
                h, how = exact[0], "tumor_hit_exact"
            elif sum(hit.values()) == 2:
                h, how = "hapX", "ambiguous"
            else:
                h, how = ok[0], "unresolved"
            stat[how] += 1
            p = placed[tid][h]
            ev = dict(event=tid, chrom=t["chrom"], grch38_pos=t["pos"], ref=t["ref"], alt=t["alt"], hap=h, how=how,
                      contig=p["contig"], pos0=p["pos0"], strand=p["strand"], ref_ok=p["ref_ok"], flank=p["flank"],
                      grch38_records=f"{t['chrom']}:{t['pos']}:{t['ref']}>{t['alt']}")
            events.append(ev)
            for tag in ("N", "T"):
                s = seqs[f"{tid}_{h}_{tag}"]
                fq.write(f"@{tid}_{tag} {p['contig']}:{p['pos0']} flank={p['flank']}\n{s}\n+\n{'I' * len(s)}\n")
    cols = list(events[0])
    with open(out / "events.tsv", "w") as f:
        w = csv.DictWriter(f, fieldnames=cols, delimiter="\t")
        w.writeheader()
        w.writerows(events)
    print("event haplotype:", dict(stat), "| events", len(events), "| REF ok", sum(e["ref_ok"] == "True" for e in events))


if __name__ == "__main__":
    {"grch38": grch38, "place": place, "pick": pick}[sys.argv[1]](Path(sys.argv[2]))
