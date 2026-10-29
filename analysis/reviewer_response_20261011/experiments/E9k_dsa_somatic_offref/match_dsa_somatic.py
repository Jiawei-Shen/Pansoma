"""E9k: COLO829T off-reference (graph-coordinate) PASS SNV calls against somatic SNVs called on the COLO829BL donor-specific
assembly (DSA; user 2026-10-29: /wanglab/jshen/scratch_jshen/data/COLO829T/somatic_snvs.vcf.gz, DeepVariant 1.6.0, contigs
haplotype1-* = verkko 2.1 hapX, haplotype2-* = hapY).

Inputs (all from E9c / E9d, nothing realigned):
  calls       E9c DATA/runs/<run>/calls.tsv: every unplaced PASS SNV call (all, or the first 3,000 for the ONT runs), the
              gold calls and the 150 / 150 / 150 controls (ctrl_tp: placed calls that match the GRCh38 truth after the PoN;
              ctrl_fp: placed false positives; ctrl_pon: placed calls removed by the PoN), each with its ALT read context
              (alt_read, alt_read_site = 0-based half-open site interval in the context).
  alignments  E9d DATA/asm/COLO829T/COLO829BL_DSA_hap{X,Y}.paf: minimap2 -cx sr --cs --secondary=yes -N 50 of every
              context to each DSA haplotype (2026-10-23).
  rescue off  pansoma_net_v2_runs/reviewer_cpu_20261013/<model>__norescue/vcf_chr1/*.linear.unplaced.vcf.gz: PASS ids =
              the unplaced PASS calls without the off-reference rescue (now the default, 91a8261).

Per call, every alignment of its alt_read context that covers >= 90% of the context, has no indel within 5 bp of the site
and at most two mismatches in all is followed base by base (cs) to the site: the DSA position and base there. The call is
  dsa_somatic  : the DSA base differs from the call's ALT and the DSA has the call's ALT as a somatic SNV there (same
                 contig, position, REF = DSA base, ALT = call ALT) in the VCF (any FILTER; PASS counted separately);
  dsa_ref_only : the DSA base differs from the ALT (the normal does not carry it) but the VCF has no such SNV there;
  dsa_alt      : the DSA base equals the call's ALT on some hit (the normal carries the allele: germline);
  unplaced     : no usable alignment.
A call with hits of several kinds is counted by the first kind in the order above.
usage: /wanglab/jshen/anaconda3/bin/python match_dsa_somatic.py
"""
import csv
import gzip
import hashlib
import json
import re
import sys
from collections import Counter, defaultdict
from pathlib import Path

HERE = Path(__file__).resolve().parent
E9C = Path("/scratch/jshen/data/reviewer_experiments_20261013/E9c_offref_filter_validation")
RUNS = Path("/scratch/jshen/data/pansoma_net_v2_runs/reviewer_cpu_20261013")
VCF = Path("/wanglab/jshen/scratch_jshen/data/COLO829T/somatic_snvs.vcf.gz")
RUN_LIST = [("v1.1 (d9)", "Illumina", "rcpu__COLO829_Illumina_SNV_d9_normbq_sweep_e06__vcf_chr1", "COLO829_Illumina_SNV_d9_normbq_sweep_e06"),
            ("v1.1 (d9)", "Fiber-seq", "rcpu__COLO829_fiberseq_SNV_7fc7396_ddp4_b1024__vcf_chr1", "COLO829_fiberseq_SNV_7fc7396_ddp4_b1024"),
            ("v1.1 (d9)", "ONT", "rcpu__COLO829_ONT_SNV_7fc7396_ddp4_b1024__vcf_chr1", "COLO829_ONT_SNV_7fc7396_ddp4_b1024"),
            ("v2.1 (d46)", "Illumina", "rcpu__COLO829_Illumina_SNV_d46_7fc7396__vcf_chr1", "COLO829_Illumina_SNV_d46_7fc7396"),
            ("v2.1 (d46)", "Fiber-seq", "rcpu__COLO829_fiberseq_SNV_d46_7fc7396__vcf_chr1", "COLO829_fiberseq_SNV_d46_7fc7396"),
            ("v2.1 (d46)", "ONT", "rcpu__COLO829_ONT_SNV_d46_7fc7396__vcf_chr1", "COLO829_ONT_SNV_d46_7fc7396")]
GROUPS = ("unplaced_pass", "gold", "ctrl_tp", "ctrl_fp", "ctrl_pon")
# part 3: every rescue-off unplaced PASS call (norescue_contexts.py: E9c's run on the __norescue re-calls, contexts aligned
# to the DSA in E9k's own data folder); group "norescue_pass"
E9K = Path("/scratch/jshen/data/reviewer_experiments_20261013/E9k_dsa_somatic_offref")
KINDS = ("dsa_somatic", "dsa_ref_only", "dsa_alt", "unplaced")
CS = re.compile(r"(:\d+|\*[a-z][a-z]|[+-][a-z]+)")
COMP = str.maketrans("ACGTN", "TGCAN")


def qid(seq, site):  # as E9d e9d_asm.qid
    return "q" + hashlib.sha1(f"{seq}|{site}".encode()).hexdigest()[:16]


def load_vcf():
    out = {}
    with gzip.open(VCF, "rt") as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            f = line.split("\t", 8)
            for alt in f[4].split(","):
                out[(f[0], int(f[1]), f[3].upper(), alt.upper())] = f[6]
    return out


def site_on_target(paf, site0):
    """(contig, 1-based position, target base, query base in target orientation, near-site indel) of query position site0."""
    f = paf.rstrip("\n").split("\t")
    qlen, qs, qe, strand, tname, ts = int(f[1]), int(f[2]), int(f[3]), f[4], f[5], int(f[7])
    cs = next((x[5:] for x in f[12:] if x.startswith("cs:Z:")), None)
    if cs is None:
        return None, []
    # walk target-forward; query position runs forward ('+') or backward ('-')
    t, q, step = ts, (qs if strand == "+" else qe - 1), (1 if strand == "+" else -1)
    events = []  # (target pos, query pos, kind, tbase, qbase)
    for tok in CS.findall(cs):
        if tok[0] == ":":
            n = int(tok[1:])
            if (q - site0) * step <= 0 <= (q + step * (n - 1) - site0) * step:
                k = (site0 - q) * step
                events.append(("match", (tname, t + k + 1, None, None)))
            t += n
            q += step * n
        elif tok[0] == "*":
            if q == site0:
                events.append(("sub", (tname, t + 1, tok[1].upper(), tok[2].upper())))
            t += 1
            q += step
        elif tok[0] == "+":  # insertion in the query
            n = len(tok) - 1
            if abs(q - site0) <= 5 or abs(q + step * n - site0) <= 5:
                events.append(("indel", None))
            q += step * n
        else:  # deletion from the query
            n = len(tok) - 1
            if abs(q - site0) <= 5:
                events.append(("indel", None))
            t += n
    return strand, events


def classify(strand, ev, alt):
    """[(kind, contig, pos, DSA base, ALT)] in DSA-forward orientation (the call's ALT complemented on '-' hits)."""
    alt_t = alt if strand == "+" else alt.translate(COMP)
    if any(e[0] == "indel" for e in ev):
        return []
    site = [e for e in ev if e[0] in ("match", "sub")]
    if not site:
        return []
    kind, (contig, pos, tbase, qbase) = site[0]
    if kind == "match":  # the DSA has the query base (the ALT) at the site
        return [("dsa_alt", contig, pos, alt_t, alt_t)]
    return [("ref", contig, pos, tbase, alt_t)]


def main():
    sys.path.insert(0, str(HERE.parents[0] / "E9c_offref_filter_validation"))
    csv.field_size_limit(sys.maxsize)
    vcf = load_vcf()
    print(f"{len(vcf):,} DSA somatic SNV alleles", flush=True)
    calls, wanted = [], set()
    for graph, plat, run, model in RUN_LIST:
        norescue = set()
        nr = list((RUNS / f"{model}__norescue" / "vcf_chr1").glob("*.linear.unplaced.vcf.gz"))
        if nr:
            with gzip.open(nr[0], "rt") as fh:
                norescue = {l.split("\t")[2] for l in fh if not l.startswith("#") and l.split("\t")[6] == "PASS"}
        rj = json.loads((E9C / "runs" / run / "run.json").read_text())
        for r in csv.DictReader(open(E9C / "runs" / run / "calls.tsv"), delimiter="\t"):
            if r["group"] not in GROUPS or r.get("kind_event") != "SNP" or not r.get("alt_read"):
                continue
            s0 = int(r["alt_read_site"].split("-")[0])
            q = qid(r["alt_read"], r["alt_read_site"])
            calls.append(dict(graph=graph, platform=plat, run=run, group=r["group"], candidate_id=r["candidate_id"],
                              alt=r["evalt"].upper(), qid=q, site0=s0, rescue_off_pass=r["candidate_id"] in norescue,
                              unplaced_pass_total=rj.get("unplaced_pass"), norescue_pass_total=len(norescue)))
            wanted.add(q)
        nrun = E9K / "runs" / f"norescue__{model}" / "calls.tsv"
        if nrun.exists():
            for r in csv.DictReader(open(nrun), delimiter="\t"):
                if r["group"] != "unplaced_pass" or r.get("kind_event") != "SNP" or not r.get("alt_read"):
                    continue
                q = qid(r["alt_read"], r["alt_read_site"])
                calls.append(dict(graph=graph, platform=plat, run=run, group="norescue_pass", candidate_id=r["candidate_id"],
                                  alt=r["evalt"].upper(), qid=q, site0=int(r["alt_read_site"].split("-")[0]),
                                  rescue_off_pass=True, unplaced_pass_total=rj.get("unplaced_pass"),
                                  norescue_pass_total=len(norescue)))
                wanted.add(q)
    hits = defaultdict(list)
    pafs = [E9C / "asm" / "COLO829T" / f"COLO829BL_DSA_hap{h}.paf" for h in ("X", "Y")]
    pafs += [p for p in (E9K / "asm" / f"COLO829BL_DSA_hap{h}.paf" for h in ("X", "Y")) if p.exists()]
    for paf in pafs:
        with open(paf) as fh:
            for line in fh:
                f = line.split("\t", 12)
                if f[0] not in wanted:
                    continue
                qlen, qs, qe, nm = int(f[1]), int(f[2]), int(f[3]), None
                if (qe - qs) < 0.9 * qlen:
                    continue
                nm = next((int(x[5:]) for x in line.split("\t")[12:] if x.startswith("NM:i:")), 99)
                if nm > 3:  # site substitution + at most two other differences
                    continue
                hits[f[0]].append(line)
    rows = []
    for c in calls:
        kinds = []
        for line in hits.get(c["qid"], []):
            strand, ev = site_on_target(line, c["site0"])
            if ev:
                kinds += classify(strand, ev, c["alt"])
        best, where = "unplaced", ""
        dsa_ref = [k for k in kinds if k[0] == "ref" and k[3] != k[4]]
        som = [k for k in dsa_ref if (k[1], k[2], k[3], k[4]) in vcf]
        if som:
            best, k = "dsa_somatic", som[0]
            where = f"{k[1]}:{k[2]}:{k[3]}>{k[4]}:{vcf[(k[1], k[2], k[3], k[4])]}"
        elif dsa_ref:
            best, k = "dsa_ref_only", dsa_ref[0]
            where = f"{k[1]}:{k[2]}:{k[3]}>{k[4]}"
        elif any(k[0] == "dsa_alt" or (k[0] == "ref" and k[3] == k[4]) for k in kinds):
            best = "dsa_alt"
        c.update(kind=best, dsa=where, n_hits=len(kinds), multi=len({(k[1], k[2]) for k in kinds}) > 2)
        rows.append(c)
    cols = ["graph", "platform", "run", "group", "candidate_id", "alt", "rescue_off_pass", "kind", "dsa", "n_hits", "multi"]
    with gzip.open(HERE / "per_call.tsv.gz", "wt") as w:
        w.write("\t".join(cols) + "\n")
        for r in rows:
            w.write("\t".join(str(r[c]) for c in cols) + "\n")
    L = ["# E9k: COLO829T graph-coordinate PASS SNV calls against the COLO829BL-DSA somatic SNVs (match_dsa_somatic.py)", "",
         f"DSA somatic VCF: {VCF} ({len(vcf):,} SNV alleles). Calls: E9c contexts aligned to DSA hapX / hapY (E9d PAF).", "",
         "| graph | platform | group | calls (of unplaced PASS) | DSA somatic (PASS) | normal REF, not in VCF | normal carries ALT | not placed |",
         "|---|---|---|---:|---|---:|---:|---:|"]
    summary = []
    for graph, plat, run, model in RUN_LIST:
        for group in GROUPS + ("unplaced_pass, rescue off", "norescue_pass"):
            if group.startswith("unplaced_pass,"):
                sel = [r for r in rows if r["run"] == run and r["group"] == "unplaced_pass" and r["rescue_off_pass"]]
            else:
                sel = [r for r in rows if r["run"] == run and r["group"] == group]
            if not sel:
                continue
            n = Counter(r["kind"] for r in sel)
            som_pass = sum(1 for r in sel if r["kind"] == "dsa_somatic" and r["dsa"].endswith(":PASS"))
            tot = ""
            if group == "unplaced_pass":
                tot = f" (of {sel[0]['unplaced_pass_total']:,})"
            elif group.startswith("unplaced_pass,") or group == "norescue_pass":
                tot = f" (of {sel[0]['norescue_pass_total']:,})"
            summary.append(dict(graph=graph, platform=plat, group=group, n=len(sel), **{k: n.get(k, 0) for k in KINDS},
                                dsa_somatic_pass=som_pass))
            L.append(f"| {graph} | {plat} | {group} | {len(sel):,}{tot} | {n.get('dsa_somatic', 0)} ({som_pass}) | "
                     f"{n.get('dsa_ref_only', 0)} | {n.get('dsa_alt', 0)} | {n.get('unplaced', 0)} |")
    som = [r for r in rows if r["group"] in ("unplaced_pass", "gold", "norescue_pass") and r["kind"] == "dsa_somatic"]
    L += ["", f"## Graph-coordinate calls that match a DSA somatic SNV ({len(som)} call records)", "",
          "| graph | platform | group | candidate id | rescue off PASS | DSA site (contig:pos:REF>ALT:FILTER) |", "|---|---|---|---|---|---|"]
    for r in som:
        L.append(f"| {r['graph']} | {r['platform']} | {r['group']} | `{r['candidate_id']}` | {r['rescue_off_pass']} | {r['dsa']} |")
    (HERE / "tables.md").write_text("\n".join(L) + "\n")
    (HERE / "summary.json").write_text(json.dumps(summary, indent=1) + "\n")
    print("\n".join(L[:60]))


if __name__ == "__main__":
    main()
