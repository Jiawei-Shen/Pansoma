"""ver4 offref (1): Response Table 3.10, validation of the graph-coordinate SNV PASS calls of the nine six-tumor call
sets (snv_newpipe_20261109/multicancer6_<P>_SNV_newpipe, seed 0; HPRC v2.1, chr1, default output).

One builder for the nine rows, with the definitions of the two ver3 builders (tmp/e9i_hg008t_20261104/r33_results.py for
HG008T / COLO829T, tmp/r33_hcc1395_20261104/r33_results.py for HCC1395T; no rescue-on mode, no v1.0-only column):
  assembly (E9d, code/collect.py `asm` of EXP/per_call.tsv.gz; ALT context = read consensus when >= 40 bp, else the
  graph path; minimap2 -x sr hit covering >= 90 % of the context, NM <= 1, no indel, exact at the site):
      in a normal assembly  asm == germline    (ALT context in a matched-normal assembly)
      tumor-only            asm == tumor_only  (in the tumor assembly, in no normal assembly)
      neither               everything else    (unsupported: only the REF context is found; unverifiable: neither
                                                context is found; not checked)
  reads (E9i, e9i/normal_reads.py: reads holding the call's ALT / REF 21-mer among the primary reads over the call's
  GRCh38 locus +- 500 bp, matched normal and tumor):
      ALT in >= 2 normal reads                    status == 'ALT in normal'
      no normal ALT read, ALT in tumor reads      status == 'absent from normal' (0 normal ALT reads, >= 5 normal REF
                                                  reads) and tumor ALT reads >= 2
  tumor-only and no normal ALT read              both of the above
  also PASS on both other platforms              the call's key (node:offset:SNP:REF>ALT) is PASS in the graph-coordinate
                                                  VCF of all three platforms of the tumor
Percentages are of the call set's graph-coordinate PASS calls.

usage: python table_3_10.py           -> table_3.10.md, table_3.10.json
       python table_3_10.py --repro   -> the same code on the ver3 inputs; every cell compared with the ver3 letter's
                                         Response Table 3.10 (prints; writes nothing)
"""
import collections
import csv
import glob
import gzip
import json
import re
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
T = Path("/scratch/jshen/Github/Pansoma/tmp")
X4 = Path("/scratch/jshen/data/reviewer_experiments_20261013/ver4_20261123/offref")
LETTER = HERE.parent / "ver3_inputs/letter.md"      # pandoc Markdown of Reviewer_Letter_Response_ver3.docx
PLAT = {"Illumina": "Illumina", "PacBio_HiFi": "PacBio HiFi", "Fiber-seq": "Fiber-seq", "ONT": "ONT"}
EXPECTED_U = {"HG008T": (77, 143, 90), "COLO829T": (110, 79, 111), "HCC1395T": (76, 106, 143)}
ASM = lambda a: "normal" if a == "germline" else "tumor" if a == "tumor_only" else "neither"
csv.field_size_limit(sys.maxsize)


def rclass(r):
    s = r["status"]
    if s == "ALT in normal":
        return "alt_normal"
    if s == "absent from normal" and int(r.get("tumor_alt") or 0) >= 2:
        return "somatic"
    return "other"


def pass_keys(path):
    out = set()
    for line in gzip.open(path, "rt"):
        if line[0] == "#":
            continue
        c = line.split("\t", 8)
        if c[6] == "PASS":
            out.add(c[2].split(";")[0])
    return out


def per_call(path):
    pc, raw = collections.defaultdict(dict), collections.defaultdict(dict)
    for r in csv.DictReader(gzip.open(path, "rt"), delimiter="\t"):
        pc[r["run"]][(r["group"], r["candidate_id"])] = ASM(r["asm"])
        raw[r["run"]][(r["group"], r["candidate_id"])] = r["asm"]
    return pc, raw


def reads(dirs):
    rd, st = collections.defaultdict(dict), collections.defaultdict(dict)
    for d in dirs:
        for f in sorted(glob.glob(str(Path(d) / "*.tsv"))):
            run = Path(f).name.split(".")[0]
            for r in csv.DictReader(open(f), delimiter="\t"):
                k = (r["group"], r["candidate_id"])
                assert k not in rd[run], (f, k)
                rd[run][k] = rclass(r)
                st[run][k] = r["status"]
    return rd, st


def sources(repro):
    """[(tumor, platform, run name, run.json, per_call, [reads dirs])] in the table's row order."""
    if repro:
        a, h = T / "assessment_20261030/e9c_norescue", T / "r33_hcc1395_20261104"
        old = [("HG008T", "Illumina", "HG008_Illumina_SNV_d46_7fc7396"), ("HG008T", "PacBio HiFi", "HG008_PacBio_SNV_d46_7fc7396"),
               ("HG008T", "ONT", "HG008_ONT_SNV_d46_7fc7396"), ("COLO829T", "Illumina", "COLO829_Illumina_SNV_d46_7fc7396"),
               ("COLO829T", "Fiber-seq", "COLO829_fiberseq_SNV_d46_7fc7396"), ("COLO829T", "ONT", "COLO829_ONT_SNV_d46_7fc7396")]
        out = [(t, p, f"rcpu__{r}__norescue__vcf_chr1", a / "DATA/runs" / f"rcpu__{r}__norescue__vcf_chr1" / "run.json",
                a / "EXP/per_call.tsv.gz", [a / "E9i_out", T / "e9i_hg008t_20261104/out"]) for t, p, r in old]
        out += [("HCC1395T", p, f"mc6__HCC1395T_{tag}__norescue__vcf_chr1",
                 h / "DATA/runs" / f"mc6__HCC1395T_{tag}__norescue__vcf_chr1" / "run.json", h / "EXP/per_call.tsv.gz", [h / "e9i/out"])
                for p, tag in (("Illumina", "Illumina"), ("PacBio HiFi", "PacBio"), ("ONT", "ONT"))]
        return out
    out = []
    for line in open(HERE / "runs.txt"):
        name, vdir, tumor, plat = line.split()
        rj = X4 / "DATA/runs" / name / "run.json"
        assert Path(json.load(open(rj))["vcfdir"]) == Path(vdir)
        out.append((tumor, PLAT[plat], name, rj, HERE / "EXP/per_call.tsv.gz", [X4 / "e9i/out"]))
    return out


def letter_table():
    """Rows of the ver3 letter's Response Table 3.10 (HTML table in the pandoc Markdown)."""
    s = LETTER.read_text()
    i = s.index("Response Table 3.10\n")
    body = s[s.index("<tbody>", i):s.index("</tbody>", i)]
    rows = []
    for tr in re.findall(r"<tr[^>]*>(.*?)</tr>", body, flags=re.S):
        rows.append([re.sub(r"<[^>]+>", "", td).strip() for td in re.findall(r"<td[^>]*>(.*?)</td>", tr, flags=re.S)])
    return rows


def build(repro):
    src = sources(repro)
    pcs, rds = {}, {}
    rows, keys, ctrl = [], collections.defaultdict(dict), []
    for tumor, plat, rn, rj_path, pc_path, rdirs in src:
        if pc_path not in pcs:
            pcs[pc_path] = per_call(pc_path)
        rk = tuple(map(str, rdirs))
        if rk not in rds:
            rds[rk] = reads(rdirs)
        (pc, raw), (rd, st) = pcs[pc_path], rds[rk]
        rj = json.load(open(rj_path))
        assert rj["graph"] == "d46" and rj["kind"] == "SNV" and rj["sample"] == tumor, (rn, rj["graph"], rj["sample"])
        vcf_keys = pass_keys(rj["unplaced_vcf"])
        calls = {k: v for k, v in pc[rn].items() if k[0] == "unplaced_pass"}
        # every graph-coordinate PASS call of the VCF is evaluated, by the assemblies and by the reads
        assert {k[1] for k in calls} == vcf_keys and len(calls) == rj["unplaced_pass"] == rj["report"]["eval_nobed"]["unplaced_pass"], \
            (rn, len(calls), len(vcf_keys), rj["unplaced_pass"])
        missing = [k for k in calls if k not in rd[rn]]
        assert not missing, (rn, "calls without a read count", len(missing))
        n = len(calls)
        a = collections.Counter(calls.values())
        araw = collections.Counter(raw[rn][k] for k in calls)
        rr = {k: rd[rn][k] for k in calls}
        b = collections.Counter(rr.values())
        both = sum(1 for k, v in calls.items() if v == "tumor" and rr[k] == "somatic")
        ctp = {k: v for k, v in pc[rn].items() if k[0] == "ctrl_tp"}
        cfp = {k: v for k, v in pc[rn].items() if k[0] == "ctrl_fp"}
        rtp = [rd[rn][k] for k in ctp if k in rd[rn]]
        rfp = [rd[rn][k] for k in cfp if k in rd[rn]]
        ctrl.append(dict(tumor=tumor, platform=plat, ctrl_tp=len(ctp), ctrl_tp_tumor_only=sum(v == "tumor" for v in ctp.values()),
                         ctrl_tp_normal=sum(v == "normal" for v in ctp.values()),
                         ctrl_tp_reads=len(rtp), ctrl_tp_somatic_reads=sum(x == "somatic" for x in rtp),
                         ctrl_tp_alt_normal_reads=sum(x == "alt_normal" for x in rtp),
                         ctrl_fp=len(cfp), ctrl_fp_normal=sum(v == "normal" for v in cfp.values()),
                         ctrl_fp_tumor_only=sum(v == "tumor" for v in cfp.values()),
                         ctrl_fp_reads=len(rfp), ctrl_fp_alt_normal_reads=sum(x == "alt_normal" for x in rfp),
                         ctrl_fp_somatic_reads=sum(x == "somatic" for x in rfp)))
        rows.append(dict(tumor=tumor, platform=plat, run=rn, calls=n, normal=a["normal"], tumor_only=a["tumor"], neither=a["neither"],
                         neither_detail={k: v for k, v in araw.items() if ASM(k) == "neither"},
                         alt_normal_reads=b["alt_normal"], somatic_reads=b["somatic"], other_reads=b["other"],
                         other_reads_status=dict(collections.Counter(st[rn][k] + ("; tumor ALT < 2" if st[rn][k] == "absent from normal" else "")
                                                                     for k in calls if rr[k] == "other")),
                         both=both))
        keys[tumor][plat] = {k[1]: (v, rr[k]) for k, v in calls.items()}
    cross = {}
    for tumor, d in keys.items():
        assert len(d) == 3, (tumor, list(d))
        common = set.intersection(*(set(v) for v in d.values()))
        det = []
        for k in sorted(common):
            asm = {d[p][k][0] for p in d}
            assert len(asm) == 1, (tumor, k, asm)   # same key and same assemblies: one class
            det.append(dict(key=k, assembly=asm.pop(), reads={p: d[p][k][1] for p in d}))
        cross[tumor] = dict(n=len(common), keys=det,
                            in_two_or_more=len({k for p in d for k in d[p] if sum(k in d[q] for q in d) >= 2}),
                            union=len(set().union(*(set(v) for v in d.values()))))
        for r in rows:
            if r["tumor"] == tumor:
                r["all_three_platforms"] = len(common)
    return rows, cross, ctrl


def cells(r):
    n = r["calls"]
    f = lambda a: f"{a} ({100 * a / n:.1f}%)"
    return [f"{r['tumor']} {r['platform']}", str(n), f(r["normal"]), f(r["tumor_only"]), f(r["neither"]), f(r["alt_normal_reads"]),
            f(r["somatic_reads"]), str(r["both"]), str(r["all_three_platforms"])]


def main():
    repro = "--repro" in sys.argv
    rows, cross, ctrl = build(repro)
    if repro:
        want = letter_table()
        got = [cells(r) for r in rows]
        ok = got == want
        for g, w in zip(got, want):
            print("OK       " if g == w else "DIFFERENT", " | ".join(g), "" if g == w else f"   letter: {' | '.join(w)}")
        print(f"reproduction of the ver3 letter's Response Table 3.10 ({len(want)} rows, {sum(len(w) for w in want)} cells):",
              "OK" if ok else "FAILED")
        return 0 if ok else 1
    for r in rows:
        exp = EXPECTED_U[r["tumor"]][[x["platform"] for x in rows if x["tumor"] == r["tumor"]].index(r["platform"])]
        assert r["calls"] == exp, (r["tumor"], r["platform"], r["calls"], exp)
    head = ["Dataset (chr1 SNV)", "Graph-coordinate PASS calls", "In a normal assembly", "Tumor-only", "Neither",
            "ALT in ≥2 normal reads", "No normal ALT read, ALT in tumor reads", "Tumor-only and no normal ALT read",
            "Also PASS on both other platforms"]
    md = ["| " + " | ".join(f"**{h}**" for h in head) + " |", "|" + "---|" * len(head)]
    md += ["| " + " | ".join(cells(r)) + " |" for r in rows]
    table = "\n".join(md) + "\n"
    tot = {}
    for label, sel in (("all", None), ("HG008T", "HG008T"), ("COLO829T", "COLO829T"), ("HCC1395T", "HCC1395T")):
        s = [r for r in rows if sel is None or r["tumor"] == sel]
        tot[label] = {k: sum(r[k] for r in s) for k in ("calls", "normal", "tumor_only", "neither", "alt_normal_reads", "somatic_reads",
                                                         "other_reads", "both")}
    p = lambda a, n: f"{a} ({100 * a / n:.1f}%)"
    a = tot["all"]
    notes = ["", f"Totals, nine call sets: {a['calls']} calls; in a normal assembly {p(a['normal'], a['calls'])}, tumor-only "
             f"{p(a['tumor_only'], a['calls'])}, neither {p(a['neither'], a['calls'])}; ALT in >= 2 normal reads "
             f"{p(a['alt_normal_reads'], a['calls'])}; no normal ALT read and ALT in tumor reads {p(a['somatic_reads'], a['calls'])}; "
             f"neither read class {p(a['other_reads'], a['calls'])}; tumor-only assembly and no normal ALT read {p(a['both'], a['calls'])}.", "",
             "Per tumor (three platforms summed):", "",
             "| Tumor | Calls | In a normal assembly | Tumor-only | Neither | ALT in ≥2 normal reads | No normal ALT read, ALT in tumor reads | Tumor-only and no normal ALT read | Keys PASS on all three platforms |",
             "|---|---|---|---|---|---|---|---|---|"]
    for t in ("HG008T", "COLO829T", "HCC1395T"):
        x = tot[t]
        notes.append(f"| {t} | {x['calls']} | {p(x['normal'], x['calls'])} | {p(x['tumor_only'], x['calls'])} | {p(x['neither'], x['calls'])} | "
                     f"{p(x['alt_normal_reads'], x['calls'])} | {p(x['somatic_reads'], x['calls'])} | {x['both']} | {cross[t]['n']} |")
    notes += ["", "Keys PASS on all three platforms of a tumor (assembly class; read class per platform):", "",
              "| Tumor | Key (node:offset:SNP:REF>ALT) | Assembly | Reads Illumina / PacBio HiFi or Fiber-seq / ONT |", "|---|---|---|---|"]
    for t in ("HG008T", "COLO829T", "HCC1395T"):
        for k in cross[t]["keys"]:
            notes.append(f"| {t} | {k['key']} | {k['assembly']} | {' / '.join(k['reads'].values())} |")
        notes.append(f"| {t} | ({cross[t]['n']} keys on all three; {cross[t]['in_two_or_more']} on two or more; {cross[t]['union']} distinct keys) | | |")
    notes += ["", "'Neither' by the E9d class (unsupported = only the REF context is found in an assembly; unverifiable = neither context):",
              "; ".join(f"{r['tumor']} {r['platform']} {r['neither_detail'] or '{}'}" for r in rows) + ".", "",
              "Read class 'other' (neither read column) by E9i status:",
              "; ".join(f"{r['tumor']} {r['platform']} {r['other_reads_status'] or '{}'}" for r in rows) + ".", "",
              "Controls (placed PASS calls of the same call sets, up to 150 each; ctrl TP = rtg true positives, ctrl FP = rtg false "
              "positives after the population filter, whole chr1): expected tumor-only / no normal ALT read for TP, normal for FP.", "",
              "| Dataset | ctrl TP: tumor-only assembly | ctrl TP: in a normal assembly | ctrl TP: no normal ALT read, ALT in tumor reads | ctrl TP: ALT in ≥2 normal reads | ctrl FP: in a normal assembly | ctrl FP: tumor-only assembly | ctrl FP: ALT in ≥2 normal reads | ctrl FP: no normal ALT read, ALT in tumor reads |",
              "|---|---|---|---|---|---|---|---|---|"]
    for c in ctrl:
        notes.append(f"| {c['tumor']} {c['platform']} | {c['ctrl_tp_tumor_only']}/{c['ctrl_tp']} | {c['ctrl_tp_normal']}/{c['ctrl_tp']} | "
                     f"{c['ctrl_tp_somatic_reads']}/{c['ctrl_tp_reads']} | {c['ctrl_tp_alt_normal_reads']}/{c['ctrl_tp_reads']} | "
                     f"{c['ctrl_fp_normal']}/{c['ctrl_fp']} | {c['ctrl_fp_tumor_only']}/{c['ctrl_fp']} | "
                     f"{c['ctrl_fp_alt_normal_reads']}/{c['ctrl_fp_reads']} | {c['ctrl_fp_somatic_reads']}/{c['ctrl_fp_reads']} |")
    print(table + "\n".join(notes))
    (HERE / "table_3.10.md").write_text(
        "# Response Table 3.10 (ver4): validation of graph-coordinate SNV PASS calls\n\n"
        "Six-tumor current SNV models (snv_newpipe_20261109/multicancer6_<P>_SNV_newpipe, seed 0), HPRC v2.1, chr1, default "
        "output. Every call was evaluated. Definitions in table_3_10.py.\n\n" + table + "\n".join(notes) + "\n")
    (HERE / "table_3.10.json").write_text(json.dumps(dict(rows=rows, totals=tot, cross_platform=cross, controls=ctrl), indent=1) + "\n")
    return 0


if __name__ == "__main__":
    sys.exit(main())
