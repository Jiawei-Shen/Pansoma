"""E17c chr1 calls: finished chr1 INDEL call sets re-scored under the PoN rules of the E17c grid; a leave-one-tumor-out
selection of the INDEL rule on calls.

    python e17c_calls_chr1.py run --shard 0/2 [--jobs 2 --threads 2]    rtg runs (Slurm); DATA/calls_chr1_v2/<call set>/
    python e17c_calls_chr1.py report                                    -> calls_chr1.md, heldout.md, calls_chr1.json

Call sets: every finished chr1 INDEL call set under RUNS (../call_level_common/calls_lib.discover: vcf.sh done): in-sample
at the run's threshold ("main") and at the validation-chosen threshold ("valbest", vcf_chr1_valbest), the joint
fiberseq model per purity (vcf_chr1_<SET>) and cross-sample.

Re-scoring, per call set and rule (nothing in the pipeline changes):
- tagging: the calls' GRCh38 VCF (<set>.INDEL.linear.vcf.gz, every record) tagged as scripts/filter_panel_of_normals.py
  tags it, under the rule (PoN features of e17c_scan.py: record_matches on the PoNs streamed over chr1);
- PASS: rtg vcfeval of the PASS records the rule does not tag -> TP / FP / FN, P, R, F1 at the run's threshold;
- curve: the PASS + LowQual records it does not tag, --all-records -> best F1 (pansoma_net_v2.vcfeval.describe).
  Both with pansoma_net_v2.vcfeval.vcfeval (imported) on the run's own eval_bed/truth.INDEL.vcf.gz and regions.bed, -e the
  run's BED, the pipeline's rtg options (vcf.sh). This is the "pon" half of the pipeline's evaluation.
- Rules: PASS for all 11 INDEL grid rules (e17c_tables.INDEL_GRID) on the selection call sets (in-sample main and
  valbest, joint P100); PASS + curve for the current rule and the rule selected on chr2-22 tensors on every call set.
  Rules that leave the same records share one rtg run. Results of the first version of this script
  (DATA/calls_chr1/<call set>/{current,selected}/) are reused for those two rules.
- Check: the current rule must reproduce each run's own eval_bed report (PASS row and best F1).

Leave-one-tumor-out selection on calls (fixed here before any grid result was computed):
- Selection call sets of a tumor: in-sample, run threshold, P100, the tumor's own model.
  - HG008T: Illumina / PacBio / ONT x linear / d9 / d46.
  - COLO829T: Illumina / fiberseq / ONT x linear / d46, and the joint fiberseq model on d9 at P100.
- Rule selected on tumor A = the INDEL grid rule with the highest mean PASS F1 over A's call sets (each call set
  weighted equally). Rules within 0.0005 of the best are tied; among them, the one with the higher mean recall is taken.
- Held-out score: the mean PASS F1 of that rule on tumor B's call sets, compared with the current rule and with B's own
  best rule. Both directions: HG008T -> COLO829T and COLO829T -> HG008T.
- Per platform: the same within Illumina, PacBio-type (HG008T PacBio <-> COLO829T fiberseq) and ONT.
- Sensitivity: the same on the valbest call sets.
"""
from __future__ import annotations

import argparse
import hashlib
import json
import statistics
import sys
import time
from argparse import Namespace
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path

HERE = Path(__file__).resolve().parent
sys.path[:0] = [str(HERE), str(HERE.parent / "call_level_common")]
import e17c_scan as scan  # noqa: E402
import e17c_tables as tables  # noqa: E402
from calls_lib import discover, sort_key  # noqa: E402
from filter_panel_of_normals import PON_NAMES, record_kind  # noqa: E402
from pansoma_net_v2 import vcfeval as pipeline_vcfeval  # noqa: E402
from pansoma_net_v2.vcfio import read_vcf, write_vcf  # noqa: E402

CHROM = "chr1"
OUT = scan.DATA / "calls_chr1_v2"
OLD = scan.DATA / "calls_chr1"   # first version: evaluate() with --pon-vcf, current / selected
SDF = "/scratch/jshen/REF/GRCh38.primary.sdf"
RTG = "/scratch/jshen/bin/rtg-tools-3.12.1/rtg"
CURVE_FILTERS = pipeline_vcfeval.CURVE_FILTERS
TIE = 0.0005
PLATFORM = {"Illumina": "Illumina", "PacBio": "PacBio-type", "fiberseq": "PacBio-type", "ONT": "ONT"}


# ------------------------------------------------------------------ rules
def pon_names(rule, pon_kind, fts):
    """PoNs (script names) that tag a record of the script's kind with ALT features fts under rule."""
    r = rule if rule["kind"] == pon_kind else tables.CURRENT[pon_kind]
    hits = [False] * 4
    for ft in fts:
        if r["kind"] == "SNV":
            fi = tables.FI[r["floor"]]
            hits[0] |= ft[tables.G_LVL] >= fi
            hits[1] |= bool(r["dbsnp"] and ft[tables.D_NONSOM])
            hits[2] |= bool(r["kg"] and ft[tables.K_ANY])
            hits[3] |= ft[tables.C_LVL] >= fi
        else:
            floor = r["hp"] if r["h"] is not None and ft[tables.HP] >= r["h"] else r["other"]
            hits[0] |= ft[tables.G_LVL] >= tables.FI[floor]
            hits[3] |= ft[tables.C_LVL] >= tables.FI[floor]
    return tuple(n for n, h in zip(PON_NAMES, hits) if h)


def alt_keys(rec):
    return [(rec[0], int(rec[1]), rec[3].upper(), a.upper()) for a in rec[4].split(",")
            if a and not a.startswith("<") and a not in ("*", ".")]


def selected_rule():
    name = json.loads((HERE / "summary.json").read_text())["selected"]["INDEL|chr2-22|d9|f1p"]
    return next(r for r in tables.INDEL_GRID if r["name"] == name)


def slug(name):
    return hashlib.sha1(name.encode()).hexdigest()[:8]


# ------------------------------------------------------------------ call sets
def callsets():
    out = [cs for cs in discover(groups=("insample", "xs", "joint"), variants=("main", "valbest"))
           if cs.kind == "INDEL" and (cs.vcf_dir / "eval_bed" / "report.json").exists()]
    return sorted(out, key=sort_key)


def is_selection_set(cs):
    return cs.group == "insample" or (cs.group == "joint" and cs.purity == "P100")


def rules_for(cs, sel):
    rules = {r["name"]: r for r in (tables.INDEL_GRID if is_selection_set(cs) else [])}
    rules[tables.CURRENT["INDEL"]["name"]] = tables.CURRENT["INDEL"]
    rules[sel["name"]] = sel
    curve = {tables.CURRENT["INDEL"]["name"], sel["name"]}
    return rules, curve


# ------------------------------------------------------------------ run
def rtg_args(report, threads):
    return Namespace(rtg=RTG, rtg_mem="4g", sdf=SDF, threads=threads, bed=report["bed"])


def run_one(task):
    """One rtg vcfeval; returns (out dir, summary, roc). A finished directory is read back, an unfinished one moved."""
    args, truth, vcf, out, regions, curve = task
    if (out / "done").exists():
        return out, pipeline_vcfeval.read_summary(out), pipeline_vcfeval.read_roc(out)
    if out.exists():
        out.rename(out.with_name(out.name + f".unfinished.{int(time.time())}"))
    r = pipeline_vcfeval.vcfeval(args, truth, vcf, out, regions, curve=curve)
    return out, r["summary"], r["roc"]


def old_result(cs, sel, name):
    tag = "current" if name == tables.CURRENT["INDEL"]["name"] else "selected" if name == sel["name"] else None
    if tag is None:
        return None
    p = OLD / cs.id / tag / "eval_bed" / "report.json"
    if not p.exists():
        return None
    r = json.loads(p.read_text())["pon"]
    return dict(at_threshold=r["at_threshold"], best=r["best"], source=str(p))


def run(args):
    sel = selected_rule()
    sets = callsets()
    k, n = map(int, args.shard.split("/"))
    mine = [cs for i, cs in enumerate(sets) if i % n == k and (args.only is None or args.only in cs.id)]
    print(f"{len(sets)} INDEL call sets, shard {args.shard}: {len(mine)}", flush=True)
    records = {}
    keys = set()
    for cs in mine:
        header, columns, recs = read_vcf(cs.linear_vcf)
        recs = list(recs)
        records[cs.id] = (header, columns, recs)
        for rec in recs:
            if rec[0] != CHROM:
                raise SystemExit(f"{cs.id}: record off {CHROM}")
            keys.update(alt_keys(rec))
    fasta_path = json.loads((scan.TENSORS / scan.D9[0] / "tensors" / "SNV" / "labels.manifest.json").read_text())["fasta"]
    repeats, fasta = scan.RepeatContext(fasta_path), scan.MemFasta(fasta_path, CHROM)
    repeats.genome = fasta
    positions = {k_[1] for k_ in keys}
    pons, contigs = [None] * 4, [None] * 4
    for p in range(4):
        _, contig, recs, _ = scan.stream_pon((p, CHROM, positions))
        pons[p], contigs[p] = scan.MemPoN(scan.PON_FILES[p], recs), contig
    feat = {k_: scan.features(pons, contigs, k_, fasta, repeats) for k_ in keys}
    print(f"{len(keys):,} call alleles, features done", flush=True)

    tasks, plans = [], {}
    for cs in mine:
        header, columns, recs = records[cs.id]
        report = json.loads((cs.vcf_dir / "eval_bed" / "report.json").read_text())
        truth, regions = cs.vcf_dir / "eval_bed" / "truth.INDEL.vcf.gz", cs.vcf_dir / "eval_bed" / "regions.bed"
        rargs = rtg_args(report, args.threads)
        d = OUT / cs.id
        d.mkdir(parents=True, exist_ok=True)
        rules, curve_rules = rules_for(cs, sel)
        filt = [set(r[6].split(";")) for r in recs]
        plan = {}
        written = set()
        for name, rule in rules.items():
            tagged = {i for i, r in enumerate(recs)
                      if pon_names(rule, record_kind(r[3], r[4].split(",")), [feat[x] for x in alt_keys(r)])}
            entry = dict(rule=name, tagged_records=len(tagged))
            old = old_result(cs, sel, name)
            for part, keep in (("pass", lambda f: f <= {"PASS", "."}),
                               ("curve", lambda f: f <= CURVE_FILTERS)):
                if part == "curve" and name not in curve_rules:
                    continue
                kept = [i for i in range(len(recs)) if keep(filt[i]) and i not in tagged]
                h = hashlib.sha1((part + ":" + ",".join(map(str, kept))).encode()).hexdigest()[:12]
                entry[part] = h
                entry[f"{part}_records"] = len(kept)
                if old is not None:
                    entry["reused"] = old["source"]
                    continue
                if h in written:
                    continue
                written.add(h)
                vcf = d / f"{part}_{h}.vcf.gz"
                if not (d / f"{part}_{h}" / "done").exists():
                    write_vcf(vcf, header, columns, [recs[i] for i in kept])
                tasks.append((rargs, truth, vcf, d / f"{part}_{h}", regions, part == "curve"))
            plan[name] = entry
        plans[cs.id] = plan
        (d / "plan.json").write_text(json.dumps(dict(id=cs.id, rules=plan), indent=1) + "\n")
    print(f"{len(tasks)} rtg runs", flush=True)
    done = 0
    with ThreadPoolExecutor(args.jobs) as pool:
        for out, _, _ in pool.map(run_one, tasks):
            done += 1
            if done % 10 == 0 or done == len(tasks):
                print(f"{done}/{len(tasks)} rtg runs done ({out.parent.name})", flush=True)


# ------------------------------------------------------------------ report
def results(cs, sel):
    """{rule name: dict(at_threshold, best or None)} from the plan and the rtg directories (or the old reports)."""
    d = OUT / cs.id
    p = d / "plan.json"
    if not p.exists():
        return None
    plan = json.loads(p.read_text())["rules"]
    out = {}
    for name, e in plan.items():
        old = old_result(cs, sel, name)
        if old is not None:
            out[name] = dict(at_threshold=old["at_threshold"], best=old["best"], tagged_records=e["tagged_records"])
            continue
        pdir = d / f"pass_{e['pass']}"
        if not (pdir / "done").exists():
            return None
        row = dict(at_threshold=next(r for r in pipeline_vcfeval.read_summary(pdir) if r["threshold"] is None),
                   best=None, tagged_records=e["tagged_records"])
        if "curve" in e:
            cdir = d / f"curve_{e['curve']}"
            if not (cdir / "done").exists():
                return None
            roc = pipeline_vcfeval.read_roc(cdir)
            row["best"] = max(roc, key=lambda r: r["f1"]) if roc else None
        out[name] = row
    return out


def f1(row):
    return row["at_threshold"]["f1"] or 0.0


def mean(xs):
    return sum(xs) / len(xs) if xs else float("nan")


def select(res, ids):
    rows = []
    for rule in tables.INDEL_GRID:
        name = rule["name"]
        rows.append(dict(rule=name, f1=mean([f1(res[i][name]) for i in ids]),
                         recall=mean([res[i][name]["at_threshold"]["recall"] or 0.0 for i in ids])))
    best = max(r["f1"] for r in rows)
    tied = [r for r in rows if r["f1"] >= best - TIE]
    chosen = max(tied, key=lambda r: (r["recall"], r["f1"]))["rule"]
    rows.sort(key=lambda r: -r["f1"])
    for k, r in enumerate(rows, 1):
        r["rank"] = k
    return chosen, rows, [r["rule"] for r in tied]


def report(args):
    sel = selected_rule()
    cur = tables.CURRENT["INDEL"]["name"]
    sets = callsets()
    res, meta, missing = {}, {}, []
    for cs in sets:
        r = results(cs, sel)
        if r is None:
            missing.append(cs.id)
            continue
        pipe = json.loads((cs.vcf_dir / "eval_bed" / "report.json").read_text())
        res[cs.id] = r
        meta[cs.id] = dict(id=cs.id, label=cs.label, group=cs.group, variant=cs.variant, sample=cs.sample,
                           platform=cs.platform, arm=cs.arm, purity=cs.purity, model_sample=cs.model_sample,
                           pipeline={k: pipe["pon"][k] for k in ("at_threshold", "best")},
                           raw={k: pipe["raw"][k] for k in ("at_threshold", "best")})
    gen = time.strftime("%Y-%m-%d %H:%M")
    # ---- calls_chr1.md: current vs chr2-22 selection on every call set
    L = ["# E17c: chr1 calls under the current PoN rule and the rule selected on chr2-22 tensors", "",
         f"Generated by `e17c_calls_chr1.py report` ({gen}). Every finished chr1 INDEL call set under `{tables.RUNS}` "
         "(vcf.sh done), inside the benchmark BED. PASS calls at the run's threshold: TP / FP / FN, P, R, F1; best F1 on "
         "the PASS + LowQual curve. *pipeline* = the run's own eval_bed report; *current (re-scored)* = this script's "
         f"tagging with the current rule (must equal pipeline); *selected* = the chr2-22 tensor selection ({sel['name']}); "
         "*no PoN* = the pipeline's raw evaluation.", ""]
    if missing:
        L += [f"Not finished here (no results yet): {', '.join(missing)}.", ""]

    def cell(d):
        a, b = d["at_threshold"], d["best"] or {}
        return (f"{a['tp_baseline']} / {a['fp']} / {a['fn']} | {a['precision'] or 0:.3f} | {a['recall'] or 0:.3f} | "
                f"**{a['f1'] or 0:.3f}** | {b.get('f1', float('nan')):.3f}")

    L += [f"## INDEL: selected on chr2-22 tensors ({sel['name']}) vs current ({cur})", "",
          "| call set | rule | chr1 TP / FP / FN | chr1 P | chr1 R | chr1 F1 | chr1 best F1 (curve) |",
          "|---|---|---|---:|---:|---:|---:|"]
    same = 0
    for i in res:
        m, r = meta[i], res[i]
        ok = r[cur]["at_threshold"] == m["pipeline"]["at_threshold"] and r[cur]["best"] == m["pipeline"]["best"]
        same += ok
        L.append(f"| {m['label']} | pipeline | {cell(m['pipeline'])} |")
        L.append(f"| | current (re-scored){'' if ok else ' **DIFFERS**'} | {cell(r[cur])} |")
        L.append(f"| | selected (chr2-22 tensors) | {cell(r[sel['name']])} |")
        L.append(f"| | no PoN | {cell(m['raw'])} |")
    L.append("")
    L.append(f"- Re-scored current rule identical to the pipeline report (PASS row and best F1): {same} of {len(res)}.")
    up = sum(f1(res[i][cur]) > (meta[i]["raw"]["at_threshold"]["f1"] or 0) for i in res)
    raw_p = [meta[i]["raw"]["at_threshold"]["precision"] or 0 for i in res]
    L.append(f"- Current rule vs no PoN: chr1 PASS F1 higher in {up} of {len(res)}. chr1 precision at the threshold "
             f"without the PoN: {min(raw_p):.3f}-{max(raw_p):.3f}.")
    L.append("")
    L += ["Selected (chr2-22 tensors) minus current:", "",
          "| call sets | n | chr1 PASS F1: better / worse / equal | median | min | max | chr1 best F1: better / worse / "
          "equal | median | min | max |", "|---|---:|---|---:|---:|---:|---|---:|---:|---:|"]
    groups = [("all", lambda m: True),
              ("in-sample, run threshold", lambda m: m["group"] == "insample" and m["variant"] == "main"),
              ("in-sample, valbest threshold", lambda m: m["variant"] == "valbest"),
              ("joint fiberseq (per purity)", lambda m: m["group"] == "joint"),
              ("cross-sample", lambda m: m["group"] == "xs")]
    for name, keep in groups:
        ids = [i for i in res if keep(meta[i])]
        if not ids:
            continue
        cells = []
        for get in (lambda d: d["at_threshold"]["f1"] or 0.0, lambda d: (d["best"] or {}).get("f1", 0.0)):
            dd = [round(get(res[i][sel["name"]]) - get(res[i][cur]), 4) for i in ids]
            cells.append(f"{sum(x > 0 for x in dd)} / {sum(x < 0 for x in dd)} / {sum(x == 0 for x in dd)} | "
                         f"{statistics.median(dd):+.3f} | {min(dd):+.3f} | {max(dd):+.3f}")
        L.append(f"| {name} | {len(ids)} | " + " | ".join(cells) + " |")
    L.append("")
    # SNV: the chr2-22 selection is the current rule; the first version re-scored the in-sample d9 call sets with it
    snv = [cs for cs in discover(groups=("insample",), variants=("main",)) if cs.kind == "SNV" and cs.arm == "d9"
           and (OLD / cs.id / "current" / "eval_bed" / "report.json").exists()]
    if snv:
        L += ["## SNV: the chr2-22 selection is the current rule; reproduction check only", "",
              "Re-scored by the first version of this script (`pansoma_net_v2.vcfeval.evaluate` with `--pon-vcf`), "
              f"`{OLD}/<call set>/current/`.", "",
              "| call set | rule | chr1 TP / FP / FN | chr1 P | chr1 R | chr1 F1 | chr1 best F1 (curve) |",
              "|---|---|---|---:|---:|---:|---:|"]
        n_same, raw_snv = 0, []
        for cs in sorted(snv, key=sort_key):
            full = json.loads((cs.vcf_dir / "eval_bed" / "report.json").read_text())
            pipe, raw = full["pon"], full["raw"]
            raw_snv.append(raw["at_threshold"]["precision"] or 0)
            mine = json.loads((OLD / cs.id / "current" / "eval_bed" / "report.json").read_text())["pon"]
            ok = mine["at_threshold"] == pipe["at_threshold"] and mine["best"] == pipe["best"]
            n_same += ok
            L.append(f"| {cs.label} | pipeline | {cell(pipe)} |")
            L.append(f"| | current (re-scored){'' if ok else ' **DIFFERS**'} | {cell(mine)} |")
            L.append(f"| | no PoN | {cell(raw)} |")
        L += ["", f"- Re-scored current rule identical to the pipeline report: {n_same} of {len(snv)}.",
              f"- chr1 precision at the threshold without the PoN: {min(raw_snv):.3f}-{max(raw_snv):.3f}.", ""]
    (HERE / "calls_chr1.md").write_text("\n".join(L) + "\n")

    # ---- heldout.md: leave-one-tumor-out selection on calls
    H = ["# E17c: leave-one-tumor-out selection of the INDEL PoN rule on chr1 calls", "",
         f"Generated by `e17c_calls_chr1.py report` ({gen}). Rules: the 11 INDEL rules of the E17c grid. Score: chr1 PASS "
         "F1 after the PoN at the run's threshold (rtg, benchmark BED), mean over the call sets (each weighted equally). "
         f"Selection on one tumor's call sets, ties within {TIE} -> higher mean recall; then scored on the other "
         "tumor's call sets (definitions: the docstring of `e17c_calls_chr1.py`).", ""]
    summary = {}
    for variant, vname in (("main", "run threshold"), ("valbest", "validation-chosen threshold")):
        def tumor_ids(sample, platform=None):
            return [i for i, m in meta.items() if m["sample"] == sample and m["purity"] == "P100" and m["variant"] == variant
                    and (m["group"] == "insample" or (m["group"] == "joint" and variant == "main"))
                    and m["model_sample"] == sample and (platform is None or PLATFORM[m["platform"]] == platform)
                    and all(r["name"] in res[i] for r in tables.INDEL_GRID)]
        H += [f"## Call sets at the {vname}", ""]
        for scope in [None, "Illumina", "PacBio-type", "ONT"]:
            ids = {s: tumor_ids(s, scope) for s in tables.SAMPLES}
            if not all(ids.values()):
                H += [f"### {scope or 'all platforms'}: not enough call sets "
                      f"({', '.join(f'{s} {len(v)}' for s, v in ids.items())})", ""]
                continue
            sels = {s: select(res, ids[s]) for s in tables.SAMPLES}
            H += [f"### {scope or 'all platforms'} (HG008T {len(ids['HG008T'])} call sets, COLO829T "
                  f"{len(ids['COLO829T'])})", ""]
            H += ["| rule | HG008T chr1 mean PASS F1 | HG008T rank | COLO829T chr1 mean PASS F1 | COLO829T rank |",
                  "|---|---:|---:|---:|---:|"]
            rank = {s: {r["rule"]: r for r in sels[s][1]} for s in tables.SAMPLES}
            for r in sels["HG008T"][1]:
                nm = r["rule"]
                marks = [x for x, c in (("current", nm == cur), ("chr2-22 tensors", nm == sel["name"]),
                                        ("selected on HG008T", nm == sels["HG008T"][0]),
                                        ("selected on COLO829T", nm == sels["COLO829T"][0])) if c]
                H.append(f"| {nm}{' (' + ', '.join(marks) + ')' if marks else ''} | {rank['HG008T'][nm]['f1']:.4f} | "
                         f"{rank['HG008T'][nm]['rank']} | {rank['COLO829T'][nm]['f1']:.4f} | {rank['COLO829T'][nm]['rank']} |")
            H.append("")
            H += ["| selected on | rule | held-out tumor | held-out chr1 mean PASS F1 | current rule there | change vs current | "
                  "held-out tumor's own best | per call set vs current: better / worse / equal |",
                  "|---|---|---|---:|---:|---:|---:|---|"]
            for a, b in (("HG008T", "COLO829T"), ("COLO829T", "HG008T")):
                chosen, _, tied = sels[a]
                got, base = rank[b][chosen]["f1"], rank[b][cur]["f1"]
                own = sels[b][1][0]
                dd = [round(f1(res[i][chosen]) - f1(res[i][cur]), 4) for i in ids[b]]
                H.append(f"| {a} | {chosen}" + (f" (tied: {'; '.join(t for t in tied if t != chosen)})" if len(tied) > 1 else "")
                         + f" | {b} | {got:.4f} | {base:.4f} | {got - base:+.4f} | {own['f1']:.4f} ({own['rule']}) | "
                         f"{sum(x > 0 for x in dd)} / {sum(x < 0 for x in dd)} / {sum(x == 0 for x in dd)} |")
                summary[f"{variant}|{scope or 'all'}|{a}->{b}"] = dict(
                    selected=chosen, tied=tied, heldout_f1=got, current_f1=base, change=got - base,
                    heldout_own_best=own["rule"], heldout_own_best_f1=own["f1"], n_select=len(ids[a]), n_heldout=len(ids[b]),
                    per_set_better=sum(x > 0 for x in dd), per_set_worse=sum(x < 0 for x in dd))
            H.append("")
            H += ["Per call set, chr1 PASS F1 (current / selected on the other tumor):", "",
                  "| call set | current | " + " | ".join(f"selected on {a}" for a in tables.SAMPLES) + " |",
                  "|---|---:|" + "---:|" * len(tables.SAMPLES)]
            for s in tables.SAMPLES:
                for i in ids[s]:
                    H.append(f"| {meta[i]['label']} | {f1(res[i][cur]):.3f} | "
                             + " | ".join(f"{f1(res[i][sels[a][0]]):.3f}" for a in tables.SAMPLES) + " |")
            H.append("")
    (HERE / "heldout.md").write_text("\n".join(H) + "\n")
    (HERE / "calls_chr1.json").write_text(json.dumps(dict(generated=gen, missing=missing, meta=meta, results=res,
                                                          heldout=summary), indent=1) + "\n")
    print(f"{len(res)} call sets reported, {len(missing)} missing; current reproduced in {same}")
    for k_, v in summary.items():
        print(k_, v["selected"], f"{v['heldout_f1']:.4f} vs current {v['current_f1']:.4f}")


def main():
    ap = argparse.ArgumentParser()
    sub = ap.add_subparsers(dest="cmd", required=True)
    r = sub.add_parser("run")
    r.add_argument("--shard", default="0/1")
    r.add_argument("--jobs", type=int, default=2, help="rtg runs at once")
    r.add_argument("--threads", type=int, default=2, help="rtg -T of each run")
    r.add_argument("--only", help="call sets whose id contains this text (a test)")
    sub.add_parser("report")
    args = ap.parse_args()
    run(args) if args.cmd == "run" else report(args)


if __name__ == "__main__":
    main()
