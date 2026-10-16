"""E17c tables: the PoN rule grid measured on chr2-22 (where the rule is selected) and on chr1 (the check).

    python e17c_tables.py          reads DATA/scan/chr{1..22}.json -> tables.md, summary.json (here)

OBJECTIVE (fixed in this file before any number of the grid was computed; select() is called on the chr2-22 rows
for the selection and, separately, on the chr1 rows for the check):
  For one tensor set and kind, and one rule,
    s  = fraction of the somatic-labelled tensors (label 1) the rule keeps (not tagged; an unplaced tensor is kept),
    n  = fraction of the negative tensors (label 0 non + label 2 germline, pooled) the rule keeps,
    F1 proxy = 2 s / (1 + s + n).
  This is the F1 after the PoN of a caller that, before the PoN, finds every somatic candidate and makes as many
  false calls as true ones (precision 0.5, F1 0.667), its true calls tagged like the somatic tensors and its false
  calls like the negative tensors. No model output is used.
  Selected rule of a kind = the grid rule with the highest mean F1 proxy over the six HPRC v1.1 d9 sets (each set
  weighted equally). Rules within 1e-4 of the best count as tied; among them, the one losing the fewest somatic
  tensors (mean over the sets) is taken. Per-set winners are reported beside it.
  Reported, not used for the selection: (a) the F1 proxy with germline tensors only as the negatives; (b) the
  candidate-level F1 with the real counts, 2TP / (2TP + FP + FN) with TP = somatic tensors kept, FN = somatic
  tensors tagged, FP = negative tensors kept (a caller that calls every candidate); (c) the same selection on the
  six HPRC v2.1 d46 sets.

RULE GRID (gnomAD / CoLoRSdb: exact allele with that ALT's AF >= floor; dbSNP: exact allele, not SAO=2; 1000G: exact
allele; all through filter_panel_of_normals.record_matches, see e17c_scan.py):
  SNV:   floor 1e-5, 1e-4, 1e-3, 0.01, each with and without dbSNP and with and without 1000G (16 rules); the current
         rule is 1e-4 + dbSNP + 1000G.
  INDEL: gnomAD / CoLoRSdb only; one floor for every INDEL (1e-3, 0.01), or a homopolymer split: floor `hp` when the
         ALT changes a homopolymer of >= H bases (the script's homopolymer_run), else floor `other`, with H 6, 7, 8 and
         (other, hp) = (1e-3, 0.01), (1e-3, 0.05), (1e-4, 0.01) (11 rules); the current rule is H 7, (1e-3, 0.01).
  While one kind's rule varies, records of the other kind keep the current rule (a record's kind is the script's
  record_kind). "no PoN" is listed for reference and is not a candidate.
Truth: somatic truth records of the kind (vcfeval's record_kinds) in the somatic benchmark BED, tagged as the
script tags a record (any ALT).
"""
from __future__ import annotations

import json
import time
from collections import Counter, defaultdict
from pathlib import Path

HERE = Path(__file__).resolve().parent
DATA = Path("/scratch/jshen/data/reviewer_experiments_20261013/E17c_pon_rule_offchr1")
RUNS = Path("/scratch/jshen/data/pansoma_net_v2_runs/reviewer_cpu_20261013")
CHROMS = [f"chr{i}" for i in range(1, 23)]
TRAIN = "chr2-22"
TEST = "chr1"
D9 = ("HG008T_PacBio", "HG008T_ONT", "HG008T_Illumina", "COLO829T_Illumina", "COLO829T_fiberseq", "COLO829T_ONT")
D46 = tuple(f"{s}_d46" for s in D9)
SHORT = {"HG008T_PacBio": "HG008T PacBio", "HG008T_ONT": "HG008T ONT", "HG008T_Illumina": "HG008T Illumina",
         "COLO829T_Illumina": "COLO Illumina", "COLO829T_fiberseq": "COLO fiberseq", "COLO829T_ONT": "COLO ONT"}
SAMPLES = ("HG008T", "COLO829T")
FLOORS = (1e-5, 1e-4, 1e-3, 1e-2, 5e-2)
FI = {f: i for i, f in enumerate(FLOORS)}
FTXT = {1e-5: "1e-5", 1e-4: "1e-4", 1e-3: "1e-3", 1e-2: "0.01", 5e-2: "0.05"}
G_LVL, C_LVL, D_NONSOM, K_ANY, HP, RC = 1, 3, 5, 6, 7, 8
TIE = 1e-4


# ------------------------------------------------------------------ rules
def snv_rule(floor, dbsnp, kg):
    name = f"AF>={FTXT[floor]}" + (" +dbSNP" if dbsnp else "") + (" +1000G" if kg else "")
    return dict(kind="SNV", name=name, floor=floor, dbsnp=dbsnp, kg=kg)


def indel_rule(other, hp=None, h=None):
    name = f"AF>={FTXT[other]} (every INDEL)" if h is None else f"HP>={h}: AF>={FTXT[hp]}, else AF>={FTXT[other]}"
    return dict(kind="INDEL", name=name, other=other, hp=hp, h=h)


SNV_GRID = [snv_rule(f, d, k) for f in (1e-5, 1e-4, 1e-3, 1e-2) for d in (False, True) for k in (False, True)]
INDEL_GRID = [indel_rule(1e-3), indel_rule(1e-2)] + [indel_rule(o, hp, h) for h in (6, 7, 8)
                                                     for o, hp in ((1e-3, 1e-2), (1e-3, 5e-2), (1e-4, 1e-2))]
GRID = {"SNV": SNV_GRID, "INDEL": INDEL_GRID}
CURRENT = {"SNV": snv_rule(1e-4, True, True), "INDEL": indel_rule(1e-3, 1e-2, 7)}
NO_PON = {"SNV": dict(kind="SNV", name="no PoN", none=True), "INDEL": dict(kind="INDEL", name="no PoN", none=True)}
assert CURRENT["SNV"]["name"] in [r["name"] for r in SNV_GRID]
assert CURRENT["INDEL"]["name"] in [r["name"] for r in INDEL_GRID]


def allele_tagged(rule, ft):
    if rule.get("none"):
        return False
    if rule["kind"] == "SNV":
        fi = FI[rule["floor"]]
        return ft[G_LVL] >= fi or ft[C_LVL] >= fi or (rule["dbsnp"] and ft[D_NONSOM] == 1) or (rule["kg"] and ft[K_ANY] == 1)
    floor = rule["hp"] if rule["h"] is not None and ft[HP] >= rule["h"] else rule["other"]
    return ft[G_LVL] >= FI[floor] or ft[C_LVL] >= FI[floor]


def record_tagged(rule, pon_kind, fts):
    """A record of the script's kind pon_kind with ALT features fts; the other kind keeps its current rule."""
    r = rule if rule["kind"] == pon_kind else CURRENT[pon_kind]
    return any(allele_tagged(r, ft) for ft in fts)


# ------------------------------------------------------------------ load
def load():
    tens = defaultdict(Counter)    # (set, kind, group, label) -> Counter(feature tuple | 'unplaced' ...)
    reasons = Counter()            # (set, kind, group, label, reason) -> n
    truth = defaultdict(Counter)   # (sample, kind, group, in_bed) -> Counter((pon_kind, fts))
    checks = dict(alleles=0, mismatches=0, distinct_alleles=0, missing=[], seconds=0.0)
    for chrom in CHROMS:
        d = json.loads((DATA / "scan" / f"{chrom}.json").read_text())
        group = TEST if chrom == TEST else TRAIN
        fts = [tuple(f) for f in d["feature_tuples"]]
        for s, kind, label, reason, fid, n in d["tensors"]:
            tens[(s, kind, group, label)][fts[fid] if isinstance(fid, int) else fid] += n
            reasons[(s, kind, group, label, reason)] += n
        for sample, pon_kind, kinds, in_bed, fids, n in d["truth"]:
            for kind in kinds:
                truth[(sample, kind, group, in_bed)][(pon_kind, tuple(fts[i] for i in fids))] += n
        checks["alleles"] += d["check"]["alleles"]
        checks["mismatches"] += d["check"]["n_mismatches"]
        checks["distinct_alleles"] += d["distinct_alleles"]
        checks["missing"] += d["missing_label_files"]
        checks["seconds"] += d["seconds"]
    return tens, reasons, truth, checks


# ------------------------------------------------------------------ metrics
def tensor_counts(tens, s, kind, group, rule):
    out = {}
    for label, name in ((1, "som"), (2, "germ"), (0, "non")):
        c = tens.get((s, kind, group, label), Counter())
        out[f"{name}_n"] = sum(c.values())
        out[f"{name}_tag"] = sum(v for k, v in c.items() if isinstance(k, tuple) and allele_tagged(rule, k))
        out[f"{name}_unplaced"] = sum(v for k, v in c.items() if not isinstance(k, tuple))
    return out


def scores(c):
    s = 1 - c["som_tag"] / c["som_n"] if c["som_n"] else float("nan")
    neg_n, neg_tag = c["germ_n"] + c["non_n"], c["germ_tag"] + c["non_tag"]
    n = 1 - neg_tag / neg_n
    ng = 1 - c["germ_tag"] / c["germ_n"]
    tp, fn, fp = c["som_n"] - c["som_tag"], c["som_tag"], neg_n - neg_tag
    return dict(som_lost=1 - s, germ_removed=c["germ_tag"] / c["germ_n"], non_removed=c["non_tag"] / c["non_n"],
                neg_removed=neg_tag / neg_n, f1p=2 * s / (1 + s + n), f1p_germ=2 * s / (1 + s + ng),
                f1_count=2 * tp / (2 * tp + fp + fn))


def truth_counts(truth, sample, kind, group, rule):
    c = truth.get((sample, kind, group, True), Counter())
    n = sum(c.values())
    tag = sum(v for (pk, fts), v in c.items() if record_tagged(rule, pk, fts))
    return dict(n=n, tag=tag, lost=tag / n if n else float("nan"))


def evaluate(tens, truth):
    """res[(kind, group, set, rule name)] = counts + scores; tru[(kind, group, sample, rule name)]."""
    res, tru = {}, {}
    for kind in ("SNV", "INDEL"):
        for rule in GRID[kind] + [NO_PON[kind]]:
            for group in (TRAIN, TEST):
                for s in D9 + D46:
                    c = tensor_counts(tens, s, kind, group, rule)
                    res[(kind, group, s, rule["name"])] = dict(c, **scores(c))
                for sample in SAMPLES:
                    tru[(kind, group, sample, rule["name"])] = truth_counts(truth, sample, kind, group, rule)
    return res, tru


def select(res, kind, group, sets, score="f1p"):
    """(selected rule name, rows sorted by mean score) under the objective in the module docstring."""
    rows = []
    for rule in GRID[kind]:
        vals = [res[(kind, group, s, rule["name"])] for s in sets]
        rows.append(dict(rule=rule["name"], mean=sum(v[score] for v in vals) / len(vals),
                         som_lost=sum(v["som_lost"] for v in vals) / len(vals)))
    best = max(r["mean"] for r in rows)
    tied = [r for r in rows if r["mean"] >= best - TIE]
    chosen = min(tied, key=lambda r: (r["som_lost"], -r["mean"]))["rule"]
    rows.sort(key=lambda r: -r["mean"])
    for rank, r in enumerate(rows, 1):
        r["rank"] = rank
    return chosen, rows, [r["rule"] for r in tied]


def per_set_winner(res, kind, group, s, score="f1p"):
    vals = [(res[(kind, group, s, r["name"])][score], -res[(kind, group, s, r["name"])]["som_lost"], r["name"])
            for r in GRID[kind]]
    best = max(v[0] for v in vals)
    return max(v for v in vals if v[0] >= best - TIE)[2]


# ------------------------------------------------------------------ pipeline cross-check (chr1 truth tags)
def pipeline_truth_check(truth):
    """The chr1 truth records and PoN-tagged truth records of every vcfeval report under RUNS vs this scan (records
    of the kind, PASS or '.', the whole chromosome: vcfeval's truth_pon_tags)."""
    rows = []
    seen = set()
    for rep in sorted(RUNS.glob("*/vcf_chr1*/eval_bed/report.json")):
        r = json.loads(rep.read_text())
        t = r.get("truth_pon_tags")
        if not t:
            continue
        sample = "HG008T" if "HG008" in r["truth"] else "COLO829T" if "COLO829T" in r["truth"] else None
        if sample is None:
            continue
        key = (sample, r["kind"], t["records"], t["tagged"])
        if key in seen:
            continue
        seen.add(key)
        c = Counter()
        for in_bed in (True, False):
            c.update(truth.get((sample, r["kind"], TEST, in_bed), Counter()))
        mine_n = sum(c.values())
        mine_tag = sum(v for (pk, fts), v in c.items() if record_tagged(CURRENT[r["kind"]], pk, fts))
        rows.append(dict(report=str(rep.parent.parent.parent.name), sample=sample, kind=r["kind"],
                         pipeline_records=t["records"], pipeline_tagged=t["tagged"], scan_records=mine_n,
                         scan_tagged=mine_tag))
    return rows


# ------------------------------------------------------------------ repeat context (INDEL)
def repeat_context(tens, truth, rules):
    """Per RepeatContext class: truth records lost and negative tensors removed (d9 sets pooled) under rules."""
    classes = ("HP>=7", "HP4-6", "STR", "VNTR", "none", "complex")
    out = {}
    for group in (TRAIN, TEST):
        for cls in classes:
            row = {}
            for sample in SAMPLES:
                c = truth.get((sample, "INDEL", group, True), Counter())
                sel = {k: v for k, v in c.items() if any(ft[RC] == cls for ft in k[1])}
                row[f"{sample}_n"] = sum(sel.values())
                for rule in rules:
                    row[f"{sample}_{rule['name']}"] = sum(v for (pk, fts), v in sel.items() if record_tagged(rule, pk, fts))
            for lab, name in ((2, "germ"), (0, "non"), (1, "som")):
                n = Counter()
                for s in D9:
                    for ft, v in tens.get((s, "INDEL", group, lab), Counter()).items():
                        if isinstance(ft, tuple) and ft[RC] == cls:
                            n[ft] += v
                row[f"{name}_n"] = sum(n.values())
                for rule in rules:
                    row[f"{name}_{rule['name']}"] = sum(v for ft, v in n.items() if allele_tagged(rule, ft))
            out[(group, cls)] = row
    # homopolymer_run >= 7 vs RepeatContext HP>=7 on the truth INDEL alleles (both samples, chr1-22, in the BED)
    agree = Counter()
    for (sample, kind, group, in_bed), c in truth.items():
        if kind != "INDEL" or not in_bed:
            continue
        for (pk, fts), v in c.items():
            for ft in fts:
                if ft[RC] != "SNV/MNP":
                    agree[(ft[HP] >= 7, ft[RC] == "HP>=7")] += v
    return out, classes, agree


# ------------------------------------------------------------------ formatting
def pct(x, d=1):
    return "-" if x != x else f"{100 * x:.{d}f}"


def f3(x):
    return "-" if x != x else f"{x:.4f}"


def mark(name, kind, chosen):
    m = []
    if name == CURRENT[kind]["name"]:
        m.append("current")
    if name == chosen:
        m.append("**selected**")
    return f"{name}" + (f" ({', '.join(m)})" if m else "")


def main():
    t0 = time.time()
    tens, reasons, truth, checks = load()
    res, tru = evaluate(tens, truth)
    sel = {}
    for kind in ("SNV", "INDEL"):
        for group in (TRAIN, TEST):
            for arm, sets in (("d9", D9), ("d46", D46)):
                for score in ("f1p", "f1p_germ", "f1_count"):
                    sel[(kind, group, arm, score)] = select(res, kind, group, sets, score)
    L = []
    w = L.append
    w("# E17c: the PoN rule re-selected on chr2-22, then checked on chr1")
    w("")
    w(f"Generated by `e17c_tables.py` ({time.strftime('%Y-%m-%d %H:%M')}) from `e17c_scan.py` outputs "
      f"(`{DATA}/scan/chr*.json`). Objective and grid: the docstring of `e17c_tables.py` (fixed before any grid number "
      "was computed). Percentages are of the tensors / truth records of the row's chromosomes.")
    w("")
    w("Terms. *som lost*: somatic-labelled tensors (label 1) the rule tags. *germ removed* / *non removed*: germline "
      "(label 2) / non (label 0) tensors it tags. *truth lost*: somatic truth records of the kind in the benchmark BED "
      "it tags (same for every platform of a sample). *F1 proxy* = 2s / (1 + s + n), s = somatic tensors kept, "
      "n = negative tensors (germline + non) kept.")
    w("")

    # ---- 0. checks
    w("## 0. Inputs and checks")
    w("")
    w(f"- Scan: chr1-22, {checks['distinct_alleles']:,} distinct GRCh38 alleles (per chromosome, summed), "
      f"{checks['seconds'] / 60:.0f} CPU-min of the scan's main process; label files missing: "
      f"{len(checks['missing'])}.")
    w(f"- In-memory PoN queries vs the tabix files (record_matches on both, random alleles, half at a PoN position): "
      f"{checks['alleles']:,} alleles, {checks['mismatches']} differences.")
    cs = json.loads((DATA / "check_script_chr22.json").read_text()) if (DATA / "check_script_chr22.json").exists() else None
    if cs:
        w(f"- The unmodified `scripts/filter_panel_of_normals.py` on {cs['alleles']:,} chr22 tensor alleles "
          f"({cs['by_kind']['SNV']:,} SNV, {cs['by_kind']['INDEL']:,} INDEL; {cs['script_tagged']:,} tagged): "
          f"{cs['differences']} differences from the current rule computed from the scan features, per PoN "
          "(`e17c_check_script.py`).")
    pc = pipeline_truth_check(truth)
    if pc:
        w("- chr1 truth records tagged, the pipeline's vcfeval (`report.json` truth_pon_tags) vs this scan, current rule:")
        w("")
        w("| sample | kind | report (first with these numbers) | pipeline records / tagged | scan records / tagged |")
        w("|---|---|---|---:|---:|")
        for r in pc:
            w(f"| {r['sample']} | {r['kind']} | {r['report']} | {r['pipeline_records']:,} / {r['pipeline_tagged']:,} | "
              f"{r['scan_records']:,} / {r['scan_tagged']:,} |")
    w("")
    w("Tensors of the six d9 sets, chr2-22 / chr1 (label: n, of them unplaced = no GRCh38 allele, never tagged):")
    w("")
    w("| set | kind | somatic chr2-22 | germline chr2-22 | non chr2-22 | somatic chr1 | germline chr1 | non chr1 |")
    w("|---|---|---:|---:|---:|---:|---:|---:|")
    for kind in ("SNV", "INDEL"):
        for s in D9 + D46:
            cells = []
            for group in (TRAIN, TEST):
                c = res[(kind, group, s, NO_PON[kind]["name"])]
                for name in ("som", "germ", "non"):
                    cells.append(f"{c[f'{name}_n']:,} ({c[f'{name}_unplaced']:,})")
            w(f"| {s} | {kind} | " + " | ".join(cells) + " |")
    w("")
    w("Somatic truth records in the benchmark BED: " + "; ".join(
        f"{sample} {kind} chr2-22 {tru[(kind, TRAIN, sample, 'no PoN')]['n']:,}, chr1 {tru[(kind, TEST, sample, 'no PoN')]['n']:,}"
        for kind in ("SNV", "INDEL") for sample in SAMPLES) + ".")
    w("")

    # ---- per kind
    sec = 1
    for kind in ("SNV", "INDEL"):
        chosen_train = sel[(kind, TRAIN, "d9", "f1p")][0]
        chosen_test = sel[(kind, TEST, "d9", "f1p")][0]
        for group, chosen in ((TRAIN, chosen_train), (TEST, chosen_test)):
            _, rows, tied = sel[(kind, group, "d9", "f1p")]
            w(f"## {sec}. {kind}, {group}: rule grid, six d9 sets" + (" (selection)" if group == TRAIN else " (check)"))
            sec += 1
            w("")
            w(f"Selected on {group}: **{chosen}**" + (f" (tied within {TIE:g}: {', '.join(tied)})" if len(tied) > 1 else "")
              + f". Current rule: {CURRENT[kind]['name']}, rank "
              f"{next(r['rank'] for r in rows if r['rule'] == CURRENT[kind]['name'])} of {len(rows)}.")
            w("")
            w("Means over the six sets (somatic, germline, non, F1 proxy); truth lost per sample.")
            w("")
            w(f"| rank | rule | {group} truth lost HG008T % | {group} truth lost COLO829T % | {group} som lost % | "
              f"{group} germ removed % | {group} non removed % | **{group} F1 proxy** | {group} F1 proxy, germline only | "
              f"{group} candidate F1 |")
            w("|---:|---|---:|---:|---:|---:|---:|---:|---:|---:|")
            order = rows + [dict(rule=NO_PON[kind]["name"], rank="-")]
            for r in order:
                vals = [res[(kind, group, s, r["rule"])] for s in D9]
                mean = {k: sum(v[k] for v in vals) / len(vals) for k in ("som_lost", "germ_removed", "non_removed",
                                                                          "f1p", "f1p_germ", "f1_count")}
                w(f"| {r['rank']} | {mark(r['rule'], kind, chosen)} | {pct(tru[(kind, group, 'HG008T', r['rule'])]['lost'])} | "
                  f"{pct(tru[(kind, group, 'COLO829T', r['rule'])]['lost'])} | {pct(mean['som_lost'])} | "
                  f"{pct(mean['germ_removed'])} | {pct(mean['non_removed'])} | **{f3(mean['f1p'])}** | "
                  f"{f3(mean['f1p_germ'])} | {f3(mean['f1_count'])} |")
            w("")
            w(f"{group} F1 proxy per set (per-set best in bold):")
            w("")
            w("| rule | " + " | ".join(f"{SHORT[s]} {group}" for s in D9) + f" | mean {group} |")
            w("|---|" + "---:|" * (len(D9) + 1))
            winners = {s: per_set_winner(res, kind, group, s) for s in D9}
            for r in rows:
                cells = []
                for s in D9:
                    v = f3(res[(kind, group, s, r["rule"])]["f1p"])
                    cells.append(f"**{v}**" if winners[s] == r["rule"] else v)
                w(f"| {mark(r['rule'], kind, chosen)} | " + " | ".join(cells) + f" | {f3(r['mean'])} |")
            w("")
            w(f"Per set, {group}: som lost % / germ removed % / non removed %:")
            w("")
            w("| rule | " + " | ".join(f"{SHORT[s]} {group}" for s in D9) + " |")
            w("|---|" + "---|" * len(D9))
            for r in rows + [dict(rule=NO_PON[kind]["name"])]:
                cells = []
                for s in D9:
                    v = res[(kind, group, s, r["rule"])]
                    cells.append(f"{pct(v['som_lost'])} / {pct(v['germ_removed'])} / {pct(v['non_removed'])}")
                w(f"| {mark(r['rule'], kind, chosen)} | " + " | ".join(cells) + " |")
            w("")

    # ---- selection summary + sensitivity
    w(f"## {sec}. Selections, sensitivity and replication")
    sec += 1
    w("")
    w("| kind | chromosomes | sets | objective | selected | current rule's rank | per-set winners |")
    w("|---|---|---|---|---|---:|---|")
    names = {"f1p": "F1 proxy (primary)", "f1p_germ": "F1 proxy, germline only", "f1_count": "candidate F1"}
    for kind in ("SNV", "INDEL"):
        for group in (TRAIN, TEST):
            for arm, sets in (("d9", D9), ("d46", D46)):
                for score in ("f1p", "f1p_germ", "f1_count"):
                    chosen, rows, _ = sel[(kind, group, arm, score)]
                    rank = next(r["rank"] for r in rows if r["rule"] == CURRENT[kind]["name"])
                    win = Counter(per_set_winner(res, kind, group, s, score) for s in sets)
                    w(f"| {kind} | {group} | {arm} | {names[score]} | {chosen} | {rank} of {len(rows)} | "
                      + "; ".join(f"{k} x{v}" for k, v in win.most_common()) + " |")
    w("")
    w("Exchange rate, chr2-22, d9: negative tensors removed per somatic tensor lost when moving from the selected rule "
      "to the current one (positive both = the current rule removes more and loses more).")
    w("")
    w("| kind | set | selected | Δ somatic lost | Δ germline removed | Δ non removed | negatives per somatic |")
    w("|---|---|---|---:|---:|---:|---:|")
    for kind in ("SNV", "INDEL"):
        chosen = sel[(kind, TRAIN, "d9", "f1p")][0]
        for s in D9:
            a, b = res[(kind, TRAIN, s, chosen)], res[(kind, TRAIN, s, CURRENT[kind]["name"])]
            ds, dg, dn = b["som_tag"] - a["som_tag"], b["germ_tag"] - a["germ_tag"], b["non_tag"] - a["non_tag"]
            rate = f"{(dg + dn) / ds:.1f}" if ds else "-"
            w(f"| {kind} | {SHORT[s]} | {chosen} | {ds:+,} | {dg:+,} | {dn:+,} | {rate} |")
    w("")

    # ---- repeat context
    chosen_i = sel[("INDEL", TRAIN, "d9", "f1p")][0]
    rules = [CURRENT["INDEL"]] + ([r for r in INDEL_GRID if r["name"] == chosen_i and r["name"] != CURRENT["INDEL"]["name"]])
    rules += [r for r in INDEL_GRID if r["name"] in ("AF>=1e-3 (every INDEL)", "AF>=0.01 (every INDEL)")]
    rc, classes, agree = repeat_context(tens, truth, rules)
    w(f"## {sec}. INDEL by repeat context (indexed_gam_pipeline_v4 RepeatContext)")
    sec += 1
    w("")
    w("Truth lost = truth INDEL records in the BED tagged; germ / non removed = INDEL tensors of the six d9 sets "
      "(pooled) tagged. Class of a record = any ALT of the class.")
    w("")
    for group in (TRAIN, TEST):
        w(f"{group}:")
        w("")
        head = ["class", "HG008T truth n", "COLO829T truth n", "germ n", "non n"]
        for rule in rules:
            head += [f"{rule['name']}: truth lost HG008T / COLO829T %, germ / non removed %"]
        w("| " + " | ".join(head) + " |")
        w("|---|" + "---:|" * 4 + "---|" * len(rules))
        for cls in classes:
            r = rc[(group, cls)]
            cells = [cls, f"{r['HG008T_n']:,}", f"{r['COLO829T_n']:,}", f"{r['germ_n']:,}", f"{r['non_n']:,}"]
            for rule in rules:
                q = lambda k, n: pct(r[f"{k}_{rule['name']}"] / r[n]) if r[n] else "-"  # noqa: E731
                cells.append(f"{q('HG008T', 'HG008T_n')} / {q('COLO829T', 'COLO829T_n')}, "
                             f"{q('germ', 'germ_n')} / {q('non', 'non_n')}")
            w("| " + " | ".join(cells) + " |")
        w("")
    w("Homopolymer of the script (`homopolymer_run` >= 7) vs RepeatContext HP>=7, truth INDEL alleles in the BED, "
      "chr1-22: " + ", ".join(f"run>=7 {a} / RC HP>=7 {b}: {agree[(a, b)]:,}" for a in (True, False) for b in (True, False)) + ".")
    w("")

    # ---- d46 replication table (compact)
    w(f"## {sec}. Replication on the HPRC v2.1 d46 sets (F1 proxy, mean over the six sets)")
    sec += 1
    w("")
    for kind in ("SNV", "INDEL"):
        w(f"| {kind} rule | chr2-22 mean | chr2-22 rank | chr1 mean | chr1 rank |")
        w("|---|---:|---:|---:|---:|")
        a = {r["rule"]: r for r in sel[(kind, TRAIN, "d46", "f1p")][1]}
        b = {r["rule"]: r for r in sel[(kind, TEST, "d46", "f1p")][1]}
        for r in sel[(kind, TRAIN, "d46", "f1p")][1]:
            nm = r["rule"]
            w(f"| {mark(nm, kind, sel[(kind, TRAIN, 'd46', 'f1p')][0])} | {f3(a[nm]['mean'])} | {a[nm]['rank']} | "
              f"{f3(b[nm]['mean'])} | {b[nm]['rank']} |")
        w("")

    (HERE / "tables.md").write_text("\n".join(L) + "\n")
    summary = dict(
        generated=time.strftime("%Y-%m-%d %H:%M"), objective="F1 proxy = 2s/(1+s+n), mean over the six d9 sets",
        current={k: v["name"] for k, v in CURRENT.items()},
        selected={f"{kind}|{group}|{arm}|{score}": v[0] for (kind, group, arm, score), v in sel.items()},
        ranks={f"{kind}|{group}|{arm}|{score}": {r["rule"]: dict(rank=r["rank"], mean=r["mean"], som_lost=r["som_lost"])
                                                 for r in v[1]} for (kind, group, arm, score), v in sel.items()},
        tensors={"|".join(k): v for k, v in res.items()},
        truth={"|".join(k): v for k, v in tru.items()},
        checks=checks, pipeline_truth_check=pc, check_script=cs,
        hp_vs_repeatcontext={f"run>=7={a}|RC_HP>=7={b}": agree[(a, b)] for a in (True, False) for b in (True, False)})
    (HERE / "summary.json").write_text(json.dumps(summary, indent=1) + "\n")
    print(f"tables.md, summary.json written in {time.time() - t0:.0f} s")
    for kind in ("SNV", "INDEL"):
        print(kind, "selected chr2-22:", sel[(kind, TRAIN, "d9", "f1p")][0], "| chr1:", sel[(kind, TEST, "d9", "f1p")][0],
              "| current:", CURRENT[kind]["name"])


if __name__ == "__main__":
    main()
