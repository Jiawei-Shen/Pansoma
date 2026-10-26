"""E9i report: merge the shards of normal_reads.py and tabulate matched-normal read evidence per COLO829T run.

Per run and group (unplaced PASS calls; placed controls ctrl_tp = true somatic, ctrl_fp = false positives after the PoN):
the share 'ALT in normal' (>= 2 COLO829BL ONT reads with the ALT 21-mer), 'one normal read', 'absent from normal'
(0 ALT reads, >= 5 REF reads) and 'low normal coverage', and among the absent ones how many have the ALT in >= 2 tumor
reads (candidate somatic). For the unplaced calls also by E9c category (population / germline / not_in_normal /
unverifiable; `../E9c_offref_filter_validation/per_call.tsv.gz`).
Writes tables.md and summary.json next to this script. usage: python3 report.py
"""
import csv
import glob
import gzip
import json
from collections import Counter, defaultdict
from pathlib import Path

HERE = Path(__file__).resolve().parent
DATA = Path("/scratch/jshen/data/reviewer_experiments_20261013/E9i_normal_reads_offref")
E9C = HERE.parent / "E9c_offref_filter_validation" / "per_call.tsv.gz"
STATUS = ("ALT in normal", "one normal read", "absent from normal", "low normal coverage", "no k-mer")
GROUPS = ("unplaced_pass", "ctrl_tp", "ctrl_fp")


def main():
    cat = {}
    with gzip.open(E9C, "rt") as f:
        for r in csv.DictReader(f, delimiter="\t"):
            cat[(r["run"], r["group"], r["candidate_id"])] = r["category"]
    runs = defaultdict(list)
    for p in sorted(glob.glob(str(DATA / "*.tsv"))):
        name = Path(p).name
        run = name.split(".")[0]
        with open(p) as f:
            runs[run] += list(csv.DictReader(f, delimiter="\t"))
    summary, L = {}, ["# E9i: matched-normal reads at COLO829T graph-coordinate SNV calls (chr1)", "",
                      "ALT / REF 21-mer of each call counted in COLO829BL ONT reads (D001 + D003) and COLO829T ONT reads over the call's "
                      "GRCh38 locus +- 500 bp. 'absent' = 0 normal ALT reads with >= 5 REF reads; 'somatic candidate' = absent from the "
                      "normal and ALT in >= 2 tumor reads.", "",
                      "| run | group | calls | ALT in normal | one normal read | absent from normal | low coverage / no k-mer | somatic candidates |",
                      "|---|---|---:|---|---|---|---|---|"]
    for run in sorted(runs):
        rows = runs[run]
        summary[run] = {}
        for g in GROUPS:
            rs = [r for r in rows if r["group"] == g]
            if not rs:
                continue
            c = Counter(r["status"] for r in rs)
            som = sum(r["status"] == "absent from normal" and int(r["tumor_alt"] or 0) >= 2 for r in rs)
            n = len(rs)
            by_cat = defaultdict(Counter)
            if g == "unplaced_pass":
                for r in rs:
                    k = cat.get((run, g, r["candidate_id"]), "no category")
                    by_cat[k][r["status"]] += 1
                    by_cat[k]["somatic candidate"] += r["status"] == "absent from normal" and int(r["tumor_alt"] or 0) >= 2
            summary[run][g] = dict(calls=n, status=dict(c), somatic_candidates=som, by_category={k: dict(v) for k, v in by_cat.items()})
            pct = lambda x: f"{x} ({100 * x / n:.1f}%)"  # noqa: E731
            L.append(f"| {run.replace('rcpu__', '').replace('__vcf_chr1', '')} | {g} | {n} | {pct(c['ALT in normal'])} | "
                     f"{pct(c['one normal read'])} | {pct(c['absent from normal'])} | "
                     f"{pct(c['low normal coverage'] + c['no k-mer'])} | {pct(som)} |")
    L += ["", "## Unplaced PASS calls by E9c category", "",
          "| run | category | calls | ALT in normal | absent from normal | somatic candidates |", "|---|---|---:|---:|---:|---:|"]
    for run in sorted(summary):
        for k, v in sorted(summary[run].get("unplaced_pass", {}).get("by_category", {}).items()):
            n = sum(v.get(s, 0) for s in STATUS)
            L.append(f"| {run.replace('rcpu__', '').replace('__vcf_chr1', '')} | {k} | {n} | {v.get('ALT in normal', 0)} | "
                     f"{v.get('absent from normal', 0)} | {v.get('somatic candidate', 0)} |")
    (HERE / "tables.md").write_text("\n".join(L) + "\n")
    (HERE / "summary.json").write_text(json.dumps(summary, indent=1) + "\n")
    print("\n".join(L))


if __name__ == "__main__":
    main()
