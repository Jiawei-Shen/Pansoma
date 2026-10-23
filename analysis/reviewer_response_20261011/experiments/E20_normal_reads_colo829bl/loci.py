"""E20 step 1: the COLO829T absorbed truth alleles on HPRC v1.1 (d9) and HPRC v2.1 (d46), and the GRCh38 regions whose
reads are extracted from the COLO829BL (normal) and COLO829T (tumor) ONT BAMs.

Loci
  d9   the 497 absorbed alleles of analysis/graph_absorbed_somatic_colo829t_20261005 (variant_set.tsv: no candidate on
       >= 1 platform and the majority of the ALT-like reads walk a graph branch / skip edge with no edit), with that
       analysis' categories (hprc_per_variant.tsv: germline filtering / pangenome-induced false negative / ambiguous)
       and patient-frame classes (normal_frame.tsv).
  d46  E5d's colo829t_d46_absorbed_loci.tsv when it exists (E5b's 'absorbed' rule, the set R1.4 quotes: no candidate
       on d46 and a tensor on the GRCh38-only graph on >= 1 platform, and a d46 path spells GRCh38 + truth); otherwise
       the E5c analogue: a truth allele lost on d46
       against the GRCh38-only graph whose read-level class is I1 (ALT reads entirely on a graph alt path, no edit) on
       >= 1 of the three COLO829T d46 sets (E5c DATA/per_locus.tsv.gz). `d46_source` in the output says which.
Regions: every locus (d9, d46, and every E5c COLO829T `lost` locus, so another d46 definition needs no new extraction)
+- FLANK bp around the truth allele's span, merged.

Writes $D/loci_d9.tsv, $D/loci_d46.tsv, $D/regions.bed (D = /scratch/jshen/data/reviewer_experiments_20261013/
E20_normal_reads_colo829bl) and loci_summary.json next to this script.
usage: /wanglab/jshen/anaconda3/bin/python loci.py
"""
import csv
import gzip
import json
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
D = Path("/scratch/jshen/data/reviewer_experiments_20261013/E20_normal_reads_colo829bl")
ABS = Path("/scratch/jshen/data/pansoma_net_v2_runs/graph_absorbed_somatic_colo829t_20261005")
E5C = Path("/scratch/jshen/data/reviewer_experiments_20261013/E5c_d46_loss_reads/per_locus.tsv.gz")
E5D = HERE.parent / "E5d_absorbed_normal_d46" / "colo829t_d46_absorbed_loci.tsv"
TRUTH = Path("/scratch/jshen/data/pansoma_v2_tensors/COLO829T_ONT_d46/truth/somatic.graph.tsv")
FLANK = 500
csv.field_size_limit(sys.maxsize)


def read_tsv(path):
    opener = gzip.open if str(path).endswith(".gz") else open
    with opener(path, "rt") as f:
        return list(csv.DictReader(f, delimiter="\t"))


def main():
    truth = {r["truth_id"]: r for r in read_tsv(TRUTH)}
    hprc = {r["truth_id"]: r for r in read_tsv(ABS / "hprc_per_variant.tsv")}
    frame = {r["truth_id"]: r for r in read_tsv(ABS / "normal_frame.tsv")}
    paths = {r["truth_id"]: r for r in read_tsv(ABS / "graph_paths.tsv")}
    fields = ["truth_id", "chrom", "pos0", "ref", "alt", "kind", "vcf_pos", "vcf_ref", "vcf_alt"]

    d9 = []
    for r in read_tsv(ABS / "variant_set.tsv"):
        t = truth[r["truth_id"]]
        assert (t["chrom"], t["vcf_pos"], t["vcf_ref"], t["vcf_alt"]) == (r["chrom"], r["vcf_pos"], r["vcf_ref"], r["vcf_alt"])
        h, f, g = hprc.get(r["truth_id"], {}), frame.get(r["truth_id"], {}), paths.get(r["truth_id"], {})
        d9.append(dict({k: t[k] for k in fields}, kind2=r["kind2"], RGN=r["RGN"], VAF_Ill=r["VAF_Ill"],
                       ONT_perfect=r["ONT_perfect"], d9_match=g.get("match", ""),
                       category=h.get("category", ""), sub_class=h.get("sub_class", ""),
                       pf_label=f.get("pf_label", "")))

    e5c = read_tsv(E5C)
    lost = [r for r in e5c if r["set"].startswith("COLO829T") and r["locus_set"] == "lost"]
    if E5D.exists():
        src = str(E5D)
        d46_ids = {r["truth_id"]: r for r in read_tsv(E5D)}
    else:
        src = "E5c I1 on >= 1 COLO829T d46 set (E5d file absent)"
        d46_ids = {}
        for r in lost:
            if r["cls"].startswith("I1"):
                d46_ids.setdefault(r["truth_id"], {}).setdefault("platforms_I1", []).append(r["set"].split("_")[1])
    d46 = []
    for tid, extra in sorted(d46_ids.items(), key=lambda kv: int(kv[0])):
        t = truth[tid]
        d46.append(dict({k: t[k] for k in fields}, kind2="SNV" if t["kind"] == "SNP" else "INDEL",
                        in_d9_set=str(tid in {x["truth_id"] for x in d9}),
                        platforms=extra.get("platforms") or ",".join(extra.get("platforms_I1", [])),
                        pon=extra.get("pon", ""), elements=extra.get("elements", ""), d46_source=src))

    spans = []
    for tid in {r["truth_id"] for r in d9} | set(d46_ids) | {r["truth_id"] for r in lost}:
        t = truth[tid]
        lo = int(t["pos0"])
        hi = lo + max(1, len(t["ref"]))
        spans.append((t["chrom"], max(0, lo - FLANK), hi + FLANK))
    spans.sort(key=lambda s: (int(s[0][3:]) if s[0][3:].isdigit() else 99, s[1]))
    merged = []
    for c, s, e in spans:
        if merged and merged[-1][0] == c and s <= merged[-1][2]:
            merged[-1][2] = max(merged[-1][2], e)
        else:
            merged.append([c, s, e])

    D.mkdir(parents=True, exist_ok=True)
    for name, rows in (("loci_d9.tsv", d9), ("loci_d46.tsv", d46)):
        with open(D / name, "w", newline="") as fh:
            w = csv.DictWriter(fh, fieldnames=list(rows[0]), delimiter="\t", lineterminator="\n")
            w.writeheader()
            w.writerows(rows)
    with open(D / "regions.bed", "w") as fh:
        for c, s, e in merged:
            fh.write(f"{c}\t{s}\t{e}\n")
    summary = dict(d9=len(d9), d9_snv=sum(r["kind2"] == "SNV" for r in d9), d46=len(d46),
                   d46_snv=sum(r["kind2"] == "SNV" for r in d46), d46_source=src,
                   d46_in_d9=sum(r["in_d9_set"] == "True" for r in d46),
                   region_loci=len(spans), regions=len(merged), region_bp=sum(e - s for _, s, e in merged))
    (HERE / "loci_summary.json").write_text(json.dumps(summary, indent=1) + "\n")
    print(json.dumps(summary, indent=1))


if __name__ == "__main__":
    main()
