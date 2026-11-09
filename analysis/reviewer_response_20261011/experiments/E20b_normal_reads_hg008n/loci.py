"""E20b step 1: the HG008T absorbed truth alleles on HPRC v2.1 (d46) and the GRCh38 regions whose reads are taken from
the HG008-N (normal) and HG008-T (tumor) PacBio HiFi BAMs.

Loci = E5d's hg008t_d46_absorbed_loci.tsv (1,607: 140 SNV / 1,467 INDEL; E5b's 'absorbed' rule, the set R1.4 / R2.3
quote: no candidate on d46 and a tensor on the GRCh38-only graph on >= 1 platform, and a d46 path spells GRCh38 +
truth), with the truth records of the HG008T d46 tensor sets (identical truth table on the three platforms).
Regions: every locus +- FLANK bp around the truth allele's span, merged (as E20's loci.py).

Writes $D/loci_d46.tsv, $D/regions.bed (D = /scratch/jshen/data/reviewer_experiments_20261013/E20b_normal_reads_hg008n)
and loci_summary.json next to this script.
usage: /wanglab/jshen/anaconda3/bin/python loci.py
"""
import csv
import json
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
D = Path("/scratch/jshen/data/reviewer_experiments_20261013/E20b_normal_reads_hg008n")
E5D = HERE.parent / "E5d_absorbed_normal_d46" / "hg008t_d46_absorbed_loci.tsv"
TRUTH = Path("/scratch/jshen/data/pansoma_v2_tensors/HG008T_PacBio_d46/truth/somatic.graph.tsv")
FLANK = 500
csv.field_size_limit(sys.maxsize)


def read_tsv(path):
    with open(path) as f:
        return list(csv.DictReader(f, delimiter="\t"))


def main():
    truth = {r["truth_id"]: r for r in read_tsv(TRUTH)}
    fields = ["truth_id", "chrom", "pos0", "ref", "alt", "kind", "vcf_pos", "vcf_ref", "vcf_alt"]
    rows = []
    for e in read_tsv(E5D):
        t = truth[e["truth_id"]]
        assert (t["chrom"], t["vcf_pos"], t["vcf_ref"], t["vcf_alt"]) == (e["chrom"], e["vcf_pos"], e["vcf_ref"], e["vcf_alt"])
        kind2 = "SNV" if t["kind"] == "SNP" else "INDEL"
        assert kind2 == e["kind2"], (e["truth_id"], t["kind"], e["kind2"])
        rows.append(dict({k: t[k] for k in fields}, kind2=kind2, platforms=e["platforms"], pon=e["pon"],
                         elements=e["elements"], d46_source=str(E5D)))
    rows.sort(key=lambda r: int(r["truth_id"]))

    spans = []
    for t in rows:
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
    with open(D / "loci_d46.tsv", "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(rows[0]), delimiter="\t", lineterminator="\n")
        w.writeheader()
        w.writerows(rows)
    with open(D / "regions.bed", "w") as fh:
        for c, s, e in merged:
            fh.write(f"{c}\t{s}\t{e}\n")
    summary = dict(d46=len(rows), d46_snv=sum(r["kind2"] == "SNV" for r in rows),
                   d46_indel=sum(r["kind2"] == "INDEL" for r in rows), d46_source=str(E5D),
                   regions=len(merged), region_bp=sum(e - s for _, s, e in merged))
    (HERE / "loci_summary.json").write_text(json.dumps(summary, indent=1) + "\n")
    print(json.dumps(summary, indent=1))


if __name__ == "__main__":
    main()
