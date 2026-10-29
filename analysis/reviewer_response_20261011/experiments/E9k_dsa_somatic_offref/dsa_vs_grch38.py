"""E9k part 2 (user 2026-10-29): which COLO829BL-DSA somatic SNVs (somatic_snvs.vcf.gz) have a GRCh38 position, which exist
only on the DSA, and how they relate to the SMaHT GRCh38 truth and to Pansoma's graph-coordinate calls (match_dsa_somatic.py).

Lift: qfu's DSA -> hg38 chain `/scratch/qfu/COLO829BL_DSA/COLO829BL_hg38.chain` (chain/net; t = DSA contig, q = hg38),
parsed here (no liftOver binary on the cluster). Per SNV (DSA contig, position, REF, ALT):
  grch38_same      the position lifts to GRCh38 and GRCh38 has the DSA's REF base there (complemented on '-' chains)
  grch38_diff_ref  it lifts, but GRCh38 has another base: the normal's haplotype differs from GRCh38 at that base, so the
                   somatic SNV sits on a germline allele (on the graph: a non-GRCh38 node, if the allele is a graph path)
  dsa_only         the position is in no chain block (DSA sequence the chain does not align to GRCh38)
For grch38_same / grch38_diff_ref the lifted allele (GRCh38 position, DSA REF -> GRCh38 orientation, ALT likewise) is looked
up in the SMaHT GRCh38 truth (pansoma_v2_tensors/COLO829T_truth/COLO829T_somatic_snv_indel.vcf.gz): same position and ALT.
Writes dsa_vs_grch38.md / .json and dsa_snvs.tsv.gz here.
usage: /wanglab/jshen/anaconda3/bin/python -I dsa_vs_grch38.py
"""
import bisect
import csv
import gzip
import json
from collections import Counter, defaultdict
from pathlib import Path

import pysam

HERE = Path(__file__).resolve().parent
CHAIN = "/scratch/qfu/COLO829BL_DSA/COLO829BL_hg38.chain"
DSA_VCF = "/wanglab/jshen/scratch_jshen/data/COLO829T/somatic_snvs.vcf.gz"
TRUTH = "/scratch/jshen/data/pansoma_v2_tensors/COLO829T_truth/COLO829T_somatic_snv_indel.vcf.gz"
FASTA = "/scratch/jshen/data/HapMap/GCA_000001405.15_GRCh38_no_alt_analysis_set.fasta"
COMP = str.maketrans("ACGTN", "TGCAN")
AUTOSOMES = {f"chr{i}" for i in range(1, 23)}


def load_chain():
    """{DSA contig: sorted [(t_start, t_end, q_name, q_start, q_strand, q_size)]} of ungapped blocks (0-based, half-open)."""
    blocks = defaultdict(list)
    with open(CHAIN) as fh:
        head = None
        for line in fh:
            f = line.split()
            if not f:
                continue
            if f[0] == "chain":
                head = dict(t=f[2], tpos=int(f[5]), q=f[7], qsize=int(f[8]), qstrand=f[9], qpos=int(f[10]))
                continue
            size = int(f[0])
            blocks[head["t"]].append((head["tpos"], head["tpos"] + size, head["q"], head["qpos"], head["qstrand"], head["qsize"]))
            if len(f) == 3:
                head["tpos"] += size + int(f[1])
                head["qpos"] += size + int(f[2])
    for v in blocks.values():
        v.sort()
    return blocks


def lift(blocks, starts, contig, pos1):
    """(hg38 chrom, 1-based pos, strand) of DSA contig:pos1, or None."""
    v = blocks.get(contig)
    if not v:
        return None
    i = bisect.bisect_right(starts[contig], pos1 - 1) - 1
    if i < 0:
        return None
    t0, t1, q, q0, strand, qsize = v[i]
    if not (t0 <= pos1 - 1 < t1):
        return None
    off = pos1 - 1 - t0
    qp = q0 + off  # 0-based on the q strand
    if strand == "-":
        qp = qsize - 1 - qp
    return q, qp + 1, strand


def main():
    blocks = load_chain()
    starts = {k: [b[0] for b in v] for k, v in blocks.items()}
    fa = pysam.FastaFile(FASTA)
    truth = set()
    with pysam.VariantFile(TRUTH) as vcf:
        for r in vcf:
            if r.chrom in AUTOSOMES and len(r.ref) == 1 and r.alts and len(r.alts[0]) == 1:
                truth.add((r.chrom, r.pos, r.alts[0].upper()))
    rows = []
    with pysam.VariantFile(DSA_VCF) as vcf:
        for r in vcf:
            ref, alt = r.ref.upper(), r.alts[0].upper()
            filt = ";".join(r.filter.keys()) or "."
            L = lift(blocks, starts, r.chrom, r.pos)
            row = dict(contig=r.chrom, pos=r.pos, ref=ref, alt=alt, filter=filt, grch38="", grch38_ref="", cls="dsa_only",
                       in_grch38_truth=False)
            if L:
                q, qp, strand = L
                gref = fa.fetch(q, qp - 1, qp).upper() if q in fa.references else "N"
                ref_g = ref if strand == "+" else ref.translate(COMP)
                alt_g = alt if strand == "+" else alt.translate(COMP)
                row.update(grch38=f"{q}:{qp}:{strand}", grch38_ref=gref, cls="grch38_same" if gref == ref_g else "grch38_diff_ref",
                           in_grch38_truth=(q, qp, alt_g) in truth, grch38_chrom=q)
            rows.append(row)
    cols = ["contig", "pos", "ref", "alt", "filter", "cls", "grch38", "grch38_ref", "in_grch38_truth"]
    with gzip.open(HERE / "dsa_snvs.tsv.gz", "wt") as w:
        w.write("\t".join(cols) + "\n")
        for r in rows:
            w.write("\t".join(str(r[c]) for c in cols) + "\n")
    n = len(rows)
    by = Counter(r["cls"] for r in rows)
    in_truth = Counter(r["cls"] for r in rows if r["in_grch38_truth"])
    auto = [r for r in rows if r.get("grch38_chrom", "") in AUTOSOMES]
    # the Pansoma graph-coordinate calls matched in part 1
    off = defaultdict(list)
    with gzip.open(HERE / "per_call.tsv.gz", "rt") as fh:
        for c in csv.DictReader(fh, delimiter="\t"):
            if c["group"] in ("unplaced_pass", "norescue_pass") and c["kind"] == "dsa_somatic":
                contig, pos, change, _ = c["dsa"].split(":")
                tag = f"{c['graph'].split()[0]} {c['platform']}{'*' if c['rescue_off_pass'] == 'True' else ''}"
                if tag not in off[(contig, int(pos))]:
                    off[(contig, int(pos))].append(tag)
    idx = {(r["contig"], r["pos"]): r for r in rows}
    L = ["# E9k part 2: the COLO829BL-DSA somatic SNVs on GRCh38 (dsa_vs_grch38.py)", "",
         f"DSA somatic SNVs: {n:,} ({DSA_VCF}); lift: {CHAIN}; GRCh38 truth: SMaHT COLO829T ({len(truth):,} autosomal SNVs).", "",
         "| class | DSA somatic SNVs | in the SMaHT GRCh38 truth |", "|---|---:|---:|"]
    lab = {"grch38_same": "lifts to GRCh38, same reference base", "grch38_diff_ref": "lifts, GRCh38 base differs (on a germline allele)",
           "dsa_only": "no GRCh38 position in the chain (DSA only)"}
    for k in ("grch38_same", "grch38_diff_ref", "dsa_only"):
        L.append(f"| {lab[k]} | {by[k]:,} ({100 * by[k] / n:.1f}%) | {in_truth[k]:,} |")
    L.append(f"| all | {n:,} | {sum(in_truth.values()):,} |")
    L += ["", f"Lifted to an autosome: {len(auto):,}; the SMaHT GRCh38 truth has {len(truth):,} autosomal SNVs, "
          f"{sum(in_truth.values()):,} of which ({100 * sum(in_truth.values()) / len(truth):.1f}%) are DSA somatic SNVs by this lift.", "",
          f"## The {len(off)} DSA somatic sites that Pansoma calls in graph coordinates (part 1; * = also PASS without the rescue)", "",
          "| DSA site | change | class | GRCh38 | in the SMaHT GRCh38 truth | calls |", "|---|---|---|---|---|---|"]
    oc = Counter()
    for (contig, pos), v in sorted(off.items()):
        r = idx[(contig, pos)]
        oc[(r["cls"], r["in_grch38_truth"])] += 1
        L.append(f"| {contig}:{pos} | {r['ref']}>{r['alt']} | {lab[r['cls']]} | {r['grch38'] or '-'} | {r['in_grch38_truth']} | {', '.join(v)} |")
    L += ["", "Summary of these sites: " + "; ".join(f"{lab[k]}, {'in' if t else 'not in'} the GRCh38 truth: {c}" for (k, t), c in sorted(oc.items()))]
    (HERE / "dsa_vs_grch38.md").write_text("\n".join(L) + "\n")
    (HERE / "dsa_vs_grch38.json").write_text(json.dumps(dict(n=n, classes=by, in_truth=in_truth, truth_autosomal_snvs=len(truth),
                                                               offref_sites={f"{k[0]}|{k[1]}": c for k, c in oc.items()}), indent=1) + "\n")
    print("\n".join(L))


if __name__ == "__main__":
    main()
