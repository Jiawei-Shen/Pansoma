"""Why do d46's own false positives pass the PoN? Population AF of each PASS FP (its site features) and examples of the
"own only: PoN removed" group (the GRCh38-only arm had a candidate there, removed by the PoN).

    ris/submit.sh -g 0 -n 2 -m 16 pnv2_fp_popaf py.sh $PWD/analysis/indel_ppr8_9919d4a_20261026/fp_context/fp_popaf.py

Groups as fp_context.py (seeds 0-2, chr1, after the PoN, in the BED; 10 bp window): shared / other cand. / own only:
PoN removed / own only: no candidate; GRCh38-only's FPs for contrast.
Per FP, the site features of its tensor in the run's set (<SET>.INDEL.npz, matched on the VCF allele vcf_pos / vcf_ref /
vcf_alt that linear_vcf writes): popaf_exact (largest gnomAD / CoLoRSdb AF of the same allele after normalisation),
popaf_pos (largest AF of any INDEL allele at the left-aligned position), pon_tag (the PoN script's exact-record match),
hp_len. The PoN's INDEL floor: AF >= 0.01 when hp_len >= 8, else >= 1e-4 (scripts/filter_panel_of_normals.py).
  same allele, missed   popaf_exact >= floor, pon_tag 0: the allele is in a database after normalisation, but no
                        record matches the call's POS / REF / ALT exactly (another representation)
  other allele here     popaf_exact < floor, popaf_pos >= floor: another INDEL allele at that position is common
  no population signal  both < floor
Examples ("own only: PoN removed", seed 0): the homopolymer tract on GRCh38, d46's call as a length change of the tract,
the GRCh38-only arm's PoN-tagged candidates within 10 bp and the gnomAD / CoLoRSdb alleles of the tract (same base,
INS / DEL of that base only) with their AF.
Writes popaf.md next to this script.
"""
import bisect
import glob
import gzip
import json
import os
import sys
from collections import Counter, defaultdict
from pathlib import Path

import numpy as np
import pysam

HERE = Path(__file__).resolve().parent
REPO = HERE.parents[2]
sys.path.insert(0, str(REPO))
from indexed_gam_pipeline_v4.tensor_postprocessing.truth_labels import split_alleles  # noqa: E402

RUNS = Path(os.environ.get("PANSOMA_RUNS", "/storage2/fs1/epigenome/Active/jiawei.s/pansoma_net_v2_runs"))
DATA = Path(os.environ.get("PANSOMA_DATA", "/storage2/fs1/epigenome/Active/jiawei.s/Dataset"))
SF = Path(os.environ.get("SITE_FEATURES", RUNS / "site_features"))
FASTA = pysam.FastaFile(str(DATA / "REF/GCA_000001405.15_GRCh38_no_alt_analysis_set.fa"))
DBS = {"gnomAD": pysam.VariantFile(str(DATA / "panel_of_normal_VCFs/af-only-gnomad.hg38.vcf.gz")),
       "CoLoRSdb": pysam.VariantFile(str(DATA / "panel_of_normal_VCFs/CoLoRSdb.GRCh38.v1.1.0.deepvariant.glnexus.vcf.gz"))}
SETS = (("COLO829", "Illumina"), ("COLO829", "fiberseq"), ("HG008", "Illumina"), ("HG008", "PacBio"))
SEEDS = ("", "_s1", "_s2")
W = 10
CATS = ("same allele, missed", "other allele here", "no population signal", "no tensor match")
N_EXAMPLES = 12


def records(path):
    out = []
    with gzip.open(path, "rt") as f:
        for line in f:
            if not line.startswith("#"):
                c = line.split("\t", 7)
                out.append((c[0], int(c[1]), c[3].upper(), c[4].split(",")[0].upper(), c[6]))
    return out


def near(pos, positions):
    i = bisect.bisect_left(positions, pos - W)
    return i < len(positions) and positions[i] <= pos + W


def features(set_name):
    """(vcf_pos, vcf_ref, vcf_alt) -> (popaf_exact, popaf_pos, pon_tag, hp_len), chr1 tensors, the largest per allele."""
    z = np.load(SF / f"{set_name}.INDEL.npz", allow_pickle=False)
    chroms, counts = [str(c) for c in z["chroms"]], z["chrom_counts"]
    i = chroms.index("chr1")
    lo, hi = int(counts[:i].sum()), int(counts[:i + 1].sum())
    a = {k: z[k][lo:hi] for k in ("vcf_pos", "vcf_ref", "vcf_alt", "popaf_exact", "popaf_pos", "pon_tag", "hp_len")}
    out = {}
    for k in np.flatnonzero(a["pon_tag"] >= 0):
        key = (int(a["vcf_pos"][k]), a["vcf_ref"][k].decode(), a["vcf_alt"][k].decode())
        v = (float(a["popaf_exact"][k]), float(a["popaf_pos"][k]), int(a["pon_tag"][k]), int(a["hp_len"][k]))
        out[key] = max(out.get(key, v), v)
    return out


def category(f):
    if f is None:
        return CATS[3]
    exact, pos, tag, hp = f
    floor = 0.01 if hp >= 8 else 1e-4
    return CATS[0] if exact >= floor and tag == 0 else CATS[1] if pos >= floor else CATS[2]


def event(pos, ref, alt):
    """(pos0, kind, inserted / deleted bases) of a pure INS / DEL, else None."""
    for pos0, r, a, kind, _ in split_alleles(pos, ref, [alt]):
        if kind in ("INS", "DEL"):
            return pos0, kind, a if kind == "INS" else r
    return None


def tract(pos0, base):
    """(start0, length) of the run of `base` on chr1 touching pos0 (an INS before pos0, a DEL from pos0)."""
    seq = FASTA.fetch("chr1", pos0 - 200, pos0 + 200).upper()
    i = 200
    a = i
    while a > 0 and seq[a - 1] == base:
        a -= 1
    b = i
    while b < len(seq) and seq[b] == base:
        b += 1
    return pos0 - 200 + a, b - a


def hp_delta(pos, ref, alt):
    """(base, tract start0, tract length, signed length change) when the call inserts / deletes copies of one base."""
    e = event(pos, ref, alt)
    if not e or len(set(e[2])) != 1:
        return None
    pos0, kind, bases = e
    start, length = tract(pos0, bases[0])
    return bases[0], start, length, len(bases) if kind == "INS" else -len(bases)


def population(base, start, length):
    """{signed change: (largest AF, databases)} of the database INS / DEL alleles of `base` touching the tract."""
    out = {}
    for name, vf in DBS.items():
        for rec in vf.fetch("chr1", max(0, start - 5), start + length + 5):
            afs = rec.info.get("AF")
            for k, alt in enumerate(rec.alts or ()):
                h = hp_delta(rec.pos, rec.ref.upper(), alt.upper()) if alt and not alt.startswith("<") else None
                if not h or h[0] != base or h[1] != start:
                    continue
                af = float(afs[k]) if afs is not None and k < len(afs) and afs[k] is not None else 0.0
                old = out.get(h[3], (0.0, set()))
                out[h[3]] = (max(old[0], af), old[1] | {name})
    return out


def main():
    L = ["# Why d46's false positives pass the PoN (population AF of each PASS FP)\n",
         "(`fp_popaf.py`; chr1, after the PoN, in the benchmark BED, seeds 0-2 summed; groups as `tables.md`. "
         "same allele, missed = the call's allele is in gnomAD / CoLoRSdb at AF >= the PoN floor after normalisation, "
         "but the PoN's exact POS / REF / ALT match does not find it; other allele here = another INDEL allele at the "
         "position is at AF >= the floor; floor 0.01 in homopolymers >= 8 bp, else 1e-4.)\n"]
    summary = {}
    for s, p in SETS:
        cnt = defaultdict(Counter)
        examples = []
        for x in SEEDS:
            run = {arm: RUNS / f"{s}_{p}_INDEL_{arm}_ppr8{x}" for arm in ("d46", "linear")}
            feat = {arm: features(Path(json.loads((run[arm] / "args.json").read_text())["tensors"][0]).parent.name)
                    for arm in run}
            fp, cand, raw = {}, {}, {}
            for arm, r in run.items():
                v = r / "vcf_chr1" / "eval_bed" / "vcfeval"
                fp[arm] = records(v / "pon_calls" / "fp.vcf.gz")
                cand[arm] = sorted(q[1] for q in records(v / "pon_curve_input.vcf.gz"))
                raw[arm] = sorted(q[1] for q in records(v / "raw_curve_input.vcf.gz"))
            lin_pon = [q for q in records(glob.glob(str(run["linear"] / "vcf_chr1/eval_bed/*.linear.pon.vcf.gz"))[0])
                       if "PanelOfNormals" in q[4]]
            lin_fp = sorted(q[1] for q in fp["linear"])
            d46_fp = sorted(q[1] for q in fp["d46"])
            for q in fp["d46"]:
                g = ("shared" if near(q[1], lin_fp) else "other cand." if near(q[1], cand["linear"])
                     else "own only: PoN removed" if near(q[1], raw["linear"]) else "own only: no candidate")
                f = feat["d46"].get(q[1:4])
                cnt[("d46", g)][category(f)] += 1
                if g == "own only: PoN removed" and x == "" and len(examples) < N_EXAMPLES:
                    examples.append((q, f, [r for r in lin_pon if abs(r[1] - q[1]) <= W]))
            for q in fp["linear"]:
                g = "own only" if not near(q[1], d46_fp) and not near(q[1], cand["d46"]) else "other"
                cnt[("linear", g)][category(feat["linear"].get(q[1:4]))] += 1
        L += [f"\n## {s}T {p}\n", "| FPs | n | " + " | ".join(CATS) + " |", "|---" * (len(CATS) + 2) + "|"]
        rows = [("d46", g) for g in ("shared", "other cand.", "own only: PoN removed", "own only: no candidate")] + \
               [("linear", "own only"), ("linear", "other")]
        for arm, g in rows:
            c = cnt[(arm, g)]
            n = sum(c.values())
            label = f"{'d46' if arm == 'd46' else 'GRCh38-only'} {g}"
            L.append(f"| {label} | {n} | " + " | ".join(f"{c[k]} ({100 * c[k] / n:.0f}%)" if n else "-" for k in CATS) + " |")
            summary[f"{s}_{p}/{arm}/{g}"] = dict(n=n, **c)
        L += [f"\nExamples, d46 own only: PoN removed (seed 0, first {N_EXAMPLES}): tract = GRCh38 homopolymer "
              "(base x length); changes are signed lengths of the tract; population = gnomAD / CoLoRSdb alleles of "
              "the tract (largest AF)\n",
              "| d46 call | tract | d46 change | popaf exact / pos | GRCh38-only PoN-tagged candidates | population alleles |",
              "|---|---|---|---|---|---|"]
        for q, f, lin in examples:
            h = hp_delta(q[1], q[2], q[3])
            pops = population(h[0], h[1], h[2]) if h else {}
            lin_txt = ", ".join(
                (f"{d[3]:+d}" if (d := hp_delta(r[1], r[2], r[3])) and h and d[0] == h[0] and d[1] == h[1]
                 else f"{r[1]} {r[2]}>{r[3]}") for r in lin) or "-"
            pop_txt = ", ".join(f"{k:+d}: {v[0]:.3g}" for k, v in sorted(pops.items())) or "-"
            call = f"{q[1]} {q[2][:12]}>{q[3][:12]}"
            L.append(f"| {call} | {f'{h[0]} x {h[2]}' if h else 'not a homopolymer'} | {f'{h[3]:+d}' if h else '-'} | "
                     f"{f'{f[0]:.2g} / {f[1]:.2g}' if f else 'no tensor'} | {lin_txt} | {pop_txt} |")
    (HERE / "popaf.md").write_text("\n".join(L) + "\n")
    (HERE / "popaf.json").write_text(json.dumps(summary, indent=1) + "\n")
    print("\n".join(L))


if __name__ == "__main__":
    main()
