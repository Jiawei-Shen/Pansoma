"""E17c check: the scan's features reproduce scripts/filter_panel_of_normals.py itself.

    python e17c_check_script.py chr22 [--n 4000]     -> DATA/check_script_<chrom>.json (+ printed summary)

A random sample of the chromosome's labelled-tensor alleles (as e17c_scan.py builds them, half of them at a
position some PoN has a record at) is written to a VCF and tagged by the unmodified script (subprocess, its CLI);
its PANSOMA_PON field is compared, per PoN, with what the current rule gives from the scan features:
  SNV    gnomAD g_lvl >= 1e-4; dbSNP d_nonsom; 1000G k_any; CoLoRSdb c_lvl >= 1e-4
  INDEL  gnomAD / CoLoRSdb level >= 0.01 if hp >= 7 else >= 0.001
"""
from __future__ import annotations

import argparse
import json
import random
import subprocess
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
import e17c_scan as scan  # noqa: E402
from filter_panel_of_normals import PON_NAMES  # noqa: E402

import pysam  # noqa: E402

FI = {f: i for i, f in enumerate(scan.FLOORS)}


def expected(kind, f):
    g_lvl, c_lvl, d_nonsom, k_any, hp = f[1], f[3], f[5], f[6], f[7]
    if kind == "SNV":
        hits = (g_lvl >= FI[1e-4], bool(d_nonsom), bool(k_any), c_lvl >= FI[1e-4])
    else:
        floor = FI[1e-2] if hp >= 7 else FI[1e-3]
        hits = (g_lvl >= floor, False, False, c_lvl >= floor)
    return tuple(n for n, h in zip(PON_NAMES, hits) if h)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("chrom")
    ap.add_argument("--n", type=int, default=4000)
    args = ap.parse_args()
    chrom = args.chrom
    m = json.loads((scan.TENSORS / scan.D9[0] / "tensors" / "SNV" / "labels.manifest.json").read_text())
    fasta_path = m["fasta"]
    genome, repeats, fasta = scan.Genome(fasta_path), scan.RepeatContext(fasta_path), scan.MemFasta(fasta_path, chrom)
    genome.fasta = repeats.genome = fasta
    keys = set()
    for s in scan.D9:
        for kind in ("SNV", "INDEL"):
            _, counts = scan.read_labels((s, kind, chrom))
            for (_, _, lin) in counts:
                if lin is None or lin[0] != chrom or genome.get(lin[0], lin[1], lin[1] + len(lin[2])) != lin[2]:
                    continue
                c, pos0, ref, alt, ev = lin
                if ev in ("INS", "DEL"):
                    pos0, ref, alt, _ = genome.left_align(c, pos0, ref, alt, ev)
                keys.add((c, *genome.vcf_allele(c, pos0, ref, alt, ev)))
    positions = {k[1] for k in keys}
    pons, contigs = [None] * 4, [None] * 4
    for p in range(4):
        _, contig, recs, _ = scan.stream_pon((p, chrom, positions))
        pons[p], contigs[p] = scan.MemPoN(scan.PON_FILES[p], recs), contig
    rng = random.Random(f"E17c-script:{chrom}")
    with_rec = sorted(k for k in keys if any(k[1] in pon.at for pon in pons))
    sample = set(rng.sample(sorted(keys), min(args.n // 2, len(keys))))
    sample |= set(rng.sample(with_rec, min(args.n // 2, len(with_rec))))
    sample = sorted(sample, key=lambda k: (k[1], k[2], k[3]))

    out_dir = scan.DATA / "check_script"
    out_dir.mkdir(parents=True, exist_ok=True)
    vcf_in, vcf_out = out_dir / f"{chrom}.alleles.vcf", out_dir / f"{chrom}.alleles.pon.vcf.gz"
    with pysam.FastaFile(fasta_path) as fa:
        contig_lines = [f"##contig=<ID={n},length={l}>" for n, l in zip(fa.references, fa.lengths)]
    with open(vcf_in, "w") as fh:
        fh.write("##fileformat=VCFv4.2\n" + "\n".join(contig_lines) + f"\n##reference=file://{fasta_path}\n")
        fh.write("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n")
        for i, (c, pos, ref, alt) in enumerate(sample):
            fh.write(f"{c}\t{pos}\ta{i}\t{ref}\t{alt}\t.\tPASS\t.\n")
    vcf_gz = pysam.tabix_index(str(vcf_in), preset="vcf", force=True)  # bgzip + index
    pon_paths = [str(scan.PON_DIR / f) for f in scan.PON_FILES]
    run = subprocess.run([sys.executable, str(scan.REPO / "scripts" / "filter_panel_of_normals.py"), vcf_gz,
                          str(vcf_out), "--pon", *pon_paths, "--fasta", fasta_path],
                         check=True, capture_output=True, text=True)
    tags = {}
    with pysam.VariantFile(str(vcf_out)) as v:
        for rec in v:
            tags[rec.id] = tuple(rec.info.get("PANSOMA_PON") or ())
    diffs, n_tagged, by_kind = [], 0, {"SNV": 0, "INDEL": 0}
    for i, key in enumerate(sample):
        kind = "SNV" if len(key[2]) == len(key[3]) else "INDEL"
        by_kind[kind] += 1
        f = scan.features(pons, contigs, key, fasta, repeats)
        exp = expected(kind, f)
        got = tags.get(f"a{i}", ())
        n_tagged += bool(got)
        if exp != got:
            diffs.append(dict(key=list(key), expected=exp, script=got, features=list(f)))
    result = dict(chrom=chrom, alleles=len(sample), by_kind=by_kind, script_tagged=n_tagged, differences=len(diffs),
                  examples=diffs[:20], script_stdout=run.stdout.strip().splitlines())
    (scan.DATA / f"check_script_{chrom}.json").write_text(json.dumps(result, indent=1))
    print(json.dumps({k: v for k, v in result.items() if k != "examples"}, indent=1))


if __name__ == "__main__":
    main()
