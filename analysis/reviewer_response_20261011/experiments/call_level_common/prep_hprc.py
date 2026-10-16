"""Normalized chr1 alleles of the HPRC v1.1 haplotypes -> DATA/hprc_v1.1_chr1_alleles.norm.tsv.gz (+ .tbi).

Source: hprc-v1.1-mc-grch38.vcfbub.a100k.wave.vcf.gz (vg deconstruct of HPRC v1.1 MC GRCh38, bubbles > 100 kb popped
by vcfbub, complex alleles decomposed by vcfwave; 44 diploid samples + CHM13). Every ALT with INFO AC > 0 is
left-aligned and trimmed on the GRCh38 FASTA (calls_lib.normalize); an MNP also adds its single-base substitutions.
Rows: chrom, pos, ref, alt, AC (the largest AC of the records that give the key). Skipped when the table exists.
usage: python prep_hprc.py [FASTA]
"""
import sys
import time
from pathlib import Path

import pysam

sys.path.insert(0, str(Path(__file__).resolve().parent))
import calls_lib as L  # noqa: E402

FASTA = sys.argv[1] if len(sys.argv) > 1 else "/scratch/jshen/data/HapMap/GCA_000001405.15_GRCh38_no_alt_analysis_set.fasta"


def main():
    out = L.HPRC_TABLE
    if out.exists() and Path(str(out) + ".tbi").exists():
        L.log(f"{out} exists")
        return
    t0 = time.time()
    seq = L.chrom_seq(FASTA)
    keys = {}
    n = skipped = 0
    tb = pysam.TabixFile(str(L.HPRC_VCF))
    for line in tb.fetch(L.CHROM):
        c = line.split("\t", 8)
        n += 1
        ac = None
        for item in c[7].split(";"):
            if item.startswith("AC="):
                ac = item[3:].split(",")
        alts = c[4].split(",")
        for i, alt in enumerate(alts):
            a = int(ac[i]) if ac and i < len(ac) and ac[i] not in (".", "") else 0
            if a <= 0:
                skipped += 1
                continue
            for k in L.keys_of(seq, int(c[1]), c[3], alt):
                if keys.get(k, 0) < a:
                    keys[k] = a
        if n % 500000 == 0:
            L.log(f"{n} records, {len(keys)} keys, {time.time() - t0:.0f} s")
    L.log(f"{n} records, {skipped} ALTs with AC 0, {len(keys)} keys; writing")
    plain = out.with_suffix("")  # .tsv
    tmp = Path(str(plain) + ".tmp")
    with open(tmp, "w") as f:
        for (p, r, a), ac in sorted(keys.items()):
            f.write(f"{L.CHROM}\t{p}\t{r}\t{a}\t{ac}\n")
    pysam.tabix_compress(str(tmp), str(out), force=True)
    tmp.unlink()
    pysam.tabix_index(str(out), seq_col=0, start_col=1, end_col=1, zerobased=False, force=True)
    L.log(f"wrote {out} in {time.time() - t0:.0f} s")


if __name__ == "__main__":
    main()
