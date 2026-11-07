"""E20b step 2a: the HG008-N normal reads at the loci. No HG008-N read data is local, so the HG008-N-D PacBio HiFi BAM
(GIAB, BCM Revio 68x, GRCh38-GIABv3, pbmm2) is read remotely over regions.bed (loci.py) with one connection at a time
and retries, as tmp/e9i_hg008t_20261104/fetch_normal.py does for R3.3 (18 parallel remote readers were refused by the
server). Every record overlapping a region is written once; extract.sh then keeps the primary records.
Output: $D/normal_subset.bam (sorted, indexed). Index: a local copy of the GIAB .bai in $D/idx.
usage: /wanglab/jshen/anaconda3/bin/python fetch_normal.py
"""
import sys
import time
from pathlib import Path

import pysam

D = Path("/scratch/jshen/data/reviewer_experiments_20261013/E20b_normal_reads_hg008n")
URL = ("https://ftp-trace.ncbi.nlm.nih.gov/ReferenceSamples/giab/data_somatic/HG008/Liss_lab/BCM_Revio_20240313/"
       "HG008-N-D_PacBio-HiFi-Revio_20240313_68x_GRCh38-GIABv3.bam")
BAI = D / "idx/HG008-N-D_PacBio-HiFi-Revio_20240313_68x_GRCh38-GIABv3.bam.bai"
OUT = D / "normal_subset"


def open_remote():
    return pysam.AlignmentFile(URL, index_filename=str(BAI))


def main():
    regions = [line.split() for line in open(D / "regions.bed")]
    bam = open_remote()
    tmp = Path(f"{OUT}.unsorted.bam")
    out = pysam.AlignmentFile(str(tmp), "wb", template=bam)
    seen, n = set(), 0
    for i, (c, a, b) in enumerate(regions):
        for attempt in range(8):
            try:
                recs = list(bam.fetch(c, int(a), int(b)))
                break
            except (OSError, ValueError) as e:
                print("retry", i, c, a, b, attempt, e, file=sys.stderr, flush=True)
                time.sleep(5 * (attempt + 1))
                try:
                    bam.close()
                except Exception:
                    pass
                bam = open_remote()
        else:
            raise SystemExit(f"failed region {c}:{a}-{b}")
        for r in recs:
            k = (r.query_name, r.flag, r.reference_start)
            if k not in seen:
                seen.add(k)
                out.write(r)
                n += 1
        if i % 100 == 0:
            print(i, "regions", n, "records", flush=True)
    out.close()
    pysam.sort("-o", f"{OUT}.bam", str(tmp))
    pysam.index(f"{OUT}.bam")
    tmp.unlink()
    print("done", len(regions), "regions", n, "records")


if __name__ == "__main__":
    main()
