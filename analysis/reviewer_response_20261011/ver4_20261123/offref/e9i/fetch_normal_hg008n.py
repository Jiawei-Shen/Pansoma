"""Fetch the HG008-N-D PacBio HiFi reads (GIAB BCM Revio 68x, GRCh38-GIABv3) over regions.bed (HG008T default-run call
loci +- 600 bp, merged) from GIAB in ONE connection at a time (18 parallel remote readers were refused by the server),
with retries, into normal_subset.bam (sorted, indexed). Reads spanning two regions are written once."""
import time, sys
import pysam
from pathlib import Path
import os
HERE = Path(os.environ.get("WORK", "/scratch/jshen/data/reviewer_experiments_20261013/ver4_20261123/offref/e9i/HG008T"))   # ver4 offref: data dir (regions.bed, idx/, normal_subset.bam)
URL = ("https://ftp-trace.ncbi.nlm.nih.gov/ReferenceSamples/giab/data_somatic/HG008/Liss_lab/BCM_Revio_20240313/"
       "HG008-N-D_PacBio-HiFi-Revio_20240313_68x_GRCh38-GIABv3.bam")
BAI = HERE / "idx/HG008-N-D_PacBio-HiFi-Revio_20240313_68x_GRCh38-GIABv3.bam.bai"
import os
REG = os.environ.get("REGIONS", "regions.bed"); OUT = os.environ.get("OUT_PREFIX", "normal_subset")
regions = [l.split() for l in open(HERE / REG)]
def open_remote():
    return pysam.AlignmentFile(URL, index_filename=str(BAI))
bam = open_remote()
out = pysam.AlignmentFile(str(HERE / f"{OUT}.unsorted.bam"), "wb", template=bam)
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
            seen.add(k); out.write(r); n += 1
    if i % 100 == 0:
        print(i, "regions", n, "reads", flush=True)
out.close()
pysam.sort("-o", str(HERE / f"{OUT}.bam"), str(HERE / f"{OUT}.unsorted.bam"))
pysam.index(str(HERE / f"{OUT}.bam"))
Path(HERE / f"{OUT}.unsorted.bam").unlink()
print("done", len(regions), "regions", n, "reads")
