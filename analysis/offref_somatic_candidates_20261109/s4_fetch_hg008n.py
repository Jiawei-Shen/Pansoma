"""s4c: HG008-N reads over the HG008T s4a regions, fetched from GIAB in one connection (parallel remote readers were
refused in E9i), with retries and resume, into DATA/s4/HG008T/normal_<LABEL>.bam (sorted, indexed). Adapted from
tmp/e9i_hg008t_20261104/fetch_normal.py (unchanged there). The job must export CURL_CA_BUNDLE.

usage: s4_fetch_hg008n.py LABEL      (HiFi | Illumina)

Regions (regions.bed) are fetched in chunks of CHUNK; each finished chunk is written as parts/<LABEL>.<k>.bam (via a .tmp
name), so a rerun skips finished chunks. Within a chunk a record fetched by two adjacent regions is written once (key:
name, flag, start); at the end the parts are concatenated, sorted (SLURM_CPUS_PER_TASK threads) and indexed. A long read
fetched in two chunks stays twice in the BAM; s4_count counts reads by name, so it is not counted twice.
"""
import glob
import os
import sys
import time
from pathlib import Path

import pysam

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from common import DATA  # noqa: E402

GIAB = 'https://ftp-trace.ncbi.nlm.nih.gov/ReferenceSamples/giab/data_somatic/HG008/Liss_lab'
SRC = {
    'HiFi': (f'{GIAB}/BCM_Revio_20240313/HG008-N-D_PacBio-HiFi-Revio_20240313_68x_GRCh38-GIABv3.bam',
             '/scratch/jshen/Github/Pansoma/tmp/e9i_hg008t_20261104/idx/'
             'HG008-N-D_PacBio-HiFi-Revio_20240313_68x_GRCh38-GIABv3.bam.bai'),
    'Illumina': (f'{GIAB}/BCM_Illumina-WGS_20240313/HG008-N-D_Illumina_169x_GRCh38-GIABv3.bam',
                 f'{DATA}/s4/HG008T/idx/HG008-N-D_Illumina_169x_GRCh38-GIABv3.bam.bai'),
}
CHUNK = 1000
TRIES = 8


def main():
    label = sys.argv[1]
    url, bai = SRC[label]
    d = Path(f'{DATA}/s4/HG008T')
    parts = d / 'parts'
    parts.mkdir(exist_ok=True)
    regions = [line.split() for line in open(d / 'regions.bed')]
    chunks = [regions[i:i + CHUNK] for i in range(0, len(regions), CHUNK)]

    def open_remote():
        return pysam.AlignmentFile(url, index_filename=bai)

    bam = open_remote()
    header = bam.header
    t0 = time.time()
    n_total = 0
    for k, chunk in enumerate(chunks):
        part = parts / f'{label}.{k:05d}.bam'
        if part.exists():
            continue
        tmp = parts / f'{label}.{k:05d}.bam.tmp'
        out = pysam.AlignmentFile(str(tmp), 'wb', header=header)
        seen, n = set(), 0
        for c, a, b in chunk:
            for attempt in range(TRIES):
                try:
                    if bam is None:
                        bam = open_remote()
                    recs = list(bam.fetch(c, int(a), int(b)))
                    break
                except (OSError, ValueError) as e:
                    print('retry', k, c, a, b, attempt, e, file=sys.stderr, flush=True)
                    try:
                        if bam is not None:
                            bam.close()
                    except Exception:
                        pass
                    bam = None
                    time.sleep(5 * (attempt + 1))
            else:
                out.close()
                raise SystemExit(f'failed region {c}:{a}-{b} after {TRIES} tries')
            for r in recs:
                key = (r.query_name, r.flag, r.reference_start)
                if key not in seen:
                    seen.add(key)
                    out.write(r)
                    n += 1
        out.close()
        os.replace(tmp, part)
        n_total += n
        print(f'chunk {k + 1}/{len(chunks)}', n, 'records', f'{time.time() - t0:.0f}s', flush=True)
    files = sorted(glob.glob(str(parts / f'{label}.*.bam')))
    if len(files) != len(chunks):
        raise SystemExit(f'{len(files)} parts for {len(chunks)} chunks')
    merged = d / f'normal_{label}.merged.bam'
    pysam.cat('-o', str(merged), *files)
    final = d / f'normal_{label}.bam'
    threads = os.environ.get('SLURM_CPUS_PER_TASK', '4')
    pysam.sort('-@', threads, '-m', '1G', '-o', str(final) + '.tmp.bam', str(merged))
    os.replace(str(final) + '.tmp.bam', final)
    pysam.index('-@', threads, str(final))
    merged.unlink()
    print('done', len(regions), 'regions', len(chunks), 'chunks', f'{time.time() - t0:.0f}s')


if __name__ == '__main__':
    main()
