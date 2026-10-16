#!/bin/bash
# E22 check 2: which vg versions write a PARAMS_JSON group into a giraffe GAM, and with which payload?
# Builds a 300-kb linear test graph (the first N-free 300-kb window of GRCh38 chr20 from 40 Mb) with each vg version, maps simulated
# paired short reads (default preset) and long reads (-b hifi, -b r10), and scans the GAM heads with gam_tag_scan.py.
# The payloads are compared with the PARAMS_JSON of the COLO829T / HG008T d9 GAMs (gam_params/, written by
# gam_params.sh). Light: 2 threads, a few minutes on the login node. Never touches a real GAM.
set -uo pipefail
HERE=/scratch/jshen/Github/Pansoma/analysis/reviewer_response_20261011/experiments/E22_truth_and_version_checks
W=/scratch/jshen/data/reviewer_experiments_20261013/E22_truth_and_version_checks/vg_fingerprint
PY=/wanglab/jshen/anaconda3/bin/python
FASTA=/scratch/jshen/data/HapMap/GCA_000001405.15_GRCh38_no_alt_analysis_set.fasta
mkdir -p $W
cd $W
if [ ! -s tiny.fa ]; then
  $PY - "$FASTA" <<'PY'
import random, sys
import pysam
fa = pysam.FastaFile(sys.argv[1])
for start in range(40_000_000, 60_000_000, 1_000_000):  # first N-free 300-kb window of chr20 from 40 Mb
    seq = fa.fetch("chr20", start, start + 300_000).upper()
    if "N" not in seq:
        break
print("test region chr20", start, start + 300_000, file=sys.stderr)
open("tiny.fa", "w").write(">chr20sub\n" + "\n".join(seq[i:i + 80] for i in range(0, len(seq), 80)) + "\n")
rnd = random.Random(7)
comp = str.maketrans("ACGT", "TGCA")
def mutate(s, rate):
    return "".join(rnd.choice("ACGT".replace(c, "")) if rnd.random() < rate else c for c in s)
with open("r1.fq", "w") as f1, open("r2.fq", "w") as f2:
    for i in range(2000):
        ins = int(rnd.gauss(400, 40))
        p = rnd.randrange(0, len(seq) - ins)
        frag = seq[p:p + ins]
        a, b = mutate(frag[:150], 0.002), mutate(frag[-150:].translate(comp)[::-1], 0.002)
        f1.write(f"@p{i}/1\n{a}\n+\n{'I' * 150}\n")
        f2.write(f"@p{i}/2\n{b}\n+\n{'I' * 150}\n")
with open("long.fq", "w") as f:
    for i in range(30):
        n = rnd.randrange(5000, 12000)
        p = rnd.randrange(0, len(seq) - n)
        s = mutate(seq[p:p + n], 0.005)
        if i % 2:
            s = s.translate(comp)[::-1]
        f.write(f"@l{i}\n{s}\n+\n{'I' * len(s)}\n")
PY
fi
VERSIONS="${VERSIONS:-/scratch/wzhang/bin/vg_v1.50.0 /scratch/wzhang/bin/vg_v1.53.0 /scratch/wzhang/bin/vg_v1.55.0 /scratch/wzhang/bin/vg_v1.56.0
/scratch/wzhang/bin/vg_v1.57.0 /scratch/wzhang/bin/vg_v1.58.0 /opt/apps/vg/1.59.0/vg /scratch/wzhang/bin/vg_v1.60.0 /scratch/wzhang/bin/vg_v1.61.0 /scratch/wzhang/bin/vg_v1.62.0
/scratch/wzhang/bin/vg_v1.63.0 /scratch/wzhang/bin/vg_v1.63.1 /scratch/wzhang/bin/vg_v1.64.0 /scratch/wzhang/bin/vg_v1.64.1
/scratch/wzhang/bin/vg_v1.65.0 /scratch/wzhang/bin/vg_v1.66.0 /scratch/jshen/bin/vg_v1.77.0}"
for VG in $VERSIONS; do
  v=$($VG version 2>/dev/null | head -1 | awk '{print $3}')
  d=$W/$v
  mkdir -p $d
  echo "=== $v ($VG)"
  ( cd $d
    [ -s tiny.giraffe.gbz ] || $VG autoindex --workflow giraffe -r $W/tiny.fa -p tiny -t 2 > autoindex_sr.log 2>&1
    if [ -s tiny.shortread.zipcodes ]; then  # vg >= 1.63 autoindex: minimizers with zipcodes
      $VG giraffe -Z tiny.giraffe.gbz -d tiny.dist -m tiny.shortread.withzip.min -z tiny.shortread.zipcodes \
        -f $W/r1.fq -f $W/r2.fq -t 2 > short.gam 2> short.err
      echo "short exit $?"
      # the May-2025 d9 Illumina run could only use the 2023 HPRC .min (no zipcodes): an old-format .min, from 1.62
      $VG giraffe -Z tiny.giraffe.gbz -d tiny.dist -m $W/v1.62.0/tiny.min \
        -f $W/r1.fq -f $W/r2.fq -t 2 > short_oldmin.gam 2> short_oldmin.err
      echo "short_oldmin exit $?"
    else
      $VG giraffe -Z tiny.giraffe.gbz -d tiny.dist -m tiny.min -f $W/r1.fq -f $W/r2.fq -t 2 > short.gam 2> short.err
      echo "short exit $?"
    fi
    [ -s tiny.longread.zipcodes ] || $VG autoindex --workflow lr-giraffe -r $W/tiny.fa -p tiny -t 2 > autoindex_lr.log 2>&1
    if [ -s tiny.longread.zipcodes ]; then
      for b in hifi r10; do
        $VG giraffe -b $b -Z tiny.giraffe.gbz -d tiny.dist -m tiny.longread.withzip.min -z tiny.longread.zipcodes \
          -f $W/long.fq -t 2 > $b.gam 2> $b.err
        echo "$b exit $?"
      done
    else
      echo "no lr-giraffe indexes"
    fi )
done
