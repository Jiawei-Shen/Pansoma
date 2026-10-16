#!/bin/bash
# Collect the E8b Slurm outputs (jobs in /scratch/jshen/Github/Pansoma/tmp/reviewer_cpu_20261013/E8b_literature_facts/jobs.txt)
# into the small result tables of this folder, with percentages. Safe to re-run; prints which inputs are missing.
set -u
E=/scratch/jshen/Github/Pansoma/analysis/reviewer_response_20261011/experiments/E8b_literature_facts
O=/scratch/jshen/data/reviewer_experiments_20261013/E8b_literature_facts
SETS="HG008T_GIAB_v0.2 COLO829T_SMaHT HCC1395_SEQC2_v1.2.1 HCC1954_CASTLE H1437_CASTLE H2009_CASTLE HapMap_positives_chr20"
out=$E/somatic_truth_in_resources.tsv
first=1
for s in $SETS; do
  f=$O/somatic_truth_in_resources.$s.tsv
  if [ ! -s "$f" ]; then echo "missing: $f"; continue; fi
  if [ $first = 1 ]; then head -1 "$f" > $out; first=0; fi
  tail -n +2 "$f" >> $out
done
[ $first = 0 ] && /wanglab/jshen/anaconda3/bin/python - "$out" "$E/somatic_truth_in_resources.pct.tsv" <<'PY'
import csv, sys
rows = list(csv.DictReader(open(sys.argv[1]), delimiter="\t"))
cols = [c for c in rows[0] if c not in ("truth_set", "kind", "n")]
with open(sys.argv[2], "w") as fh:
    fh.write("truth_set\tkind\tn\t" + "\t".join(f"{c}_pct" for c in cols) + "\n")
    for r in rows:
        n = int(r["n"])
        fh.write(f"{r['truth_set']}\t{r['kind']}\t{n}\t" +
                 "\t".join(f"{100 * int(r[c]) / n:.2f}" if n else "NA" for c in cols) + "\n")
PY
m=$O/membership_private_alleles.chr20.tsv
if [ -s $m ]; then cp $m $E/membership_private_alleles.chr20.tsv; else echo "missing: $m"; fi
ls -la $E/*.tsv
