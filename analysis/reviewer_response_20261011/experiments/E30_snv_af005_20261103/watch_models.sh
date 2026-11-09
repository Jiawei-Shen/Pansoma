#!/bin/bash
# Exits (one notification) when a training's four vcf jobs have all ended and that run is not in reported.txt yet, or
# when any job of the Illumina / ONT / multicancer4 chains (470391-470506, 470902-470912) ends badly (reported_bad.txt: known). Re-run after each exit.
W=/scratch/jshen/Github/Pansoma/tmp/snv_af05_20261103; S=$W/reported.txt; touch $S
J=$(seq -s, 470391 470506),$(seq -s, 470902 470912),$(seq -s, 470962 470966),$(seq -s, 471037 471066)
while true; do
  st=$(sacct -n -X -P -j $J -o JobID,JobName,State,Elapsed 2>/dev/null)
  bad=$(echo "$st" | grep -E "\|(FAILED|CANCELLED|TIMEOUT|OUT_OF_ME|NODE_FAIL)" | grep -v "_fallback|CANCELLED" | grep -vxFf $W/reported_bad.txt 2>/dev/null)
  never=$(squeue -h -j $J -o "%i %j %r" 2>/dev/null | grep -i "DependencyNeverSatisfied")
  if [ -n "$bad$never" ]; then echo "PROBLEM"; echo "$bad"; echo "$never"; echo "$bad" >> $W/reported_bad.txt; exit 1; fi
  for run in $(echo "$st" | grep -o "|v_multicancer[^|]*" | sed -E 's/^\|v_//; s/_(HCC1395T|HG008T|COLO829T)_.*//' | sort -u); do
    grep -qx "$run" $S && continue
    v=$(echo "$st" | grep "|v_${run}_" | grep -v "|CANCELLED")   # a cancelled (superseded) chain does not count
    [ "$(echo "$v" | grep -c "|COMPLETED|")" = 4 ] && { echo "MODEL DONE $run"; echo "$v"; echo "$run" >> $S; exit 0; }
  done
  sleep 300
done
