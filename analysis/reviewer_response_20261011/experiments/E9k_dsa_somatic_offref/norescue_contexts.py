"""E9k part 3: contexts of the unplaced PASS calls made WITHOUT the off-reference rescue (the default since 91a8261), which
E9c did not process (E3b's `<model>__norescue/vcf_chr1` re-calls), so that match_dsa_somatic.py can place them on the DSA.

  norescue_contexts.py run VCFDIR NAME   E9c's e9c_run.py on VCFDIR with its DATA redirected to E9k's data folder (no
                                         controls; nothing written to E9c's folder) -> DATA/runs/NAME/calls.tsv
  norescue_contexts.py queries           every alt_read context of DATA/runs/*/calls.tsv -> DATA/asm/queries.fa (qid as E9d)
Job script: tmp/reviewer_cpu_20261029/E9k/norescue.sh (minimap2 -cx sr --cs --secondary=yes -N 50 to DSA hapX / hapY, as E9d).
"""
import csv
import hashlib
import os
import sys
from pathlib import Path

E9C = Path(__file__).resolve().parents[1] / "E9c_offref_filter_validation"
DATA = "/scratch/jshen/data/reviewer_experiments_20261013/E9k_dsa_somatic_offref"


def qid(seq, site):  # as E9d e9d_asm.qid
    return "q" + hashlib.sha1(f"{seq}|{site}".encode()).hexdigest()[:16]


def run(vcfdir, name):
    sys.path.insert(0, str(E9C))
    import e9c_lib as L
    L.DATA = DATA
    import e9c_run
    sys.argv = ["e9c_run.py", vcfdir, name, "--controls", "0"]
    e9c_run.main()


def queries():
    csv.field_size_limit(sys.maxsize)
    seen = {}
    for p in sorted(Path(DATA, "runs").glob("*/calls.tsv")):
        for r in csv.DictReader(open(p), delimiter="\t"):
            s, site = r.get("alt_read", ""), r.get("alt_read_site", "")
            if s and site and len(s) >= 25:
                seen.setdefault(qid(s, site), s)
    os.makedirs(f"{DATA}/asm", exist_ok=True)
    with open(f"{DATA}/asm/queries.fa", "w") as fa:
        for k, s in sorted(seen.items()):
            fa.write(f">{k}\n{s}\n")
    print(f"{len(seen):,} alt_read contexts -> {DATA}/asm/queries.fa")


if __name__ == "__main__":
    if sys.argv[1] == "run":
        run(sys.argv[2], sys.argv[3])
    else:
        queries()
