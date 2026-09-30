"""Pattern files for filtered_scan.sh: the candidate keys of every truth allele whose recall status is `filtered`.

One file per sample and kind directory (SNV keys -> SNV/filtered_candidates.ndjson, INS/DEL keys -> INDEL/...),
one line per key: {"candidate_id": "KEY"  (the closing quote makes the match exact).
"""
import csv, sys
from pathlib import Path

csv.field_size_limit(sys.maxsize)
DATA = Path("/scratch/jshen/data/pansoma_v2_tensors")
SAMPLES = ["HG008T_PacBio", "HG008T_ONT", "HG008T_Illumina", "COLO829T_Illumina", "COLO829T_fiberseq", "COLO829T_ONT"]
HERE = Path(__file__).resolve().parent
for sample in SAMPLES:
    status = {r["truth_id"]: r["status"] for r in csv.DictReader(open(DATA / sample / "tensors/somatic.recall.tsv"), delimiter="\t")}
    keys = {"SNV": set(), "INDEL": set()}
    for t in csv.DictReader(open(DATA / sample / "truth/somatic.graph.tsv"), delimiter="\t"):
        if status[t["truth_id"]] == "filtered":
            keys["SNV" if t["kind"] == "SNP" else "INDEL"].update(k for k in t["keys"].split(",") if k)
    for kind, found in keys.items():
        with open(HERE / f"{sample}.{kind}.patterns", "w") as out:
            out.writelines(f'{{"candidate_id": "{k}"\n' for k in sorted(found))
        print(sample, kind, len(found))
