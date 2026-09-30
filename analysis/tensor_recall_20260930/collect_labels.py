"""Label of the representative tensor of every truth allele with recall status tensor_representative.

The label record's candidate_id is the tensor's representative allele (A1); a truth allele is represented when A1 is
one of its keys. (Records labelled below_*_min_af carry no truth list, so the match goes through the keys.)
Output: labels/<sample>.tsv  truth_id, candidate_id, label, reason
"""
import csv, re, sys
from concurrent.futures import ProcessPoolExecutor
from pathlib import Path

csv.field_size_limit(sys.maxsize)
DATA = Path("/scratch/jshen/data/pansoma_v2_tensors")
SAMPLES = ["HG008T_PacBio", "HG008T_ONT", "HG008T_Illumina", "COLO829T_Illumina", "COLO829T_fiberseq", "COLO829T_ONT"]
HERE = Path(__file__).resolve().parent
CANDIDATE = re.compile(rb'^\{"candidate_id": "([^"]+)"')
LABEL = re.compile(rb'"label": (-?\d+), "label_name": "[^"]*", "reason": "([^"]+)"')


def collect(sample):
    status = {r["truth_id"]: r["status"] for r in csv.DictReader(open(DATA / sample / "tensors/somatic.recall.tsv"), delimiter="\t")}
    owner = {}
    for t in csv.DictReader(open(DATA / sample / "truth/somatic.graph.tsv"), delimiter="\t"):
        if status[t["truth_id"]] == "tensor_representative":
            for k in t["keys"].split(","):
                owner.setdefault(k, []).append(t["truth_id"])
    rows = []
    for kind in ("SNV", "INDEL"):
        for f in sorted((DATA / sample / "tensors" / kind).glob("chr*_labels.ndjson")):
            with f.open("rb") as stream:
                for line in stream:
                    c = CANDIDATE.match(line).group(1).decode()
                    if c in owner:
                        m = LABEL.search(line)
                        for tid in owner[c]:
                            rows.append((tid, c, m.group(1).decode(), m.group(2).decode()))
    (HERE / "labels").mkdir(exist_ok=True)
    with open(HERE / "labels" / f"{sample}.tsv", "w") as out:
        out.write("truth_id\tcandidate_id\tlabel\treason\n")
        out.writelines("\t".join(r) + "\n" for r in rows)
    found = {r[0] for r in rows}
    missing = sum(1 for s in status.values() if s == "tensor_representative") - len(found)
    return f"{sample}: {len(rows)} rows, {len(found)} truths, representative truths without a label row: {missing}"


if __name__ == "__main__":
    with ProcessPoolExecutor(len(SAMPLES)) as pool:
        for line in pool.map(collect, SAMPLES):
            print(line, flush=True)
