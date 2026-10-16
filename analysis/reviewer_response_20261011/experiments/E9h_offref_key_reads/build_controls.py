"""E9h control queries: the hit rules and query lengths of classify.py applied to sequences of known origin, the events'
own assembly alleles. For every event with a key (windows.tsv of build_queries.py), from its 201-bp windows:
  T context  the somatic allele with HALF bases of HG008-N flank on each side (somatic by construction: expected
             tumor_only)
  N context  the normal allele with the same flanks (the patient's germline haplotype: expected normal)
HALF = 30 and 50, as the key queries; site = the event's first allele base (the anchor base of an indel).

usage: python build_controls.py DATA_DIR   (reads windows.tsv and queries.fa; writes controls.fa, controls.tsv)"""
import csv
import sys
from pathlib import Path

HALVES = (30, 50)


def main():
    d = Path(sys.argv[1])
    seqs, name = {}, None
    for line in open(d / "queries.fa"):
        if line[0] == ">":
            name = line[1:].strip()
        else:
            seqs[name] = line.strip()
    n = 0
    with open(d / "controls.fa.tmp", "w") as fa, open(d / "controls.tsv.tmp", "w") as tsv:
        tsv.write("qid\tevent\twhich\thalf\tqlen\tsite\n")
        for w in csv.DictReader(open(d / "windows.tsv"), delimiter="\t"):
            seq, flank = seqs[w["qid"]], int(w["flank"])
            allele = w["alt"] if w["which"] == "T" else w["ref"]
            for h in HALVES:
                q = seq[max(0, flank - h):flank + len(allele) + h]
                site = flank - max(0, flank - h)
                qid = f"c_{w['event']}_{w['which']}_h{h}"
                fa.write(f">{qid}\n{q}\n")
                tsv.write(f"{qid}\t{w['event']}\t{w['which']}\t{h}\t{len(q)}\t{site}\n")
                n += 1
    for stem in ("controls.fa", "controls.tsv"):
        (d / f"{stem}.tmp").replace(d / stem)
    print(f"{n} control queries -> {d}/controls.fa")


if __name__ == "__main__":
    main()
