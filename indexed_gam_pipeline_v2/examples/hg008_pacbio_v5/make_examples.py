#!/usr/bin/env python3
"""Pick 100 random HG008 PacBio v5 production tensors (34 SNP, 33 DEL, 33 INS) and render them.

Writes NNN_<event>_<candidate>.png plus index.tsv (identity, counts, AF, and the strand
of the rows whose pixels show the candidate ALT: SNP = X op with the ALT base and BQ >= 10
at the anchor, INS = I ops spelling the ALT over the candidate columns, DEL = D ops over them).
Run with the matplotlib-capable interpreter (see ../../README.md, "Rendering PNGs").
"""
import json, random, subprocess, sys
from pathlib import Path
import numpy as np

TENSORS = Path("/scratch/jshen/data/HG008_GIAB/pansoma_v2_tensors/Liss_lab_PacBio_Revio_20240125/v5_tensors")
REPO = Path(__file__).resolve().parents[3]
OUT = Path(__file__).resolve().parent
WANT = {"SNP": 34, "DEL": 33, "INS": 33}
rng = random.Random(20260923)

picked = {k: [] for k in WANT}
for task in rng.sample(range(1024), 1024):
    if all(len(v) >= WANT[k] for k, v in picked.items()):
        break
    for kind in ("SNV", "INDEL"):
        folder = TENSORS / kind / f"task_{task:04d}"
        records = [json.loads(l) for l in open(folder / "variant_summary.ndjson")]
        for m in rng.sample(records, min(2, len(records))):
            if len(picked[m["event_type"]]) < WANT[m["event_type"]]:
                picked[m["event_type"]].append((task, folder, m))

def alt_rows(x, m):
    n, lo, hi = m["selected_alignments"], *m["candidate_columns"]
    if m["event_type"] == "SNP":
        hit = (x[4, :n, 50] == 2) & (x[0, :n, 50] == x[2, :n, 50]) & (x[1, :n, 50] >= 10)
    elif m["event_type"] == "INS":
        hit = ((x[4, :n, lo:hi] == 3) & (x[0, :n, lo:hi] == x[2, :n, lo:hi])).all(axis=1)
    else:
        hit = (x[4, :n, lo:hi] == 4).all(axis=1)
    strand = x[7, :n].max(axis=1)
    return int(hit.sum()), int((hit & (strand == 2)).sum())

rows = ["index\tpng\tcandidate_id\ttask\tcoverage\talt_count\tref_count\tother_count\taf\tselected\t"
        "alt_pixel_rows\talt_pixel_rows_reverse\talt_reverse_fraction\tall_rows_reverse_fraction"]
i = 0
for event in ("SNP", "DEL", "INS"):
    for task, folder, m in picked[event]:
        shard = folder / f"shard_{m['shard_index']:05d}_data.npy"
        x = np.load(shard, mmap_mode="r")[m["index_within_shard"]].astype(np.int64)
        alt, alt_rev = alt_rows(x, m)
        n = m["selected_alignments"]
        rev_all = float((x[7, :n].max(axis=1) == 2).mean()) if n else 0.0
        name = f"{i:03d}_{event}_{m['candidate_id'].replace(':', '_').replace('>', '-')}"[:120] + ".png"
        subprocess.run([sys.executable, str(REPO / "scripts/visualize_tensor.py"), str(shard),
                        "-i", str(m["index_within_shard"]), "-o", str(OUT / name)], check=True,
                       stdout=subprocess.DEVNULL)
        rows.append("\t".join(map(str, [i, name, m["candidate_id"], task, m["coverage"], m["alt_count"],
                                        m["ref_count"], m["other_count"], round(m["af"], 4), n, alt, alt_rev,
                                        round(alt_rev / alt, 3) if alt else "", round(rev_all, 3)])))
        i += 1
(OUT / "index.tsv").write_text("\n".join(rows) + "\n")
print(f"rendered {i} tensors into {OUT}")
