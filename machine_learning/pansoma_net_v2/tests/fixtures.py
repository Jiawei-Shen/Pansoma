"""Small merged tensor sets in the real layout (manifest, shards, labels, summaries with row_groups)."""
import json
from pathlib import Path

import numpy as np

CODES = {"A": 1, "C": 2, "G": 3, "T": 4}


def make_tensor(rng, blocks, alt="A", ref="C", width=101, height=200):
    """One (8, 200, 101) int8 tensor: rows below blocks[3] carry reads (some cells uncovered at read ends),
    the rest is padding. A1 rows carry `alt` at column 50, REF rows `ref`, ALT rows another base, OTHER 0."""
    x = np.zeros((8, height, width), np.int8)
    a1, alt_end, ref_end, rows = blocks
    for r in range(rows):
        lo, hi = sorted(rng.integers(0, width, 2))
        lo, hi = min(lo, 40), max(hi, 60)  # every read covers the site column
        cells = slice(lo, hi + 1)
        x[4, r, cells] = rng.integers(1, 7, hi + 1 - lo)             # operation
        x[0, r, cells] = rng.integers(1, 7, hi + 1 - lo)             # read base
        x[5, r, cells] = rng.integers(1, 5, hi + 1 - lo)             # graph base
        x[1, r, cells] = rng.integers(-1, 61, hi + 1 - lo)           # BQ, -1 = none
        x[3, r, cells] = rng.integers(0, 61)                          # MAPQ (0 is legal)
        x[6, r, cells] = rng.integers(1, 128, hi + 1 - lo)           # path count
        x[7, r, cells] = rng.integers(1, 3)                           # strand
        x[5, r, 50] = CODES[ref]
        if r < a1:
            x[2, r, 50] = x[0, r, 50] = CODES[alt]
        elif r < alt_end:
            x[2, r, 50] = x[0, r, 50] = CODES["G"]
        elif r < ref_end:
            x[2, r, 50] = x[0, r, 50] = CODES[ref]
    return x


def make_tensor_set(root, spec, shard_size=4, seed=0, labels_created="2026-09-28T00:00:00+00:00"):
    """root/<KIND>/ merged directories from spec {kind: {chrom: n}}; returns {kind: [dict per tensor]}."""
    rng = np.random.default_rng(seed)
    truth = {}
    for kind, chroms in spec.items():
        d = Path(root) / kind
        d.mkdir(parents=True, exist_ok=True)
        manifest = dict(status="complete", layout="chromosome-shards", kind=kind, shape=[8, 200, 101], dtype="int8",
                        tensors=sum(chroms.values()), chromosomes={})
        truth[kind] = []
        for chrom, n in chroms.items():
            records, shards = [], []
            for s in range(0, n, shard_size):
                count = min(shard_size, n - s)
                name = f"{chrom}_shard_{s // shard_size:05d}_data.npy"
                data, labels = np.zeros((count, 8, 200, 101), np.int8), np.zeros(count, np.int8)
                for i in range(count):
                    rows = int(rng.integers(3, 40))
                    a1 = int(rng.integers(1, rows))
                    alt_end = a1 + int(rng.integers(0, 3)) if a1 + 2 < rows else a1
                    ref_end = int(rng.integers(alt_end, rows + 1))
                    blocks = [a1, alt_end, ref_end, rows]
                    data[i] = make_tensor(rng, blocks)
                    labels[i] = rng.choice([-1, 0, 0, 1, 2])
                    groups = [dict(start_row=0, end_row=a1, allele="A1")]
                    if alt_end > a1:
                        groups.append(dict(start_row=a1, end_row=alt_end, allele="A2"))
                    if ref_end > alt_end:
                        groups.append(dict(start_row=alt_end, end_row=ref_end, allele="REF"))
                    if rows > ref_end:
                        groups.append(dict(start_row=ref_end, end_row=rows, allele="OTHER"))
                    node = 1000 + 7 * len(truth[kind])
                    cid = f"{node}:5:SNP:C>A"
                    af = round(float(rng.random()), 4)
                    nums = dict(coverage=rows + 3, alt_count=a1, ref_count=ref_end - alt_end, other_count=rows - ref_end,
                                af=af, site_coverage=rows + 5, allele_count=1 + (alt_end > a1),
                                second_allele_af=0.1 if alt_end > a1 else 0.0, event_length=1 + i % 3)
                    records.append(json.dumps(dict(
                        candidate_id=cid, node_id=node, start=5, ref="C", alt="A", event_type="SNP",
                        event_length=nums["event_length"], row_groups=groups,
                        coverage=nums["coverage"], alt_count=nums["alt_count"], ref_count=nums["ref_count"],
                        other_count=nums["other_count"], af=af, site_coverage=nums["site_coverage"], site_id=f"{cid}:SNV",
                        alleles=[dict(candidate_id="9:9:SNP:C>T", event_length=99, coverage=999, alt_count=99,
                                      ref_count=99, other_count=99, af=0.5)],
                        allele_count=nums["allele_count"], second_allele_af=nums["second_allele_af"],
                        parameters=dict(min_af=0.06), shard_index=s // shard_size, index_within_shard=i, chrom=chrom,
                        shard_file=name, source_task=0, source_shard_index=0, source_index_within_shard=7)))
                    truth[kind].append(dict(chrom=chrom, file=name, row=i, label=int(labels[i]), blocks=blocks,
                                            af=af, candidate_id=cid, node=node, numbers=nums, x=data[i]))
                np.save(d / name, data)
                np.save(d / name.replace("_data.npy", "_labels.npy"), labels)
                shards.append(dict(file=name, tensors=count))
            (d / f"{chrom}_variant_summary.ndjson").write_text("\n".join(records) + "\n")
            manifest["chromosomes"][chrom] = dict(tensors=n, shards=shards, summary=f"{chrom}_variant_summary.ndjson")
        (d / "manifest.json").write_text(json.dumps(manifest))
        (d / "labels.manifest.json").write_text(json.dumps(dict(
            version="truth-labels-v6", created=labels_created, tensors=manifest["tensors"],
            labels=dict(ignore=-1, non=0, somatic=1, germline=2))))
    return truth
