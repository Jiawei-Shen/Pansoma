"""Copy the two tensors drawn in fig_R1.m2_tensor_across_nodes.png and indel/fig_R1.m2_INDEL_span_01.png out of their
chr1 shards, so the figures can be redrawn without the tensor sets (read-only on the sets). Per tensor:
  <name>.npy             the tensor as a one-sample shard, shape (1, 8, 200, 101) int8, a byte copy of the shard row
  <name>.summary.ndjson  its chr1_variant_summary.ndjson record, with shard_index / index_within_shard set to 0 (the
                         source values kept as source_shard_file / source_shard_index / source_index_within_shard)
  <name>.label.json      its chr1_labels.ndjson record and the value in chr1_shard_00000_labels.npy
Redraw: see render.sh. Run: /wanglab/jshen/anaconda3/envs/hunyuanvideo15/bin/python extract.py
"""
import json
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
D = Path("/scratch/jshen/data/pansoma_v2_tensors")
TENSORS = (  # name, tensor directory, shard, index within the shard, candidate id (checked)
    ("R1.m2_tensor_across_nodes", "HG008T_PacBio/tensors/SNV", "chr1_shard_00000_data.npy", 32257, "3130240:147:SNP:G>A"),
    ("R1.m2_INDEL_span_01", "HG008T_Illumina_d46/tensors/INDEL", "chr1_shard_00000_data.npy", 271,
     "1811743:0:DEL:TATATATA>@1811744+1811746+1811747+1811748+1811749"),
)


def find(path, key, value, index):
    """The record with r[key] == value and index_within_shard == index (summary: shard_index, labels: shard_file)."""
    for line in open(path):
        r = json.loads(line)
        if r.get(key) == value and r.get("index_within_shard") == index:
            return r
    raise LookupError((path, value, index))


for name, sub, shard, index, cid in TENSORS:
    src = D / sub
    shard_index = int(shard.split("_")[2])
    data = np.load(src / shard, mmap_mode="r")
    tensor = np.ascontiguousarray(data[index:index + 1])
    rec = find(src / "chr1_variant_summary.ndjson", "shard_index", shard_index, index)
    lab = find(src / "chr1_labels.ndjson", "shard_file", shard, index)
    assert rec["candidate_id"] == cid == lab["candidate_id"], (name, rec["candidate_id"], lab["candidate_id"])
    rec.update(source_set=str(src), source_shard_file=shard, source_shard_index=shard_index,
               source_index_within_shard=index, shard_index=0, index_within_shard=0)
    np.save(HERE / f"{name}.npy", tensor)
    (HERE / f"{name}.summary.ndjson").write_text(json.dumps(rec) + "\n")
    label = np.load(src / shard.replace("_data.npy", "_labels.npy"), mmap_mode="r")[index].item()
    (HERE / f"{name}.label.json").write_text(json.dumps(dict(labels_npy=label, labels_ndjson=lab), indent=1) + "\n")
    print(name, tensor.shape, tensor.dtype, "label", label, lab.get("reason"))
