"""Check that CPU inference reproduces the logged v1 GPU run: e053 checkpoint on v1 chr1 *testing* tensors vs the PROB
in the May-14 e053 VCF (SHARD/IDX link), and that val tensors equal testing tensors for the same (node, variant_key)."""
import gzip, json, sys, time
import numpy as np
import torch
torch.set_num_threads(4)
sys.path.insert(0, "/scratch/jshen/Github/Pansoma/machine_learning/pansoma_net")
from mynet import ConvNeXtCBAMClassifier  # noqa
MEAN = torch.tensor([18.41781616, 12.64912987, -0.54525274, 24.72385406, 4.69061136]).view(1, 5, 1, 1)
STD = torch.tensor([25.02832222, 14.80963230, 0.61813378, 29.97283554, 7.92317915]).view(1, 5, 1, 1)
E053 = "/scratch/jshen/Github/GoogleNet/HG008T_GIAB_AF-HPRC_CE_Large_Model_V2_weight200/model_e053_f1_0.1733.pth"
TEST = "/scratch/jshen/data/Pansoma/HG008_GIAB/AF_HPRC/5ch_testing_data_SNV/5ch_testing_data_SNV_chr1"
VAL = "/scratch/jshen/data/Pansoma/HG008_GIAB/AF_HPRC/5ch_training_data_SNV/val/SNV_chr1_5chan_tensor_dataset"
VCF = "/scratch/jshen/Pansoma_testing_results_V2/HG008_GIAB/AF_HPRC/pansoma_HG008T_WGS_ALL_chr/pansoma_HG008T_WGS_chr1/pansoma-to_SNV_pansoma_HG008T_WGS_chr1.linear.vcf.gz"
rng = np.random.default_rng(7)
recs = []
for line in gzip.open(VCF, "rt"):
    if line.startswith("#"):
        continue
    f = line.split("\t")
    info = dict(kv.split("=", 1) for kv in f[7].split(";"))
    recs.append((int(info["SHARD"]), int(info["IDX"]), float(info["PROB"]), info["NID"], f[1], f[3], f[4]))
pick = [recs[i] for i in rng.choice(len(recs), 40, replace=False)]
# 16 random test tensors not in the VCF (expect PROB < 0.5)
invcf = {(s, i) for s, i, *_ in recs}
test_meta = [json.loads(l) for l in open(f"{TEST}/variant_summary.ndjson")]
others = [m for m in (test_meta[i] for i in rng.choice(len(test_meta), 200, replace=False))
          if (m["shard_index"], m["index_within_shard"]) not in invcf][:16]
ck = torch.load(E053, map_location="cpu", mmap=True, weights_only=False)
m = ConvNeXtCBAMClassifier(in_channels=5, class_num=2, depths=[3, 3, 27, 3], dims=[192, 384, 768, 1536])
m.load_state_dict(ck["model_state_dict"], strict=True); m.eval()
print("e053 epoch", ck.get("epoch"), "best_f1", ck.get("best_f1_true"))
sh = {}
def get(s, i):
    if s not in sh:
        sh[s] = np.load(f"{TEST}/shard_{s:05d}_data.npy", mmap_mode="r")
    return np.asarray(sh[s][i])
xs = [get(s, i) for s, i, *_ in pick] + [get(m_["shard_index"], m_["index_within_shard"]) for m_ in others]
t = time.time()
with torch.inference_mode():
    p = []
    for b in range(0, len(xs), 16):
        xb = (torch.from_numpy(np.stack(xs[b:b + 16])).float() - MEAN) / STD
        p.append(torch.softmax(m(xb), 1)[:, 1].numpy())
p = np.concatenate(p)
print("infer s", round(time.time() - t, 1))
vp = np.array([r[2] for r in pick])
d = np.abs(p[:40] - vp)
print("VCF PROB vs CPU: max abs diff", d.max().round(5), "median", np.median(d).round(6), "CPU>=0.5:", int((p[:40] >= 0.5).sum()), "/40")
print("non-VCF tensors CPU p>=0.5:", int((p[40:] >= 0.5).sum()), "/", len(others), "max", p[40:].max().round(4))
# val vs test tensor identity on (node, variant_key)
val_meta = {}
for l in open(f"{VAL}/variant_summary_classified.ndjson"):
    r = json.loads(l); val_meta[(r["node_id"], r["variant_key"])] = (r["shard_index"], r["index_within_shard"])
same = tot = 0
vsh = {}
for mm in (test_meta[i] for i in rng.choice(len(test_meta), 300, replace=False)):
    k = (mm["node_id"], mm["variant_key"])
    if k not in val_meta:
        continue
    s, i = val_meta[k]
    if s not in vsh:
        vsh[s] = np.load(f"{VAL}/HG008T_Illumina_chr1_{s:05d}_data.npy", mmap_mode="r")
    tot += 1; same += bool(np.array_equal(np.asarray(vsh[s][i]), get(mm["shard_index"], mm["index_within_shard"])))
print("val tensor == test tensor:", same, "/", tot)
