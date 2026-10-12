"""pansoma_net_v2: PansomaNetV2 for merged Pansoma tensors (8 x 200 x 101 int8, labels 0/1/2, -1 ignored).

encode.py      int8 channels -> 36 float planes on the GPU (one-hot, masked z-score, coverage, derived planes)
data.py        merged tensor sets -> cached index -> samples (x, row blocks, scalars, label); block split, truth ids
bed.py         BED regions and the labeller's confident-region test
chunks.py      zstd chunk cache of the training shards and ChunkLoader (python -m pansoma_net_v2.chunks)
model.py       encoder -> 1x1 front -> ConvNeXt-CBAM backbone (mynet's ConvNeXtCBAMClassifier, blocks v2 / v1) -> 3 classes
train.py       training (one GPU or torchrun DDP)
predict.py     per-tensor probabilities, the SNV off-reference rescue and the test metrics
metrics.py     per-tensor and truth-level metrics, decision rules (numpy only)
summary.py     per-epoch table, test results and vcfeval reports of training runs
env.py         Triton libcuda fix for torch.compile on the cluster's GPU nodes
graph_vcf.py   predictions -> a VCF in graph coordinates
linear_vcf.py  graph VCF -> GRCh38 VCF (plus the records without a GRCh38 position)
vcfeval.py     panel-of-normals filter and rtg vcfeval against a somatic truth VCF
vcfio.py       plain-text VCF helpers of the three VCF steps
"""
