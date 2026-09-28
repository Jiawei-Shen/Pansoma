"""pansoma_net_v2: PansomaNetV2 for merged Pansoma tensors (8 x 200 x 101 int8, labels 0/1/2, -1 ignored).

encode.py  int8 channels -> 36 float planes on the GPU (one-hot, masked z-score, coverage, derived planes)
data.py    merged tensor sets -> cached index -> samples (x, row blocks, label); chromosome splits
model.py   encoder -> 1x1 front -> ConvNeXt-CBAM backbone (from pansoma_net/mynet.py) -> 3 classes
train.py   training (one GPU or torchrun DDP); predict.py  per-tensor probabilities and metrics
"""
