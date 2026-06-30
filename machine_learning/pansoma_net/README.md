# Pansoma ML Module

This folder contains the bundled PyTorch model code maintained as the `pansoma_net` model module. It is used for Pansoma tensor training and inference.

## Main Files

```text
mynet.py                                      model definitions used by Pansoma scripts
pansoma_net.py                                 original GoogLeNet model implementation
scripts/train_5channels_npy_pansoma.py       5-channel training
scripts/train_6channels_npy_pansoma.py       6-channel training
scripts/train_6channels_npy_pansoma_DDP.py   distributed 6-channel training
scripts/test_5channels_npy_pansoma.py        5-channel inference
scripts/test_6channels_npy_pansoma*.py       6-channel inference variants
scripts/utilities/                           conversion/evaluation helpers
```

## Environment

Use the ML-specific files here when training or inference requires CUDA/PyTorch:

```bash
conda env create -f "machine_learning/pansoma_net/environment.yml"
conda activate pansoma_net_env
```

or install with:

```bash
pip install -r "machine_learning/pansoma_net/requirements.txt"
```

## Docker

The ML Dockerfile is available at:

```text
machine_learning/pansoma_net/Dockerfile
docker/Dockerfile.ml
```

Build from the repository root:

```bash
docker build -f docker/Dockerfile.ml -t pansoma-ml .
```

## SLURM Wrappers

The top-level wrappers default to this folder:

```text
scripts/slurm/train_pansoma_net.sh
scripts/slurm/infer_pansoma_net.sh
```

Override `PANSOMA_NET_DIR` if you want to use an external checkout instead.

