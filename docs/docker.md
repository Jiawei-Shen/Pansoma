# End-to-End Docker Workflow

The primary Docker image runs the complete 5-channel Pansoma workflow on
Linux: FASTQ alignment, GAM preprocessing, tensor generation, model training,
model inference, and indexed linear VCF output.

## Required Inputs

FASTQ files are sample inputs, but the pipeline also requires graph/reference
resources that must match each other:

| Input | Training | Inference | Purpose |
| --- | --- | --- | --- |
| FASTQ R1 and optional R2 | Required | Required | Sequencing reads |
| `.gbz`, `.min`, `.dist` | Required | Required | `vg giraffe` graph indexes |
| `.zipcodes` | Long reads | Long reads | Long-read Giraffe distance hints |
| GFA | Required | Required | Graph node sequences |
| Truth VCF (`.vcf` or `.vcf.gz`) | Required | No | Supervised training labels |
| PansomaNet checkpoint | No | Required | Trained model used for inference |

If a matching GFA is not already available, generate it once from the GBZ as
shown in [Step 0 of the pipeline guide](pipeline.md#0-gbz-to-gfa).

Pansoma generates these intermediate resources automatically:

```text
WORK_DIR/
├── graph_resources/<graph-fingerprint>/
│   ├── chromosome_filters/
│   ├── coordinates/
│   └── manifest.json
├── nodes/
│   ├── SAMPLE.<graph-id>.candidate_nodes.raw.json
│   └── SAMPLE.<graph-id>.candidate_nodes.json
└── truth/
    ├── truth.<truth-id>.vcf.gz
    └── truth.<truth-id>.vcf.gz.tbi
```

The optional `--resource-cache` argument places graph intermediates in a
shared persistent directory. Without it, they are stored under
`WORK_DIR/graph_resources`. The graph fingerprint prevents resources from a
different GBZ/GFA pair from being reused. The `build-node-filters` command is
still available for diagnostics, but it is not required before `train` or
`infer`.

## Build

Build from the repository root:

```bash
docker build --platform linux/amd64 \
  -f docker/Dockerfile \
  -t pansoma:latest .
```

The image is based on CUDA 12.1 and installs PyTorch, the Pansoma Python/C++
code, `vg`, `samtools`, `bcftools`, `tabix`, and `jq`. The pinned `vg` binary is
the official Linux x86-64 release, so the production image targets
`linux/amd64`. NVIDIA Container Toolkit is required for GPU use.

Validate the image:

```bash
docker run --rm --gpus all pansoma:latest doctor
```

## WUSTL RIS/LSF Interactive Use

Use a fixed image version for cluster jobs so an older cached `latest` image
cannot be selected unexpectedly. This opens an interactive shell with six
GPUs:

```bash
bsub -Is \
  -G compute-epigenome-condo \
  -n 4 \
  -R "select[gpuhost] rusage[mem=200G]" \
  -gpu "num=6" \
  -q epigenome-interactive \
  -a 'docker(jiaweiwustk/pansoma:0.1.3)' \
  /bin/bash
```

Inside the container, verify the isolated Python runtime and external tools:

```bash
which python
python -c 'import numpy, torch; print(numpy.__version__, torch.__version__)'
pansoma doctor
nvidia-smi
```

`python` should resolve to `/opt/venv/bin/python`. The `pansoma` command can be
used from an interactive shell for `doctor`, `build-node-filters`, `train`, and
`infer`.

The cluster must expose input and output directories inside the container.
Use the RIS-supported mount or volume options for your environment, then pass
the resulting container paths to Pansoma. Confirm visibility before a long
run with `ls` on the FASTQ, graph, model, and output paths.

## Train From FASTQ

The truth VCF can be uncompressed, gzip-compressed, or already BGZF-compressed
and indexed. When needed, Pansoma writes a BGZF copy and `.tbi` index under
`--work-dir`. The example mounts inputs read-only while keeping work, cache,
and model outputs on persistent host volumes:

```bash
docker run --rm --gpus all --shm-size=32g \
  -v /path/to/reads:/inputs/reads:ro \
  -v /path/to/graph:/inputs/graph:ro \
  -v /path/to/truth:/inputs/truth:ro \
  -v /path/to/pansoma-work:/work \
  -v /path/to/models:/output \
  pansoma:latest train \
  --sample HG008T \
  --fastq1 /inputs/reads/HG008T_R1.fastq.gz \
  --fastq2 /inputs/reads/HG008T_R2.fastq.gz \
  --mapper-preset default \
  --gbz /inputs/graph/hprc.gbz \
  --min-index /inputs/graph/hprc.min \
  --dist-index /inputs/graph/hprc.dist \
  --gfa /inputs/graph/hprc.gfa \
  --truth-vcf /inputs/truth/HG008T.truth.vcf \
  --resource-cache /work/graph-cache \
  --work-dir /work/HG008T-train \
  --output-dir /output/HG008T-pansoma-net \
  --variant-type snp \
  --val-chromosomes chr1 \
  --epochs 70 \
  --batch-size 32 \
  --learning-rate 0.0001 \
  --nproc-per-node 1
```

By default, chromosomes 1-22 are processed. `chr1` is routed to validation
and the remaining chromosomes are routed to training. With multiple GPUs,
pass all devices to Docker and set `--nproc-per-node` to the GPU count.

The stable model path is:

```text
/output/HG008T-pansoma-net/best_model.pth
```

It points to the best F1 checkpoint retained by the trainer.

## Infer From FASTQ

```bash
docker run --rm --gpus all --shm-size=32g \
  -v /path/to/reads:/inputs/reads:ro \
  -v /path/to/graph:/inputs/graph:ro \
  -v /path/to/models:/inputs/models:ro \
  -v /path/to/pansoma-work:/work \
  -v /path/to/results:/output \
  pansoma:latest infer \
  --sample HG008T \
  --fastq1 /inputs/reads/HG008T_R1.fastq.gz \
  --fastq2 /inputs/reads/HG008T_R2.fastq.gz \
  --mapper-preset default \
  --gbz /inputs/graph/hprc.gbz \
  --min-index /inputs/graph/hprc.min \
  --dist-index /inputs/graph/hprc.dist \
  --gfa /inputs/graph/hprc.gfa \
  --checkpoint /inputs/models/best_model.pth \
  --resource-cache /work/graph-cache \
  --work-dir /work/HG008T-infer \
  --output-vcf /output/HG008T.pansoma.vcf.gz \
  --variant-type snp \
  --batch-size 64
```

For long-read samples, provide one FASTQ and select the matching native
Giraffe preset:

```bash
# PacBio HiFi: vg giraffe -b hifi
pansoma infer ... \
  --fastq1 sample.hifi.fastq.gz \
  --mapper-preset hifi \
  --min-index graph.longread.withzip.min \
  --zipcode-index graph.longread.zipcodes

# ONT R10: vg giraffe -b r10
pansoma infer ... \
  --fastq1 sample.ont-r10.fastq.gz \
  --mapper-preset r10 \
  --min-index graph.longread.withzip.min \
  --zipcode-index graph.longread.zipcodes
```

Do not pass `--fastq2` for long-read data. PacBio CLR and ONT R9 are not
accepted. The workflow calls `vg giraffe` directly. Short reads use
`--mapper-preset default` unless another native Giraffe preset is selected.

Inference writes:

```text
HG008T.pansoma.vcf.gz
HG008T.pansoma.vcf.gz.tbi
```

Per-chromosome linear VCFs are produced first and then streamed into one
BGZF-compressed, Tabix-indexed VCF.

## Restart Behavior

All expensive intermediate files are stored under `--work-dir`. If a command
is restarted with the same arguments and work directory, completed GAM,
`.pkl`, `.dat/.idx`, candidate node-map, tensor, label, and per-chromosome
inference stages are reused. Graph resources are keyed by a GBZ/GFA
fingerprint, and truth labels are keyed by a truth-VCF fingerprint. Use a new
work directory when changing FASTQ inputs or tensor-generation parameters.

Run command-specific help with:

```bash
docker run --rm pansoma:latest train --help
docker run --rm pansoma:latest infer --help
```
