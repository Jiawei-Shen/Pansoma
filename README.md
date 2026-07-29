# Pansoma

Pansoma is a research pipeline for generating machine-learning-ready variant tensors from pangenome graph alignments. It takes sequencing reads aligned to a pangenome graph, extracts graph-node pileups around candidate variants, writes sharded NumPy tensors, labels them against truth sets when available, and provides bundled model training and inference code.

The project is organized so that reproducible pipeline entry points live in `scripts/`, reusable code lives in `src/`, model code lives in `machine_learning/`, and older exploratory scripts remain available under `experiments/legacy/`.

![Pansoma workflow](./figures/Pansoma_Fig1.png)

## How to use Pansoma
To use Pansoma, it follows two major stage: 1.Data preprocessing: Create Node Pile Up data. 2. Data Inference: Use the model we provided to inference.
We also provide a model training so that you can use your own data to train Pansoma to get your Pansoma model.
We provide a Docker that contains everything, running commands, and the step to step explanations.
###

## Repository Layout

```text
Pansoma/
├── configs/                    # Example YAML configs for graph resources and tensor runs
├── cpp/                        # C++ extension code and legacy C++ experiments
├── data/                       # Placeholder folders for local raw/interim/processed data
├── docker/                     # Docker build files
├── docs/                       # Pipeline, data format, and migration documentation
├── experiments/legacy/         # Preserved old scripts and research versions
├── machine_learning/pansoma_net/              # Bundled PyTorch model training and inference code
├── scripts/                    # Main command-line entry points
├── src/pangenome_ml_data_generation/
│   ├── alignment/              # GAM/alignment filtering helpers
│   ├── analysis/               # Statistics and normalization utilities
│   ├── graph/                  # Graph path, node, and coordinate helpers
│   ├── io/                     # Format-specific I/O helpers
│   ├── pileup/                 # Pileup/candidate extraction logic
│   ├── plotting/               # Figure and tensor visualization helpers
│   ├── tensors/                # Tensor builders and visualization
│   └── variants/               # VCF, truth-set, and AF utilities
└── tests/                      # Unit/integration test placeholders
```

Large files such as FASTQ, GAM, GBZ/GFA, VCF, BAM, `.dat/.idx`, and `.npy` shards should stay outside git. The `data/` tree is present as a local staging convention, not as a place to commit large artifacts.

## Installation

Create the data-generation environment:

```bash
conda env create -f environment.yml
conda activate pangenome-ml-data-generation
```

If the environment already exists, synchronize its dependencies with the
repository before continuing:

```bash
conda env update -n pangenome-ml-data-generation -f environment.yml
conda activate pangenome-ml-data-generation
```

Do not use the Conda `base` environment for Pansoma. Confirm that Python and
the legacy VG protobuf dependency come from the project environment:

```bash
which python
python -c "import google.protobuf; print(google.protobuf.__version__)"
```

The protobuf version must be `3.20.3`.

External command-line tools are also required for the full pipeline:

```text
vg
jq
bgzip/tabix, usually through htslib or pysam
```

Platform notes are in [docs/platform_support.md](docs/platform_support.md).
In short, Linux/HPC is the recommended target for full production runs, while
macOS is supported for development, Python utilities, and local extension
builds. The FASTQ-to-GAM stage on any platform requires `vg` plus matching
`.gbz`, `.min`, and `.dist` graph indexes. Long-read Giraffe mapping also
requires a matching `.zipcodes` index.

Build the C++ `fast_writer` extension before generating `.dat/.idx` files:

```bash
bash scripts/build_fast_writer.sh
```

`fast_writer` is platform- and Python-version-specific. Build it separately on
each Linux or macOS machine; compiled `.so` files are not stored in git. The
build script verifies the resulting import and prints its installed path.

Check a local machine:

```bash
python scripts/check_environment.py
```

For a full production node where `vg` and other command-line tools must be
available:

```bash
python scripts/check_environment.py --strict-external
```

## Docker

The primary container runs the complete 5-channel workflow, including `vg`
alignment, GAM preprocessing, tensor generation, training, inference, and VCF
output:

```text
docker/Dockerfile
```

Build it from the repository root:

```bash
docker build --platform linux/amd64 -f docker/Dockerfile -t pansoma:latest .
```

Validate the runtime:

```bash
docker run --rm --gpus all pansoma:latest doctor
```

For WUSTL RIS/LSF, the image supports `/bin/bash` as the submitted command and
keeps all Python dependencies in `/opt/venv`. A complete interactive `bsub`
example and cluster checks are included in [docs/docker.md](docs/docker.md).

The image exposes two end-to-end commands:

```bash
docker run ... pansoma:latest train [inputs and training parameters]
docker run ... pansoma:latest infer [inputs, checkpoint, and output VCF]
```

FASTQ files are not sufficient by themselves: both commands also require
matching `.gbz`, `.min`, `.dist`, and GFA resources. Pansoma now generates the
GRCh38 coordinate map and chromosome node filters automatically and reuses
them through a graph-fingerprinted cache. Training additionally requires a
truth VCF in `.vcf` or `.vcf.gz` format; compression and Tabix indexing are
automatic. Inference requires a model checkpoint. Complete mount layouts and
commands are documented in [docs/docker.md](docs/docker.md).

The previous ML-only image remains at `docker/Dockerfile.ml` for legacy use.

## End-To-End Pipeline

Full command examples are in [docs/pipeline.md](docs/pipeline.md).

### 0. GBZ To GFA

If you have a GBZ graph but not its matching GFA, export one with `vg`:

```bash
vg convert -f --no-translation \
  /path/to/graph.gbz \
  > /path/to/graph.gfa
```

Run this once for each graph. `--no-translation` preserves the internal `vg`
node IDs so the GFA remains compatible with alignments produced against the
GBZ. See the official
[`vg giraffe` graph guidance](https://github.com/vgteam/vg/wiki/Giraffe-best-practices#graphs).

### 1. FASTQ To GAM

Align reads to the pangenome graph with `vg giraffe`.

Short reads (paired-end Illumina):

```bash
vg giraffe \
  -Z /path/to/graph.gbz \
  -m /path/to/graph.min \
  -d /path/to/graph.dist \
  -f /path/to/read_1.fq.gz \
  -f /path/to/read_2.fq.gz \
  -b default \
  -t 12 \
  -p \
  > /path/to/sample.gam
```

Long reads (single-end PacBio HiFi):

```bash
vg giraffe \
  -Z /path/to/graph.gbz \
  -m /path/to/graph.longread.withzip.min \
  -z /path/to/graph.longread.zipcodes \
  -d /path/to/graph.dist \
  -f /path/to/long_reads.fq.gz \
  -b hifi \
  -t 12 \
  -p \
  > /path/to/sample.gam
```

Use `-b hifi` for PacBio HiFi or `-b r10` for Oxford Nanopore R10. Long reads
use one FASTQ file and require long-read minimizer and zipcode indexes. See the
official [`vg giraffe` documentation](https://github.com/vgteam/vg/wiki/Giraffe-best-practices)
for additional options.

### 2. GAM To `.dat/.idx`

Find graph nodes that contain imperfect read alignments:

```bash
python -u scripts/find_unperfect_nodes.py sample.gam \
  --output sample.unperfect_nodes.pkl \
  --output_format pickle \
  --milestone 10000000 \
  --threads 12
```

Build the packed node-read store:

```bash
python -u scripts/build_dat_idx.py \
  sample.gam \
  sample.unperfect_nodes.pkl \
  sample.unperfect_nodes \
  --milestone 1000000 \
  --threads 12
```

This writes:

```text
sample.unperfect_nodes.dat
sample.unperfect_nodes.idx
```

### 3. Graph Node Mapping

Build per-chromosome component and GRCh38 path filters:

```bash
cd /path/to/Pansoma

GBZ=/path/to/graph.gbz \
GFA=/path/to/graph.gfa \
OUTDIR=/path/to/chr_node_filters \
CHROMOSOMES="1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22" \
bash ./scripts/build_chr_node_filters.sh
```

`GBZ` and `GFA` are required and must describe the same graph. `OUTDIR`
contains one directory per chromosome plus `summary.tsv`. If `OUTDIR` is
omitted, output is written under
`tmp/chr_component_vs_GRCh38_summary` in the repository.

Build or filter node JSON resources:

```bash
python -u scripts/build_grch38_path_json.py graph.gfa '^W\tGRCh38\t0\tchr1' \
  -o chr1.GRCh38.nodes.json

python -u scripts/filter_node_json.py \
  chr1.GRCh38.nodes.json \
  sample.unperfect_nodes.idx \
  chr1.filtered.nodes.json
```

For whole-testing-set generation from `.idx` plus GFA:

```bash
python -u scripts/build_node_json.py \
  --gfa graph.gfa \
  --idx sample.unperfect_nodes.idx \
  --out candidate_nodes.json
```

### 4. Tensor Generation And Labeling

Generate sharded variant-centered tensors:

```bash
python -u scripts/generate_testing_tensors.py \
  sample.unperfect_nodes.dat \
  sample.unperfect_nodes.idx \
  tensors_chr1 \
  candidate_nodes.json \
  --chr_nodes chr1.component.nodes.raw.txt chr1.GRCh38_path.nodes.raw.txt \
  --num_workers 8 \
  --variant_type snp \
  --view 0 \
  --min_af 0.08 \
  --shard_size 32768
```

Label tensors against a truth VCF:

```bash
python -u scripts/label_tensors.py \
  tensors_chr1/variant_summary.ndjson \
  candidate_nodes.json \
  truth.vcf.gz \
  --chr chr1 \
  --data-dir tensors_chr1
```

### 5. Model Training

Model code is bundled in:

```text
machine_learning/pansoma_net/
```

Training wrapper:

```bash
sbatch scripts/slurm/train_pansoma_net.sh
```

The wrapper defaults to `machine_learning/pansoma_net` but accepts overrides:

```bash
PANSOMA_NET_DIR=/path/to/pansoma_net \
TRAIN_DATA_PATHS_FILE=train_data_dir.txt \
VAL_DATA_PATHS_FILE=val_data_dir.txt \
OUT_DIR=/path/to/output_model \
sbatch scripts/slurm/train_pansoma_net.sh
```

### 6. Model Inference

Inference wrapper:

```bash
INPUT_DIR=/path/to/tensor_shards \
CKPT=/path/to/model.pth \
OUT_PREFIX=/path/to/results/pansoma_sample \
MAP_JSON=/path/to/candidate_nodes.json \
VARIANT_SUMMARY=/path/to/variant_summary.ndjson \
sbatch scripts/slurm/infer_pansoma_net.sh
```

## Key Scripts

```text
scripts/find_unperfect_nodes.py            find nodes with imperfect reads
scripts/build_fast_writer.sh               build C++ writer extension
scripts/build_dat_idx.py                   GAM + node set -> .dat/.idx
scripts/build_chr_node_filters.sh          chromosome component/path node filters
scripts/build_grch38_path_json.py          GRCh38 path node JSON
scripts/filter_node_json.py                filter node JSON by idx/chrom/truth VCF
scripts/build_node_json.py                 candidate node JSON from idx + GFA
scripts/generate_testing_tensors.py        sharded tensor generation
scripts/label_tensors.py                   truth-VCF tensor labeling
scripts/classify_tensors.py                organize true/false tensor datasets
scripts/visualize_tensor.py                tensor visualization
scripts/pansoma_workflow.py                end-to-end Docker train/infer CLI
```

## Data Formats

See [docs/data_formats.md](docs/data_formats.md) for `.idx`, candidate-node JSON, and chromosome node-filter formats.

## Development Notes

- `experiments/legacy/` preserves older script versions for traceability.
- `docs/legacy_mapping.md` maps old script names to the new entry points.
- The current package modules under `src/` are seeded from proven scripts and are intended for progressive refactoring.

## License

See [LICENSE](LICENSE).
