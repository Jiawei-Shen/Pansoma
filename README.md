# Pansoma

Pansoma is a research pipeline for generating machine-learning-ready variant tensors from pangenome graph alignments. It takes sequencing reads aligned to a pangenome graph, extracts graph-node pileups around candidate variants, writes sharded NumPy tensors, labels them against truth sets when available, and provides bundled model training and inference code.

The project is organized so that reproducible pipeline entry points live in `scripts/`, reusable code lives in `src/`, model code lives in `machine_learning/`, and older exploratory scripts remain available under `experiments/legacy/`.

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
`.gbz`, `.min`, and `.dist` graph indexes.

Build the C++ `fast_writer` extension before generating `.dat/.idx` files:

```bash
bash scripts/build_fast_writer.sh
```

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

The bundled Dockerfile is adapted from the original model repository and is focused on the PyTorch/ML environment:

```text
docker/Dockerfile.ml
machine_learning/pansoma_net/Dockerfile
```

Build it from the repository root:

```bash
docker build -f docker/Dockerfile.ml -t pansoma-ml .
```

Run it with GPU support and mount the project plus data directories:

```bash
docker run --gpus all -it --rm \
  -v "$PWD":/workspace/Pansoma \
  -v /scratch:/scratch \
  pansoma-ml
```

Inside the container:

```bash
cd /workspace/Pansoma
```

The Dockerfile installs the ML requirements from `machine_learning/pansoma_net/requirements.txt`.

## End-To-End Pipeline

Full command examples are in [docs/pipeline.md](docs/pipeline.md).

### 1. FASTQ To GAM

Align reads to the pangenome graph with `vg giraffe`:

```bash
GBZ=/path/to/graph.gbz \
MIN_INDEX=/path/to/graph.min \
DIST_INDEX=/path/to/graph.dist \
FASTQ1=/path/to/read_1.fq.gz \
FASTQ2=/path/to/read_2.fq.gz \
READ_TYPE=illumina \
OUT_GAM=/path/to/sample.gam \
bash scripts/run_giraffe.sh
```

Use `READ_TYPE=hifi` for PacBio HiFi and `READ_TYPE=ont` or `READ_TYPE=r10` for ONT.

### 2. GAM To `.dat/.idx`

Find graph nodes that contain imperfect read alignments:

```bash
python -u scripts/find_unperfect_nodes.py sample.gam \
  --output sample.unperfect_nodes.pkl \
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
bash scripts/build_chr_node_filters.sh
```

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
scripts/run_giraffe.sh                     FASTQ -> GAM
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
```

## Data Formats

See [docs/data_formats.md](docs/data_formats.md) for `.idx`, candidate-node JSON, and chromosome node-filter formats.

## Development Notes

- `experiments/legacy/` preserves older script versions for traceability.
- `docs/legacy_mapping.md` maps old script names to the new entry points.
- The current package modules under `src/` are seeded from proven scripts and are intended for progressive refactoring.

## License

See [LICENSE](LICENSE).
