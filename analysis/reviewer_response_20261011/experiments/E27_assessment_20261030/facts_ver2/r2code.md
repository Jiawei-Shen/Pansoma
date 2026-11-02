# R2 code availability: install path, minimal example, clean-environment verification (2026-10-30)

Staging dir (git-ignored): `/scratch/jshen/Github/Pansoma/tmp/assessment_20261030/r2code/` (below: `$S`).
Nothing in the shared working tree was edited, added or committed. All runs used fresh clones of the **public GitHub
`develop` branch at 44821f3** (`git clone --depth 1 --branch develop https://github.com/Jiawei-Shen/Pansoma.git`), in `$S/repo*`.

## 0. Findings that matter before anything else

1. **The public default branch is the old code.** GitHub `HEAD` = `main` = e1c0447: the retired `.dat/.idx` pipeline
   (`docker/`, `setup.py`, `environment.yml` with protobuf 3.20.3, `machine_learning/pansoma_net` only). The current
   pipeline (`indexed_gam_pipeline_v4`, `pansoma_net_v2`) is only on `develop` (267 commits ahead of `origin/main`).
   A reviewer who clones the default branch gets the v1 code. This is probably the code they could not run.
2. **The current code does not install as published.** In a fresh conda env, `python -m indexed_gam_pipeline_v4.native compile`
   fails: `compile_decoder` builds with `-static-libstdc++ -static-libgcc` first. That build compiles but its self-test
   crashes: SIGSEGV with conda-forge g++ 15, or "ImportError: NumPy: invalid buffer descriptor!" with the system g++ 11.4.
   The code then raises instead of trying the dynamic build. `run discover` **requires** the native module
   (`discovery.py`: "discover needs the native module"), so nothing can run from a GAM. Verified twice: in a conda-forge env
   (`env_cpu`) and in the pinned env (`env_cpu2`; debug script `$S/debug/native_variants.py`). The dynamic build passes the
   self-test with either compiler. The fix is `$S/patches/native_selftest_fallback.patch` (self-test each link variant,
   then fall back to dynamic). It is applied only in the clones.
3. **The retired env files** (`machine_learning/pansoma_net/{Dockerfile,environment.yml,requirements.txt}`) don't match the
   current code. They pin torch 2.1 / torchvision 0.16 / timm 1.0.15, and the Dockerfile (Python 3.8) and environment.yml
   (Python 3.10) disagree with each other. NumPy 2.3.2 needs Python ≥ 3.11. They have no scipy, protobuf, zstandard,
   pybind11, compiler, rtg or bcftools. I did not try to build them.
4. **The response drafts are public.** The public `develop` branch has `Reviewer_Response_Assessment.md` and
   `analysis/reviewer_response_20261011/points/R*.md`. Anyone with the repo link, reviewers included, can read them.
5. **The v2 checkpoints are not deposited anywhere.** `pretrained_model/Pansoma_zenodo` holds the v1 models.

## 1. Dependency inventory (current code)

| what | used by | version in the developer env (`/wanglab/jshen/anaconda3/bin/python`) |
|---|---|---|
| Python | all | 3.11.3 (Anaconda); NumPy 2.3.2 needs ≥ 3.11 |
| numpy | all | 2.3.2 (pip) |
| scipy | v4 `candidates.py` (hierarchy.linkage) | 1.16.1 |
| protobuf | v4 `vg_pb2.py` (generated, builder API) | 5.29.3 (pip) |
| pysam | v4 BGZF/GAI/VCF/FASTA, labels, PoN script, linear_vcf, site_features | 0.23.0 (PyPI wheel, bundled htslib) |
| pybind11 + C++17 compiler | v4 `native.py` (`_fastdecode`, required by `run discover`) | 3.0.1, system g++ 11.4 |
| sqlite3 (stdlib) | graph index | SQLite 3.41.2 |
| torch (+ torchvision, pulled in by timm) | pansoma_net_v2 | 2.8.0+cu128 / 0.23.0 |
| timm | `model.py` (DropPath, trunc_normal_) | 1.0.19 |
| zstandard | `chunks.py` (training chunk cache) | 0.19.0 |
| matplotlib, Pillow | `scripts/visualize_tensor.py` + its tests | base matplotlib broken with NumPy 2.3 (README); Pillow 11.3.0 |
| rtg-tools (Java) | `pansoma_net_v2.vcfeval` | 3.12.1 / 3.13 (`$RTG`) |
| bcftools | `site_features build` (INDEL recipe) | module bcftools/1.18 |
| vg | giraffe, gamsort -i (`tools/jobs/gam_sort.sh`), stats, graph_prep (convert, chunk) | v1.77.0 (`scripts/use_vg.sh`); v1.65/1.66 for some steps (memory notes) |
| gbz_graph_index | graph index of a NEW graph (`tools/gbz_graph_index.cpp`) | built against gbwtgraph/sdsl in `/scratch/jshen/Github/gbz-tool/dependency` (not in repo); needs AVX/BMI, glibc ≥ 2.34 |
| GPU | training (bf16, torch.compile/Triton; H100) | inference runs on CPU too (the example does) |
| whole-genome resources | one node: ~48 CPUs, ~420–450 GiB for tasks; finalize ~22 GiB (v4 README) | — |

## 2. Files created (all under `$S`)

- `env_files/pansoma-cpu.yml`, `env_files/pansoma-gpu.yml`: the pinned envs, both verified. conda-forge/bioconda:
  python=3.11.3, libsqlite=3.41.2, cxx-compiler=2.0.0, bcftools=1.25, htslib=1.25, rtg-tools=3.13. pip:
  numpy 2.3.2, scipy 1.16.1, protobuf 5.29.3, pysam 0.23.0, zstandard 0.19.0, pybind11 3.0.1, torch 2.8.0,
  torchvision 0.23.0, timm 1.0.19, matplotlib 3.10.5, pillow 11.3.0. The torch wheels come from the cpu or cu128 index.
- `env_files/pansoma-{cpu,gpu}.conda-lock.txt` (`conda list --explicit --md5`) and `pansoma-{cpu,gpu}.pip-lock.txt`
  (`pip freeze`): full lock files of the verified envs. The GPU lock differs only by the nvidia-* CUDA 12.8 libraries and triton 3.4.0.
- `env_files/Dockerfile`, `env_files/pansoma.def` (micromamba + env file + vg v1.77.0 + clone + native compile).
  **Not built.** There is no Docker here, and `singularity build --fakeroot` fails ("no mapping entry found in /etc/subuid for jshen").
- `env_files/v1_conda_forge_stack/`: the first, unpinned env files. They were used only to show the failures in §4a.
- `patches/native_selftest_fallback.patch`: the `native.py` fix (§0.2).
- `example/run_example.sh`: the one command. It runs native compile if needed → `run discover` → `orchestrate prepare` +
  `run.sh` (build + merge, with the genome run's settings) → for SNV and INDEL: `predict` → `graph_vcf` → `linear_vcf` →
  `scripts/filter_panel_of_normals.py`. If `rtg` is on PATH it then runs `pansoma_net_v2.vcfeval` inside the benchmark BED.
- `example/check_example.py`: compares an output with `example/expected/expected.json`. Tensors: count plus SHA-256 of
  every shard and of the candidate ids, which must match exactly. VCF records must match exactly. P_SOMATIC must be within
  0.02, and a PASS/LowQual flip is allowed only within 0.02 of `##pansoma_threshold`.
- `example/expected/expected.json`: recorded from the clean CPU run (`example_out_cpu2`).
- `example/data/` (345 MB apparent, 31 files) and `example/data_SHA256SUMS.txt`:
  - the HG008-T Illumina (BCM 2024-03-13) records visiting the node IDs of GRCh38 chr1:242,498,000-243,002,000
    (nodes 3129540-3141988). That is 575,145 records, re-sorted with vg 1.77 `gamsort -i`: 67 MB GAM + GAI, and
    `vg stats -a` is identical before and after the sort.
  - graph index SQLite cut to 12,464 nodes (0.7 MB) and reference-path dir cut to those nodes (53 MB sparse .npy).
  - the chr table (unchanged).
  - GRCh38 no-alt chr1 with N outside 242.4-243.1 Mb (1.1 MB bgzip).
  - the 4 PoNs and the HG008-T somatic v0.2 truth (VCF + BED) cut to the window.
  - the 2 checkpoints `HG008_Illumina_{SNV,INDEL}_d9_normbq` (RIS `ris_7fc7396_20261012`, E8 "latest models") with the
    optimizer and scheduler state removed: 202 MB + 35 MB.
- `make_example/make_example.py` (how the data were cut; developer paths) and `make_example/compare_to_genome.py`
  (the example against the genome run).
- `jobs/verify.sh`, the clean-env runner. It uses `env -i HOME PATH=/usr/bin:/bin LANG`, then `conda activate <env>`,
  `PYTHONNOUSERSITE=1`, PYTHONPATH and PYTHONHOME unset, and a node-local TMPDIR. Also `jobs/create_env.sh`,
  `clean_env_tests.sh`, `run_example_clean.sh`, `make_example.sh`.

## 3. Exact commands (reviewer form)

```bash
git clone https://github.com/Jiawei-Shen/Pansoma.git && cd Pansoma            # after develop is merged into main
conda env create -p ./pansoma-env -f env/pansoma-cpu.yml                     # or env/pansoma-gpu.yml (driver >= 570)
conda activate ./pansoma-env
python -m indexed_gam_pipeline_v4.native compile                             # needs the native.py fix
python -m unittest discover -s indexed_gam_pipeline_v4/tests -t .
(cd machine_learning && python -m unittest discover -s pansoma_net_v2/tests -t .)
MPLBACKEND=Agg python -m unittest tests.test_tensor_visualization
# example data: download + unpack (Zenodo, to be deposited), then
bash example/run_example.sh pansoma_example_data out 4
python example/check_example.py out example/expected                         # "OK: matches the expected output"
```

How I created the envs here: `env -i HOME=$HOME PATH=/usr/bin:/bin CONDA_SOLVER=libmamba CONDA_PKGS_DIRS=$S/pkgs
PIP_CACHE_DIR=$S/pip_cache PYTHONNOUSERSITE=1 /wanglab/jshen/anaconda3/bin/conda env create -q -p $S/env_cpu2 -f
$S/env_files/pansoma-cpu.yml`. The package caches were fresh, so every package was downloaded. conda was 23.5.2 with libmamba.
How I verified: `ENVDIR=env_cpu2 REPODIR=repo_cpu2 TAG=cpu2 sbatch -w meow jobs/verify.sh` (all steps).
GPU: `sbatch -p gpu -w tequila --gres=gpu:12gb:1 ... STEPS=gpu`.
Final check: fresh clone `repo_final`, no compiled .so, example only, on miao.

## 4. Clean-environment results

Logs: `$S/logs/verify_cpu2.462888.out` (+ `logs/cpu2/*.log`), `$S/logs/verify_gpu2.462889.out` (+ `logs/gpu2/`),
`$S/logs/verify_final.462893.out` (+ `logs/final/`), and `$S/logs/compare_to_genome_cpu2.json`. Unpatched public code:
`$S/logs/clean_env_tests_cpu.462874.out` + `logs/test_*.log` and `$S/logs/run_example_cpu.462882.out`.

Every run: `sys.path` held only the env's paths, user site was False, and PYTHONPATH was unset. Hosts: meow and miao
(Xeon Gold 6336Y, Ubuntu 22.04, glibc 2.35) and tequila (H100 NVL, MIG 1g.12gb, driver 580.95.05).

**a. Public code as is (44821f3), unpinned conda-forge env (`env_cpu`, Python 3.11.18, SQLite 3.54, bioconda pysam): FAIL**
- native compile: exit 1 (the self-test of the static build segfaults).
- `run_example.sh`: stops at `discover` ("discover needs the native module").
- v4 tests: 207 run, 3 failures + 1 error, 18 skipped.
  - The golden failures are only "fixture inputs differ" for `graph.sqlite` / `*.gam` bytes. SQLite and the BGZF
    deflate library write different bytes; no tensor output differed.
  - The record-index error came from a rewritten GAI keeping its mtime: TMPDIR was on BeeGFS, which has 1-s mtimes.
- Python-decoder pass: OK, 207 tests (22 skipped). vg tests: 8 OK. ML: 112 OK (8 skipped). Visualizer: 14 OK.

**b. Pinned CPU env (`env_cpu2`) + native.py patch: PASS**
- env create: 247 s.
- native compile: 46 s (dynamic runtime, self-test 5,000 records).
- v4 tests: 207 OK (3 skipped: the native graph-index builder and 2 vg cases), 144 s. Every golden passes, inputs and outputs.
- Python-decoder pass: 207 OK (8 skipped), 131 s.
- `PANSOMA_VG=vg_v1.77.0` test_prep_tools: 8 OK.
- ML: 112 OK, 87 s. 8 skipped: 4 test_gpu, 4 test_real_data (`PANSOMA_TEST_TENSORS=/nonexistent`); the rtg case ran.
- Visualizer: 14 OK.
- Example: exit 0 in 166 s, peak RSS 5.5 GB, 4 threads.

**c. Pinned GPU env (`env_gpu2`) + patch, on tequila: PASS**
- env create: 300 s.
- ML: 112 OK (4 skipped, real data). The 4 test_gpu tests ran: torch.compile and bf16 on the H100 MIG slice.
- Example: 141 s.
- `check_example`: OK. Records and FILTERs are identical to the CPU run; max |ΔP_SOMATIC| is 0.0078 (SNV) and 0.0013 (INDEL).

**d. Final one-command run (fresh clone + patch + example files, `env_cpu2`, miao): PASS**
- `run_example.sh` compiled the decoder itself (49 s). Total 3 min 39 s, peak RSS 5.7 GB.
- `check_example`: "OK: matches the expected output". ΔP is 0.00000 against the meow run.

**Expected output** (`example/expected/expected.json`):
- SNV: 448 tensors, `*.SNV.pon.vcf.gz` with 66 records. FILTER: PASS 13, PanelOfNormals 15, LowQual 5, LowQual;PanelOfNormals 33.
- INDEL: 47 tensors, 31 records. FILTER: PASS 1, PanelOfNormals 12, LowQual 6, LowQual;PanelOfNormals 12.
- rtg vcfeval in the BED, after the PoN:
  - SNV: 13 calls, TP 12 / FP 1 / FN 0 (P 0.923, R 1.000, F1 0.960). Before the PoN: 28 calls, TP 12 / FP 16.
  - INDEL: 1 call, TP 1 / FP 0 / FN 5.
- First PASS records:
  `chr1 242525253 3130144:273:SNP:C>T C T QUAL 4.678 P_SOMATIC 0.65943`,
  `chr1 242528270 3130225:126:SNP:G>C ... 0.83121`, `chr1 242529195 3130238:65:SNP:G>T ... 0.66549`;
  INDEL PASS `chr1 242760124 3135065:71:DEL:A> GA G P_SOMATIC 0.48674`.
- The whole-file checksums are not portable, because the headers hold absolute paths. Use `check_example.py`.

**Representativeness** (`compare_to_genome.py`):
- All 448 SNV and 47 INDEL tensors are **byte-identical** to the HG008T_Illumina genome run's tensors at the same
  candidates, and their summaries are identical apart from the shard-layout fields.
- The window's VCF records match the genome run's chr1 calls exactly (SNV 65/65, INDEL 30/30), and so do the PASS sets
  before the PoN (28/28, 13/13).
- The genome run's vcfeval in the window gave the same result: SNV TP 12 / FP 1 / FN 0, INDEL TP 1 / FN 5.
- P_SOMATIC differs by ≤ 0.0105 (median 0.0009). The cause: predict fits the per-dataset BQ statistics on the example's
  own tensors, while the genome run used the checkpoint's training-dataset statistics.

**What is and is not covered.**
- Covered: from a sorted GAM + GAI to a PoN-filtered GRCh38 VCF and rtg vcfeval, on CPU and GPU. That is discovery,
  build, merge, predict, graph_vcf, linear_vcf and the PoN filter.
- Not covered:
  - vg giraffe alignment and gamsort of raw reads. The slice was re-sorted with vg 1.77, but alignment was not re-run.
  - graph preparation: graph_prep, gbz_graph_index, components, and the chr-index build.
  - truth labelling (`tensor_postprocessing label`, finalize with truth files) and training.
  - site features (the INDEL recipe's d46 final models).
  - whole-genome scale and memory.
  - the containers.

## 5. What a reviewer still needs

- Linux x86-64 and any conda/mamba with network access. The CPU env is 2.9 GB and the GPU env 7.1 GB.
- About 350 MB of example data. **It must first be deposited** (Zenodo) together with the two checkpoints.
- No GPU and no vg are needed for the example. A GPU (CUDA 12.8 wheels, driver ≥ 570) is needed for training and
  sensible for whole-genome inference. I did not measure whole-genome CPU throughput.
- For their own data:
  - the vg v1.77.0 static binary and the HPRC v1.1 d9 indexes (public).
  - Pansoma's graph files for that graph: graph-index SQLite 4.1 GB, reference path 1.8 GB, chr table. These are **not
    distributed**. Building them needs `gbz_graph_index`, whose gbwtgraph/sdsl dependencies come from the authors'
    `gbz-tool/dependency` and have no recipe in the repo.
  - a 48-CPU / ~450 GiB node per whole genome.

## 6. What must be committed/published for the response letter to be true

1. Merge `develop` into `main`, or make `develop` the GitHub default branch. Before that, decide what to do about
   `Reviewer_Response_Assessment.md` and the `analysis/reviewer_response_20261011/` drafts on the public branch.
2. `indexed_gam_pipeline_v4/native.py`: the patch in `$S/patches/native_selftest_fallback.patch`. Also adjust its docstring
   and v4 README §4 "native compile" ("static if it passes the self-test, else dynamic"). Without this, the published
   install path fails.
3. `env/pansoma-cpu.yml`, `env/pansoma-gpu.yml`, `env/*.conda-lock.txt`, `env/*.pip-lock.txt`. Add the
   Dockerfile / `pansoma.def` only after building and testing them, or leave them out of the letter.
4. `example/run_example.sh`, `example/check_example.py`, `example/expected/expected.json`, `example/data_SHA256SUMS.txt`,
   optionally `example/make_example.py` as provenance. The run_example usage line assumes `env/`.
5. Zenodo: the example data tarball plus the two stripped checkpoints (or the release models you choose), with a DOI or URL
   added to run_example.sh / README.
6. README:
   - an "Installation" section with the commands in §3.
   - a "Minimal example" section: download, the one command, expected output, and runtime (~3.5 min on 4 CPU threads,
     5.7 GB RAM; ~2.5 min on a GPU).
   - External tools: vg only for alignment and gamsort; gbz_graph_index only for a new graph.
   - Say that the `/scratch` and `/wanglab` paths in the READMEs (25 lines) and `scripts/use_vg.sh` are the authors' cluster.
   - Mark or remove `machine_learning/pansoma_net/{Dockerfile,environment.yml,requirements.txt}` as v1-only.
   - Update the pansoma_net_v2 README test count (it says 80; there are 112).
7. Optional test hardening: golden input hashes over SQLite/BGZF bytes (pass only with SQLite 3.41.2 + PyPI pysam,
   hence the pins), and `test_record_index`'s GAI rewrite on 1-s-mtime filesystems.

## 7. Draft response to Reviewer 2 (code comment; 116 words; true only after §6 items 1-6)

We thank the reviewer for testing Pansoma, and we apologise for the difficulty. The repository now provides pinned conda
environments (CPU and GPU, with lock files) and a minimal example. One command runs a sorted GAM of HG008-T Illumina reads
over GRCh38 chr1:242.5–243.0 Mb (a chromosome held out from training) through target discovery, tensor construction, the
SNV and INDEL models, VCF conversion and the panel-of-normals filter, and then compares the result with the expected
output. In freshly created environments without user site-packages, the example finished in 3.7 min on four CPU threads
(5.7 GB RAM) and in 2.3 min on an NVIDIA H100. The pipeline and model test suites passed. The example's tensors are
byte-identical to those of the genome-wide run.

(Precisely: both suites passed in the CPU environment, and the model suite, including its GPU tests, in the GPU
environment. The only cases skipped need the authors' full tensor sets or the GBZ index builder.)

## 8. Disk

`$S` is about 14 GB apparent: env_cpu (v1, kept as evidence; removable), env_cpu2, env_gpu2, example_out_*, and 4 clones.
I deleted my download caches (pkgs, pip_cache), the unused v1 GPU env and temp dirs, all by literal path.
