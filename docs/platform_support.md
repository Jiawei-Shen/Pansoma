# Platform Support

Pansoma is intended to run on Linux/HPC systems and macOS development machines.
The full production pipeline is easiest on Linux because `vg` release binaries
and large graph indexes are usually managed there.

## Supported Roles

| Platform | Recommended role | Notes |
| --- | --- | --- |
| Linux x86_64 / HPC | Full pipeline | Best target for `vg giraffe`, SLURM jobs, large graph indexes, tensor generation, and model inference. |
| macOS arm64 / Apple Silicon | Development and small tests | Python utilities and the `fast_writer` extension build locally. Full FASTQ-to-GAM requires a working `vg` build plus graph indexes. |
| macOS x86_64 | Development and small tests | Similar to Apple Silicon; `vg` is generally built from source. |

## Python Environment

Use the shared Conda environment on both Linux and macOS:

```bash
conda env create -f environment.yml
conda activate pangenome-ml-data-generation
```

The environment pins `protobuf==3.20.3` because the checked-in `vg_pb2.py`
was generated with an older protobuf toolchain. Regenerate `vg_pb2.py` before
upgrading protobuf.

## Build The C++ Extension

```bash
bash scripts/build_fast_writer.sh
```

The extension uses `mmap`/`msync` on both Linux and macOS. On Linux it also
uses `posix_fadvise(..., POSIX_FADV_DONTNEED)` when available to reduce page
cache pressure. macOS does not expose that constant, so the optimization is
compiled out there.

## Validate A Machine

```bash
python scripts/check_environment.py --ml
```

This checks Python and ML packages, command-line tools, and the compiled
`fast_writer` extension. By default, missing external tools are reported as
warnings so macOS development machines can still pass core validation. For a
full Linux/HPC pipeline node, use:

```bash
python scripts/check_environment.py --strict-external --ml
```

## External Tools

`jq` and `samtools` are available through Conda/Homebrew on both Linux and
macOS. `vg` is different:

- Linux: use the official `vg` Linux release binary or a cluster module.
- macOS: build `vg` from source following the upstream macOS instructions.

The full FASTQ-to-GAM stage also requires graph indexes:

```text
*.gbz
*.min
*.dist
```

Those files are intentionally not stored in this repository.
