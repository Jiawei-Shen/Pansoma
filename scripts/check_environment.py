#!/usr/bin/env python3
"""Check local dependencies needed by the Pansoma pipeline."""

from __future__ import annotations

import importlib
import argparse
import os
import platform
import shutil
import subprocess
import sys
from pathlib import Path


PYTHON_MODULES = [
    "numpy",
    "pandas",
    "matplotlib",
    "pysam",
    "pybind11",
    "google.protobuf",
]

ML_MODULES = [
    "torch",
    "torchvision",
    "timm",
    "tqdm",
]

COMMANDS = [
    ("vg", "required for FASTQ -> GAM and graph extraction"),
    ("jq", "required by graph/node shell utilities"),
    ("samtools", "recommended for BAM/VCF-adjacent workflows"),
    ("bcftools", "required to merge per-chromosome inference VCFs"),
    ("tabix", "required to index truth and output VCFs"),
]

PROJECT_ROOT = Path(__file__).resolve().parents[1]
SRC_ROOT = PROJECT_ROOT / "src"
for path in (PROJECT_ROOT, SRC_ROOT):
    path_text = str(path)
    if path_text not in sys.path:
        sys.path.insert(0, path_text)


def check_module(name: str) -> bool:
    try:
        module = importlib.import_module(name)
    except Exception as exc:
        print(f"[missing] python module {name}: {exc}")
        return False
    version = getattr(module, "__version__", "installed")
    print(f"[ok] python module {name}: {version}")
    return True


def check_command(command: str, note: str) -> bool:
    path = shutil.which(command)
    if not path:
        print(f"[missing] command {command}: {note}")
        return False
    try:
        result = subprocess.run(
            [command, "--version"],
            check=False,
            capture_output=True,
            text=True,
            timeout=10,
        )
        first_line = (result.stdout or result.stderr).strip().splitlines()
        version = first_line[0] if first_line else "installed"
    except Exception:
        version = "installed"
    print(f"[ok] command {command}: {path} ({version})")
    return True


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--strict-external",
        action="store_true",
        help="fail if external command-line tools such as vg are missing",
    )
    parser.add_argument(
        "--ml",
        action="store_true",
        help="also validate PyTorch and PansomaNet dependencies",
    )
    return parser.parse_args()


def main() -> int:
    args = parse_args()
    print(f"Platform: {platform.system()} {platform.machine()}")
    print(f"Python: {sys.executable}")
    print(f"Conda environment: {os.environ.get('CONDA_DEFAULT_ENV', 'not active')}")
    print()

    ok = True
    for module in PYTHON_MODULES:
        ok = check_module(module) and ok

    if args.ml:
        for module in ML_MODULES:
            ok = check_module(module) and ok

    try:
        import google.protobuf

        protobuf_version = google.protobuf.__version__
    except Exception:
        pass
    else:
        if protobuf_version != "3.20.3":
            print(
                "[incompatible] protobuf: expected 3.20.3 for the bundled "
                f"vg_pb2.py, found {protobuf_version}"
            )
            ok = False

    print()
    for command, note in COMMANDS:
        command_ok = check_command(command, note)
        if args.strict_external:
            ok = command_ok and ok

    print()
    try:
        import fast_writer  # noqa: F401
    except Exception as exc:
        print(f"[missing] fast_writer extension: {exc}")
        print("          run: bash scripts/build_fast_writer.sh")
        ok = False
    else:
        print("[ok] fast_writer extension")

    return 0 if ok else 1


if __name__ == "__main__":
    raise SystemExit(main())
