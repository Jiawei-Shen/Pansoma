#!/usr/bin/env python3
"""CLI for the shared legacy / candidate-v2..v5 tensor image visualizer."""
from pathlib import Path
import sys

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))
from pangenome_ml_data_generation.tensors.visualization import main

if __name__ == "__main__":
    main()
