#!/usr/bin/env python3
"""CLI for the shared tensor image visualizer: legacy five-channel, indexed_gam_pipeline candidate and candidate-v2..v6 tensors."""
from pathlib import Path
import sys

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))
from pangenome_ml_data_generation.tensors.visualization import main

if __name__ == "__main__":
    main()
