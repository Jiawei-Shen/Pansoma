"""E13 exclude files on RIS: analysis/reviewer_response_20261011/experiments/E21_training_composition/make_e13_subsets.py
run unchanged on the RIS copy of the set, with summarize.py's tequila paths pointed at RIS (tensors and the index cache
the training runs use, so the cache is read, not rebuilt). The subsets are SHA-1 draws over the seed-0 block split, so
they equal tequila's when the set does: e13_subsets_<SET>_<KIND>.tsv here must match the E21 one (all columns but
`file`). Same arguments as make_e13_subsets.py; e13_subsets.sh runs it in the container.
"""
import os
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
REPO = HERE.parents[2]
E21 = REPO / "analysis/reviewer_response_20261011/experiments/E21_training_composition"
sys.path.insert(0, str(REPO / "machine_learning"))
sys.path.insert(0, str(E21))

import summarize  # noqa: E402

summarize.TENSORS = Path(os.environ["PANSOMA_DATA"]) / "pansoma_v2_tensors"
summarize.CACHE = Path(os.environ["PANSOMA_RUNS"]) / "index_cache"

import make_e13_subsets  # noqa: E402  (imports TENSORS / CACHE from summarize: the RIS paths)

make_e13_subsets.HERE = HERE  # the subset table goes here, not over E21's

if __name__ == "__main__":
    make_e13_subsets.main()
