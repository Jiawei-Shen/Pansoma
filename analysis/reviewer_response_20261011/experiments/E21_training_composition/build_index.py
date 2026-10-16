"""E21 step 1: build the pansoma_net_v2 index of one tensor set (SNV and INDEL) into this experiment's own cache.

The index (machine_learning/pansoma_net_v2/data.KindIndex, imported as a library, not modified) holds per tensor:
chromosome, node id, label, evaluation label, label reason, the site scalars (A1 AF is column 5), the GRCh38
position and, for tensors labelled 1, the somatic truth ids they stand for. Training replays its block split from
exactly these arrays, so summarize.py can count the training / validation tensors of every run without a GPU.

    python build_index.py SET [SET ...]        # SET = a directory name under /scratch/jshen/data/pansoma_v2_tensors

Reads only the tensor set (labels.ndjson, variant_summary.ndjson, labels npy, manifests, the BEDs named in
labels.manifest.json); writes only CACHE (default below, or $E21_CACHE).
"""
import os
import resource
import sys
import time
from pathlib import Path

REPO = Path("/scratch/jshen/Github/Pansoma")
sys.path.insert(0, str(REPO / "machine_learning"))
from pansoma_net_v2.data import KindIndex  # noqa: E402

TENSORS = Path("/scratch/jshen/data/pansoma_v2_tensors")
CACHE = Path(os.environ.get("E21_CACHE", "/scratch/jshen/data/reviewer_experiments_20261013/E21_training_composition/index_cache"))


def one(s, kind):
    d = TENSORS / s / "tensors" / kind
    if not (d / "labels.manifest.json").exists():
        return f"{s} {kind}: no labels.manifest.json, skipped"
    t0 = time.time()
    index = KindIndex(d, cache_dir=CACHE)
    rss = resource.getrusage(resource.RUSAGE_SELF).ru_maxrss / 1e6
    return f"{s} {kind}: {len(index):,} tensors, {time.time() - t0:.0f} s, max RSS {rss:.1f} GB"


def main(sets):
    """SNV and INDEL of each set in two processes (the Slurm job has -c 2)."""
    from concurrent.futures import ProcessPoolExecutor
    CACHE.mkdir(parents=True, exist_ok=True)
    with ProcessPoolExecutor(max_workers=2) as pool:
        for s in sets:
            for line in pool.map(one, [s, s], ["SNV", "INDEL"]):
                print(line, flush=True)


if __name__ == "__main__":
    main(sys.argv[1:])
