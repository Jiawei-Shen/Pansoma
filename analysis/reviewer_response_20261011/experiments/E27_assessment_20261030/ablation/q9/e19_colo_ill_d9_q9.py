"""q9 (assessment 2026-10-30): E19 germline-density strata (experiments/E19_germline_density/germline_density.py, read
only, run unchanged) for COLO829T Illumina SNV with an HPRC v1.1 (d9) arm, which E19 never had. One variant per call:
  ep8  d9 = reviewer_cpu_20261013/COLO829_Illumina_SNV_7fc7396_ddp4_b1024 (validation-selected epoch 8, seed 0)
  ep6  d9 = reviewer_cpu_20261013/COLO829_Illumina_SNV_d9_normbq_sweep_e06 (chr1-chosen epoch 6, seed 0)
with the linear / d46 seed-0 runs E19 used. Per-record caches and tables go to q9/ only (calls_lib.DATA and E19's OUT /
HERE are redirected; the chr1 tensor-key caches are copied there beforehand with cp -p).
usage: python e19_colo_ill_d9_q9.py ep8|ep6
"""
import sys
from pathlib import Path

Q9 = Path(__file__).resolve().parent
EXP = Path("/scratch/jshen/Github/Pansoma/analysis/reviewer_response_20261011/experiments")
sys.path.insert(0, str(EXP / "call_level_common"))
sys.path.insert(0, str(EXP / "E19_germline_density"))
import calls_lib as L  # noqa: E402
import germline_density as G  # noqa: E402

D9 = {"ep8": "COLO829_Illumina_SNV_7fc7396_ddp4_b1024", "ep6": "COLO829_Illumina_SNV_d9_normbq_sweep_e06"}
KEEP = {"COLO829_Illumina_SNV_linear_7fc7396", "COLO829_Illumina_SNV_d46_7fc7396"}


def main():
    which = sys.argv[1]
    L.DATA = Q9 / "data"
    G.OUT = Q9 / "data" / "E19"
    G.HERE = Q9 / f"E19_{which}"
    G.HERE.mkdir(parents=True, exist_ok=True)
    orig = L.discover
    keep = KEEP | {D9[which]}
    L.discover = lambda **kw: [cs for cs in orig(groups=("insample",), variants=("main",)) if cs.run in keep and cs.kind == "SNV"]
    G.main()


if __name__ == "__main__":
    main()
