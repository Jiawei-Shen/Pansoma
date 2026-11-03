"""E8's score_comparator.py on RIS: the same protocol and code, the wangcluster paths swapped for their RIS copies.

    python score_comparator_ris.py --vcf snv.vcf.gz --vcf indel.vcf.gz --set COLO829T_Illumina_P75 --kind INDEL \
        --caller clairs_to --output DIR [--threads 4]

../../reviewer_response_20261011/experiments/E8_comparators/score_comparator.py is imported unchanged; only its module
paths are replaced: tensor sets ($PANSOMA_DATA/pansoma_v2_tensors), the four PoN VCFs ($PANSOMA_DATA/panel_of_normal_VCFs,
the files ris/vcf.sh uses), the GRCh38 no-alt FASTA / SDF of ris/vcf.sh and rtg ($RTG, ris/env.sh). The label manifests'
truth VCF / BED keep their old-cluster paths, bind-mounted by ris/submit.sh. bcftools comes from PATH (score.sh).
"""
import importlib.util
import os
import sys
from pathlib import Path

REPO = Path(__file__).resolve().parents[3]
sys.path.insert(0, str(REPO / "machine_learning"))
spec = importlib.util.spec_from_file_location(
    "score_comparator", REPO / "analysis/reviewer_response_20261011/experiments/E8_comparators/score_comparator.py")
sc = importlib.util.module_from_spec(spec)
spec.loader.exec_module(sc)

DATA = Path(os.environ.get("PANSOMA_DATA", "/storage2/fs1/epigenome/Active/jiawei.s/Dataset"))
sc.TENSORS = DATA / "pansoma_v2_tensors"
sc.PON_DIR = DATA / "panel_of_normal_VCFs"
sc.PONS = [sc.PON_DIR / p.name for p in sc.PONS]
sc.FASTA = str(DATA / "REF/GCA_000001405.15_GRCh38_no_alt_analysis_set.fa")
sc.SDF = str(DATA / "REF/GCA_000001405.15_GRCh38_no_alt_analysis_set.sdf")
sc.RTG = os.environ["RTG"]
_labels_manifest = sc.labels_manifest
sc.labels_manifest = lambda set_name, kind: (*_labels_manifest(set_name, kind)[:2], sc.FASTA)  # the labels' FASTA, RIS copy

if __name__ == "__main__":
    sc.evaluate(sc.parse_args())
