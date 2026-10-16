#!/bin/bash
# E22 check 2: read-only scan of the giraffe GAMs behind the COLO829T d9 tensor sets (and reference GAMs of known vg
# version) for the PARAMS_JSON group; compares the payloads with the test runs of vg_version_fingerprint.sh.
# Writes $W/gam_params/<name>.params.json and HERE/vg_fingerprint.tsv. Light: streams the first groups of each GAM.
set -uo pipefail
HERE=/scratch/jshen/Github/Pansoma/analysis/reviewer_response_20261011/experiments/E22_truth_and_version_checks
W=/scratch/jshen/data/reviewer_experiments_20261013/E22_truth_and_version_checks
PY=/wanglab/jshen/anaconda3/bin/python
mkdir -p $W/gam_params
# name | GAM | groups to scan
cat > $W/gam_params/gams.tsv <<'EOF'
COLO829T_Illumina_d9	/scratch/jshen/data/COLO829T/illumina/GAM/COLO829T.gam	3000
COLO829T_fiberseq_d9	/scratch/jshen/data/COLO829T/COLO829T_fiberseq/GAM/COLO829T_fiberseq.gam	100
COLO829T_ONT_d9	/scratch/jshen/data/COLO829T/COLO829T_ONT/COLO829T_ONT_merged.gam	100
COLO829BL_Illumina_d9_2025-06	/scratch/jshen/data/COLO829T/illumina/GAM/COLO829BL.gam	2000
SEQC2_SRR7890824_d9_2025-06	/scratch/jshen/data/SEQC2/SRR7890824_Tumor/SRR7890824.gam	2000
HG008T_Illumina_d9_vg1.66.0	/scratch/jshen/data/HG008_GIAB/vg_giraffe_GAM/Liss_lab_BCM_Illumina-WGS_20240313/HG008-T_p23_BCM_IlluminaWGS_20240313.gam	100
HG008T_PacBio_d9_vg1.66.0	/scratch/jshen/data/HG008_GIAB/raw_sequencing_data/Liss_lab_PacBio_Revio_20240125/HG008-T_PacBio-HiFi-Revio_20240125_116x.AF-HPRC.gam	100
HG008T_ONT_d9_vg1.66.0	/scratch/jshen/data/HG008_GIAB/raw_sequencing_data/Liss_lab_Northeastern-ONT-UL-20241216/HG008-T_GRCh38_GIABv3_ONT-UL-R10.4.1-dorado_0.8.1_sup.5mC_5hmC_54x_20241216.AF-HPRC.gam	100
COLO829_mix_ONT_p25_vg1.66.0	/scratch/jshen/data/COLO829T/Tumor_Normal_inSilico_Mix/COLO829_ONT/giraffe/MIX_p25_SMHTCOLO829_ont.gam	100
EOF
while IFS=$'\t' read -r name gam groups; do
  $PY -I $HERE/gam_tag_scan.py "$gam" --groups $groups > $W/gam_params/$name.scan.txt
  grep '^  payload:' $W/gam_params/$name.scan.txt | head -1 | sed 's/^  payload: //' > $W/gam_params/$name.params.json
done < $W/gam_params/gams.tsv
$PY -I - "$W" "$HERE/vg_fingerprint.tsv" <<'PY'
import hashlib, subprocess, sys
from pathlib import Path
W, out = Path(sys.argv[1]), Path(sys.argv[2])
SCANNER = str(out.parent / "gam_tag_scan.py")
def md5(p):
    b = p.read_bytes().strip() if p.exists() else b""
    return hashlib.md5(b + b"\n").hexdigest()[:12] if b else "none"
rows = []
for line in (W / "gam_params/gams.tsv").read_text().splitlines():
    name, gam, groups = line.split("\t")
    scan = (W / f"gam_params/{name}.scan.txt").read_text().strip().splitlines()[-1]
    rows.append(("GAM", name, scan.split(": ", 1)[1], md5(W / f"gam_params/{name}.params.json")))
fp = W / "vg_fingerprint"
for d in sorted(p for p in fp.iterdir() if p.is_dir()):
    for mode in ("short", "short_oldmin", "hifi", "r10"):
        g = d / f"{mode}.gam"
        if not g.exists() or g.stat().st_size == 0:
            continue
        scan = subprocess.run([sys.executable, "-I", SCANNER, str(g), "--groups", "1000"],
                              capture_output=True, text=True).stdout.strip().splitlines()
        payload = [l[len("  payload: "):] for l in scan if l.startswith("  payload: ")]
        h = hashlib.md5((payload[0].strip() + "\n").encode()).hexdigest()[:12] if payload else "none"
        rows.append(("test", f"{d.name} {mode}", scan[-1].split(": ", 1)[1], h))
with open(out, "w") as fh:
    fh.write("kind\tname\tscan (groups, alignments, tags)\tPARAMS_JSON md5\n")
    for r in rows:
        fh.write("\t".join(r) + "\n")
print(out.read_text())
PY
