"""Before vs after the haplotype-overlap fix (c9eb3af): label counts per reason for the six sets.
'Before' = the backup that relabel.sh (jobs 376732-376737) wrote of the truth-labels-v6 manifests."""
import glob, json, os
R = "/scratch/jshen/data/pansoma_v2_tensors"
JOBS = dict(l.split() for l in open("relabel_ovfix_jobs.txt"))
SETS = [("HG008 PacBio", "Liss_lab_PacBio_Revio_20240125", "HG008_PacBio"), ("HG008 ONT", "Liss_lab_Northeastern-ONT-UL-20241216", "HG008_ONT"),
        ("HG008 Illumina", "Liss_lab_BCM_Illumina-WGS_20240313", "HG008_Illumina"), ("COLO829T Illumina", "COLO829T_Illumina", "COLO_Illumina"),
        ("COLO829T ONT", "COLO829T_ONT", "COLO_ONT"), ("COLO829T fiberseq", "COLO829T_fiberseq", "COLO_fiberseq")]
for name, d, key in SETS:
    (backup,) = glob.glob(f"{R}/{d}/labels_backup_v3_tensors_*_{JOBS[key]}")
    for kind in ("SNV", "INDEL"):
        old = json.load(open(f"{backup}/{kind}.labels.manifest.json")); new = json.load(open(f"{R}/{d}/v3_tensors/{kind}/labels.manifest.json"))
        t0, t1 = old["totals"], new["totals"]
        print(f"{name} {kind}: " + "  ".join(f"{v}: {t0.get(k,0):,} -> {t1.get(k,0):,}" for k, v in (("somatic", 1), ("germline", 2), ("non", 0), ("ignore", -1))))
        r0, r1 = old["reasons"], new["reasons"]
        for r in sorted(set(r0) | set(r1)):
            if r0.get(r, 0) != r1.get(r, 0):
                print(f"      {r:48s} {r0.get(r,0):>9,} -> {r1.get(r,0):>9,}")
