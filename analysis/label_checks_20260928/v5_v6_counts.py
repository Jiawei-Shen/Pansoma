"""truth-labels-v5 vs v6: label counts per reason for the six sets (v5 = the backup manifest the relabel job wrote)."""
import glob, json, os, collections
R = "/scratch/jshen/data/pansoma_v2_tensors"
SETS = [("HG008 PacBio", "Liss_lab_PacBio_Revio_20240125"), ("HG008 ONT", "Liss_lab_Northeastern-ONT-UL-20241216"),
        ("HG008 Illumina", "Liss_lab_BCM_Illumina-WGS_20240313"), ("COLO829T Illumina", "COLO829T_Illumina"),
        ("COLO829T ONT", "COLO829T_ONT"), ("COLO829T fiberseq", "COLO829T_fiberseq")]
VALUE = {"somatic": 1, "germline": 2, "non": 0, "ignore": -1}
for name, d in SETS:
    backups = sorted(glob.glob(f"{R}/{d}/labels_backup_v3_tensors_truth-labels-v5_*"), key=os.path.getmtime)
    for kind in ("SNV", "INDEL"):
        new = json.load(open(f"{R}/{d}/v3_tensors/{kind}/labels.manifest.json"))
        old = json.load(open(f"{backups[-1]}/{kind}.labels.manifest.json")) if backups else None
        tot_old = old["totals"] if old else {}; tot_new = new["totals"]
        print(f"\n{name} {kind}: {new['version']} (min_overlap {new.get('min_overlap')}, snv_min_af {new.get('snv_min_af')})  tensors {new['tensors']:,}")
        print("   " + "  ".join(f"{VALUE[k]:>2}: {tot_old.get(k,0):>9,} -> {tot_new.get(k,0):>9,}" for k in ("somatic", "germline", "non", "ignore")))
        ro, rn = (old or {}).get("reasons", {}), new["reasons"]
        for r in sorted(set(ro) | set(rn)):
            if ro.get(r, 0) != rn.get(r, 0):
                print(f"      {r:52s} {ro.get(r,0):>9,} -> {rn.get(r,0):>9,}")
