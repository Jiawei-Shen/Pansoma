"""Training / validation tensor counts of the 18 final INDEL runs (release / unseen x 3 platforms x 3 seeds),
parsed from each run's train.log. Output TSV on stdout."""
import re, sys, glob, os
R = "/scratch/jshen/data/pansoma_net_v2_runs"
num = lambda s: int(s.replace(",", ""))
print("\t".join(["run", "protocol", "platform", "seed", "set", "train", "train_non", "train_somatic", "train_germline",
                 "val", "val_non", "val_somatic", "val_germline", "left_out_train", "left_out_val", "val_offref_as_non",
                 "val_truth_alleles", "val_truth_with_tensor", "best_epoch"]))
for proto in ("unseen", "release"):
    for p in ("Illumina", "PacBio", "ONT"):
        for s in (0, 1, 2):
            run = f"{proto}_s{s}_{p}_INDEL"
            t = open(f"{R}/{run}/train.log").read()
            for m in re.finditer(r"tensors/(\S+)_d46/tensors/INDEL: labels .*?train ([\d,]+) non=([\d,]+) somatic=([\d,]+) germline=([\d,]+)", t):
                pass
            for m in re.finditer(r"pansoma_v2_tensors/(\S+?)/tensors/INDEL: labels [^\n]*train ([\d,]+) non=([\d,]+) somatic=([\d,]+) germline=([\d,]+)", t):
                print("\t".join([run, proto, p, str(s), m.group(1)] + [str(num(m.group(i))) for i in range(2, 6)] + [""] * 10))
            m = re.search(r"training tensors ([\d,]+) \(\{'non': (\d+), 'somatic': (\d+), 'germline': (\d+)\}\), validation ([\d,]+) \(\{'non': (\d+), 'somatic': (\d+), 'germline': (\d+)\}\)", t)
            lo = re.search(r"left out \(--ignore-reasons[^)]*\): ([\d,]+) training, ([\d,]+) validation", t)
            off = re.search(r"validation includes ([\d,]+) off-reference", t)
            vt = re.search(r"truth-level validation: ([\d,]+) somatic truth alleles in the validation region \(([\d,]+) with", t)
            be = re.search(r"done; best epoch (\d+)", t)
            print("\t".join([run, proto, p, str(s), "ALL"] + [str(num(m.group(i))) for i in range(1, 9)] +
                            [str(num(lo.group(1))), str(num(lo.group(2))), str(num(off.group(1))),
                             str(num(vt.group(1))), str(num(vt.group(2))), be.group(1)]))
