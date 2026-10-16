"""E1.1 extension: the assembly off-reference somatic SNV keys of an HPRC v2.1 d46 set (asm_d46_classify.py's
offref_somatic_d46.tsv) scored by a d46 model under the E1.1 modes, with the LOW and training mean / median of that
run's main E1.1 report (ch6_counterfactual.py), so the edits are the same as on the run's test tensors. The shuffle
donors are random chr1 test tensors of the set. Reads only; writes OUT/report_assembly.json.

usage: python asm_scores.py --run RUN --set SET --table offref_somatic_d46.tsv --main OUT_OF_MAIN_RUN [--threads 4]"""
import argparse
import csv
import json
import sys
import time
from pathlib import Path

import numpy as np
import torch

sys.path.insert(0, str(Path(__file__).resolve().parent))
import ch6_counterfactual as cf  # noqa: E402
from pansoma_net_v2.model import PansomaNetV2  # noqa: E402


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("--run", required=True)
    ap.add_argument("--set", required=True)
    ap.add_argument("--table", required=True)
    ap.add_argument("--main", required=True, help="the run's ch6_counterfactual.py output directory (report.json)")
    ap.add_argument("--threads", type=int, default=4)
    args = ap.parse_args()
    torch.set_num_threads(args.threads)
    main_dir = Path(args.main)
    rep = json.loads((main_dir / "report.json").read_text())
    params = {k: rep["params"][k] for k in ("low", "train_mean", "train_median")}
    modes = list(rep["results"])
    ck = torch.load(cf.RUNS / args.run / "best.pth", map_location="cpu", weights_only=False)
    model = PansomaNetV2.from_checkpoint(ck).to(cf.DEVICE).eval()
    thr = float(ck["val"]["somatic_at_recall"]["0.9"]["threshold"])
    assert abs(thr - rep["threshold_val_recall_0_9"]) < 1e-9
    index = cf.open_index(cf.TENSORS / args.set / "tensors" / "SNV", main_dir.parent / "index_cache")
    a, reasons, cands, chroms = index.arrays, index.meta["reasons"], index.candidates(), index.meta["chroms"]
    rows = [r for r in csv.DictReader(open(args.table), delimiter="\t") if r["kind"] == "SNP"]
    by_key = {r["key"]: r for r in rows}
    found = sorted(p for p, c in enumerate(cands) if c in by_key)
    asm = np.array(found, np.int64)
    ignore = set(ck.get("args", {}).get("ignore_reasons") or [])
    lab = a["eval_label"].astype(np.int64).copy()
    lab[np.isin(a["reason"], [k for k, r in enumerate(reasons) if r in ignore])] = -1
    pos = index.select(["chr1"], labelled=False)
    test_pos = pos[lab[pos] >= 0]
    rng = np.random.default_rng(1)
    donor_pos = rng.choice(test_pos, len(asm), replace=False)
    ev = np.concatenate([asm, donor_pos])
    donors = np.arange(len(asm), 2 * len(asm))  # tensor k takes the ch6 of donor k
    cf.log(f"{args.run}: {len(rows)} assembly SNP keys, {len(asm)} with a tensor; params {params}; modes {modes}")
    out_rows = [dict(key=cands[p], event=by_key[cands[p]]["event"], grch38=by_key[cands[p]]["grch38_records"],
                     chrom=chroms[a["chrom"][p]], label=int(a["label"][p]), reason=reasons[a["reason"][p]],
                     off_reference=bool(a["off_reference"][p]), in_test=bool(lab[p] >= 0 and chroms[a["chrom"][p]] == "chr1"))
                for p in asm]
    summary = {}
    for mode in modes:
        ds = cf.Counterfactual([(index, ev)], mode, params, np.concatenate([donors, donors]))
        probs = []
        t0 = time.time()
        with torch.no_grad():
            for s in range(0, len(asm), 64):
                items = [ds[i] for i in range(s, min(s + 64, len(asm)))]
                x = torch.stack([t[0] for t in items])
                b = torch.stack([t[1] for t in items])
                sc = torch.stack([t[2] for t in items])
                with torch.autocast(cf.DEVICE.type, dtype=torch.bfloat16, enabled=cf.DEVICE.type == "cuda"):
                    logits = model(x.to(cf.DEVICE), b.to(cf.DEVICE), sc.to(cf.DEVICE))
                probs.append(torch.softmax(logits.float(), 1).cpu().numpy())
        probs = np.concatenate(probs)
        for r, p in zip(out_rows, probs):
            r[f"p_{mode}"] = round(float(p[1]), 4)
            r[f"argmax_{mode}"] = int(p.argmax())
        chr1 = np.array([r["chrom"] == "chr1" for r in out_rows])
        summary[mode] = dict(n=int(len(asm)), ge_thr=int((probs[:, 1] >= thr).sum()),
                             argmax_somatic=int((probs.argmax(1) == 1).sum()), p_median=float(np.median(probs[:, 1])),
                             chr1_n=int(chr1.sum()), chr1_ge_thr=int((probs[chr1, 1] >= thr).sum()),
                             chr1_argmax_somatic=int((probs[chr1].argmax(1) == 1).sum()))
        cf.log(f"  {mode}: {summary[mode]} ({time.time() - t0:.0f}s)")
    (main_dir / "report_assembly.json").write_text(json.dumps(dict(
        run=args.run, set=args.set, table=args.table, threshold_val_recall_0_9=thr, params=params,
        assembly_keys_in_table=len(rows), summary=summary, assembly_keys=out_rows), indent=1) + "\n")


if __name__ == "__main__":
    main()
