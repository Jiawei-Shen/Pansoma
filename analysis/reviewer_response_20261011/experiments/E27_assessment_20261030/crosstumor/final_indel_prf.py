"""Seed-mean precision / recall / F1 (ppr8, chr BED, after the shared PoN, calls with a GRCh38 position) of the final
INDEL batch, plus comparators, read from the rtg report.json files listed in results_ppr8.json. TSV on stdout."""
import json
D = json.load(open("/scratch/jshen/Github/Pansoma/tmp/indel_gap_20261020/final_20261023/results_ppr8.json"))
print("protocol\ttumor\tchrom\tplatform\tkind\tn_seeds\tP_mean\tR_mean\tF1_mean\tF1_range\tP_seeds\tR_seeds\tF1_seeds\ttruth\t"
      "DS_P\tDS_R\tDS_F1\tCL_P\tCL_R\tCL_F1\td_DS[CI]\td_CL[CI]")
for c in D["cells"]:
    P, R, F = [], [], []
    truth = None
    for s, v in sorted(c["seeds"].items()):
        a = json.load(open(v["report"]))["pon"]["at_threshold"]
        P.append(a["precision"]); R.append(a["recall"]); F.append(a["f1"])
        truth = a["tp_baseline"] + a["fn"]
    comp = []
    for k in ("deepsomatic", "clairs_to"):
        r = json.load(open(c["comparators"][k]["report"]))
        a = r["eval_bed"]["pass_pon"]["at_threshold"]
        comp += [a["precision"], a["recall"], a["f1"]]
    m = lambda x: sum(x) / len(x)
    f = lambda x: "/".join(f"{y:.3f}" for y in x)
    print("\t".join(map(str, [c["protocol"], c["tumor"], c["chrom"], c["platform"], c["kind"], len(F),
          f"{m(P):.3f}", f"{m(R):.3f}", f"{m(F):.3f}", f"{min(F):.3f}-{max(F):.3f}", f(P), f(R), f(F), truth] +
          [f"{x:.3f}" for x in comp] +
          [f"{c['d_deepsomatic']:+.3f} [{c['d_deepsomatic_ci'][0]:+.3f},{c['d_deepsomatic_ci'][1]:+.3f}]",
           f"{c['d_clairs_to']:+.3f} [{c['d_clairs_to_ci'][0]:+.3f},{c['d_clairs_to_ci'][1]:+.3f}]"])))
