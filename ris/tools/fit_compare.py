"""Compare rescorer fits on the tuning chromosome (chr2): per set the block-CV F1 (rescorer and p-only) and the
leave-one-sample-out F1, with means. Standard library only (login node).

    python3 ris/tools/fit_compare.py RESCORE_RUN [RESCORE_RUN ...]      (names under pansoma_net_v2_runs/ or paths)
"""
import json
import os
import sys

RUNS = "/storage2/fs1/epigenome/Active/jiawei.s/pansoma_net_v2_runs"
ORDER = ["HG008T_Illumina", "HG008T_PacBio", "HG008T_ONT", "COLO829T_Illumina", "COLO829T_fiberseq", "COLO829T_ONT"]


def main():
    print(f"{'run':18s} {'kind':9s} " + " ".join(f"{s.replace('HG008T_', 'H.').replace('COLO829T_', 'C.'):>10s}" for s in ORDER)
          + f" {'mean':>7s} {'min':>7s}")
    for r in sys.argv[1:]:
        path = r if os.path.isabs(r) else f"{RUNS}/{r}"
        m = json.load(open(f"{path}/model.json"))
        cv = m["chr2_block_cv"]
        loso = {}
        for name, v in m["variants"].items():
            if name.startswith("loso_"):
                loso.update(v["chr2_apply_f1"])
        rows = [("CV p-only", {s: cv[s]["p_only"] for s in cv}), ("CV rescor", {s: cv[s]["rescorer"] for s in cv}),
                ("LOSO", loso)]
        for kind, vals in rows:
            xs = [vals.get(s) for s in ORDER]
            ok = [x for x in xs if x is not None]
            print(f"{os.path.basename(path):18s} {kind:9s} " + " ".join(f"{x:10.3f}" if x is not None else f"{'-':>10s}" for x in xs)
                  + (f" {sum(ok) / len(ok):7.4f} {min(ok):7.3f}" if ok else ""))
        print(f"{'':18s} feature set {m.get('feature_set', 'v1')}, pooled floor {m['variants']['pooled']['floor']}, "
              f"threshold {m['variants']['pooled']['threshold']:.4f}")


if __name__ == "__main__":
    main()
