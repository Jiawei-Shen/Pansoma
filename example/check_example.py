"""Compare an example run (run_example.sh OUT) with the expected output.

    python example/check_example.py OUT EXPECTED [--write]

EXPECTED/expected.json holds, per kind (SNV, INDEL):
  tensors      the number of merged tensors and the SHA-256 of every <chrom>_shard_*_data.npy and of the summary
               records' candidate ids: the pipeline is integer-exact, so these must match exactly
  calls        every record of <name>.pon.vcf.gz: CHROM POS ID REF ALT FILTER and P_SOMATIC
Records must agree in CHROM POS ID REF ALT; P_SOMATIC within --tolerance (default 0.02: CPU / GPU and BLAS
differences move the model's probabilities by far less); FILTER may differ only for a record whose P_SOMATIC is within
the tolerance of the PASS threshold (##pansoma_threshold). --write records OUT as the expected output. Exit 1 on
any difference.
"""
import argparse
import gzip
import hashlib
import json
from pathlib import Path
import sys

KINDS = ("SNV", "INDEL")


def sha256(path):
    h = hashlib.sha256()
    with open(path, "rb") as stream:
        for block in iter(lambda: stream.read(1 << 22), b""):
            h.update(block)
    return h.hexdigest()


def tensors(out, kind):
    d = Path(out) / "HG008T_Illumina" / "tensors" / kind
    shards = sorted(d.glob("*_shard_*_data.npy"))
    ids = []
    for summary in sorted(d.glob("*_variant_summary.ndjson")):
        with open(summary) as stream:
            ids += [json.loads(line)["candidate_id"] for line in stream]
    return dict(count=len(ids), shards={p.name: sha256(p) for p in shards},
                candidate_ids_sha256=hashlib.sha256("\n".join(ids).encode()).hexdigest())


def calls(out, kind):
    path = Path(out) / "vcf" / f"HG008T_Illumina.tensors.{kind}.pon.vcf.gz"
    threshold, records = None, []
    with gzip.open(path, "rt") as stream:
        for line in stream:
            if line.startswith("##pansoma_threshold="):
                threshold = float(line.split("=", 1)[1])
            if line.startswith("#"):
                continue
            f = line.rstrip("\n").split("\t")
            info = dict(x.split("=", 1) for x in f[7].split(";") if "=" in x)
            records.append(dict(chrom=f[0], pos=int(f[1]), id=f[2], ref=f[3], alt=f[4], filter=f[6],
                                p_somatic=float(info["P_SOMATIC"])))
    return dict(threshold=threshold, records=records,
                filters={k: sum(r["filter"] == k for r in records) for k in sorted({r["filter"] for r in records})})


def summary(out):
    return {kind: dict(tensors=tensors(out, kind), calls=calls(out, kind)) for kind in KINDS}


def compare(got, want, tolerance):
    problems = []
    for kind in KINDS:
        g, w = got[kind], want[kind]
        if g["tensors"] != w["tensors"]:
            problems.append(f"{kind}: tensors differ ({g['tensors']['count']} vs expected {w['tensors']['count']}, "
                            f"shards or candidate ids not identical)")
        gr, wr = g["calls"]["records"], w["calls"]["records"]
        key = lambda r: (r["chrom"], r["pos"], r["id"], r["ref"], r["alt"])  # noqa: E731
        if [key(r) for r in gr] != [key(r) for r in wr]:
            problems.append(f"{kind}: VCF records differ ({len(gr)} vs expected {len(wr)})")
            continue
        worst = 0.0
        for a, b in zip(gr, wr):
            d = abs(a["p_somatic"] - b["p_somatic"])
            worst = max(worst, d)
            if d > tolerance:
                problems.append(f"{kind} {key(a)}: P_SOMATIC {a['p_somatic']} vs expected {b['p_somatic']}")
            elif a["filter"] != b["filter"]:
                near = {a["filter"], b["filter"]} <= {"PASS", "LowQual"} and \
                    abs(b["p_somatic"] - w["calls"]["threshold"]) <= tolerance
                if not near:  # a PASS / LowQual flip is expected only next to the threshold
                    problems.append(f"{kind} {key(a)}: FILTER {a['filter']} vs expected {b['filter']}")
                else:
                    print(f"note: {kind} {key(a)} {b['filter']} -> {a['filter']} (P_SOMATIC {a['p_somatic']}, "
                          f"threshold {w['calls']['threshold']:.4f})")
        print(f"{kind}: {g['tensors']['count']} tensors (identical: {g['tensors'] == w['tensors']}); "
              f"{len(gr)} VCF records, FILTER {g['calls']['filters']} (expected {w['calls']['filters']}); "
              f"max |P_SOMATIC - expected| {worst:.5f}")
    return problems


def main():
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("out")
    p.add_argument("expected")
    p.add_argument("--tolerance", type=float, default=0.02)
    p.add_argument("--write", action="store_true", help="record OUT as the expected output")
    args = p.parse_args()
    got = summary(args.out)
    target = Path(args.expected) / "expected.json"
    if args.write:
        target.parent.mkdir(parents=True, exist_ok=True)
        target.write_text(json.dumps(got, indent=1) + "\n")
        print(f"wrote {target}")
        return 0
    problems = compare(got, json.loads(target.read_text()), args.tolerance)
    for line in problems:
        print("DIFFERENCE", line)
    print("OK: matches the expected output" if not problems else f"FAILED: {len(problems)} differences")
    return 1 if problems else 0


if __name__ == "__main__":
    sys.exit(main())
