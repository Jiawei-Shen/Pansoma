"""Developer check: the example's tensors and calls against the HG008T Illumina genome run in the same window.

    python compare_to_genome.py EXAMPLE_OUT

Tensors: every example candidate on a node of the extracted node-ID range has all of its GAM records in the
example GAM, so its tensor bytes and summary must equal the genome run's (pansoma_v2_tensors/HG008T_Illumina).
Calls: the genome run's chr1 linear VCF (RIS runs HG008_Illumina_{SNV,INDEL}_d9_normbq, vcf_chr1/) inside the
window against the example's. The genome run z-scored BQ with the checkpoint's statistics of its training dataset;
the example fits them on its own tensors (predict, --norm dataset), so probabilities differ a little.
"""
import gzip
import json
from pathlib import Path
import sys

import numpy as np

GENOME = Path("/scratch/jshen/data/pansoma_v2_tensors/HG008T_Illumina/tensors")
RIS = Path("/scratch/jshen/data/pansoma_net_v2_runs/ris_7fc7396_20261012")
LO, HI = 3129540, 3141988
WINDOW = ("chr1", 242_500_000, 243_000_000)
SKIP = {"shard_file", "index_within_shard", "shard_index", "source_index_within_shard", "grch38_path", "source_task",
        "source_shard", "source_index"}  # positions in the shard layout only


def summaries(d, chrom="chr1"):
    out = {}
    for path in sorted(Path(d).glob(f"{chrom}_variant_summary.ndjson")):
        with open(path) as stream:
            for line in stream:
                r = json.loads(line)
                out[r["candidate_id"]] = r
    return out


def tensor(d, r, cache):
    shard = Path(d) / r["shard_file"] if "/" not in r["shard_file"] else Path(r["shard_file"])
    if shard not in cache:
        cache[shard] = np.load(shard, mmap_mode="r")
    return np.asarray(cache[shard][r["index_within_shard"]])


def vcf(path):
    out = {}
    with gzip.open(path, "rt") as stream:
        for line in stream:
            if line.startswith("#"):
                continue
            f = line.split("\t")
            if f[0] == WINDOW[0] and WINDOW[1] <= int(f[1]) < WINDOW[2]:
                info = dict(x.split("=", 1) for x in f[7].split(";") if "=" in x)
                out[(f[0], int(f[1]), f[3], f[4])] = (f[6], float(info["P_SOMATIC"]))
    return out


def main(out):
    out = Path(out)
    report = {}
    for kind in ("SNV", "INDEL"):
        ex_dir, ge_dir = out / "HG008T_Illumina" / "tensors" / kind, GENOME / kind
        ex, ge = summaries(ex_dir), summaries(ge_dir)
        inside = {c for c, r in ex.items() if LO <= r["node_id"] <= HI}
        ge_inside = {c for c, r in ge.items() if LO <= r["node_id"] <= HI}
        same_bytes = same_summary = 0
        cache_e, cache_g, diffs = {}, {}, []
        for c in sorted(inside & ge_inside):
            a, b = tensor(ex_dir, ex[c], cache_e), tensor(ge_dir, ge[c], cache_g)
            same_bytes += a.shape == b.shape and bool(np.array_equal(a, b))
            ra = {k: v for k, v in ex[c].items() if k not in SKIP and k != "grch38"}
            rb = {k: v for k, v in ge[c].items() if k not in SKIP and k != "grch38"}
            if ra == rb:
                same_summary += 1
            elif len(diffs) < 3:
                diffs.append({k: (ra.get(k), rb.get(k)) for k in set(ra) | set(rb) if ra.get(k) != rb.get(k)})
        e_calls = vcf(out / "vcf" / f"HG008T_Illumina.tensors.{kind}.linear.vcf.gz")
        g_calls = vcf(RIS / f"HG008_Illumina_{kind}_d9_normbq" / "vcf_chr1" / f"HG008T_Illumina.tensors.{kind}.linear.vcf.gz")
        both = set(e_calls) & set(g_calls)
        e_pass = {k for k, v in e_calls.items() if v[0] == "PASS"}
        g_pass = {k for k, v in g_calls.items() if v[0] == "PASS"}
        dp = [abs(e_calls[k][1] - g_calls[k][1]) for k in both]
        report[kind] = dict(
            example_tensors=len(ex), example_tensors_in_node_range=len(inside), genome_tensors_in_node_range=len(ge_inside),
            in_both=len(inside & ge_inside), only_example=sorted(inside - ge_inside)[:10],
            only_genome=sorted(ge_inside - inside)[:10], identical_tensor_bytes=same_bytes,
            identical_summaries=same_summary, summary_diffs=diffs,
            window_vcf_records=dict(example=len(e_calls), genome=len(g_calls), both=len(both)),
            window_pass=dict(example=len(e_pass), genome=len(g_pass), both=len(e_pass & g_pass),
                             only_example=sorted(e_pass - g_pass), only_genome=sorted(g_pass - e_pass)),
            p_somatic_abs_diff=dict(max=max(dp) if dp else None, median=float(np.median(dp)) if dp else None))
    print(json.dumps(report, indent=1, default=str))


if __name__ == "__main__":
    main(sys.argv[1])
