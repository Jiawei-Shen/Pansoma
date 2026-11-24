"""R3.2 (2026-11-04, user): the GRCh38-only graph's false positives BEFORE population filtering, split by whether the
allele is the normal's germline allele and whether the HPRC v2.1 graph carries it, and what happens to the same
allele on HPRC v2.1 with the same reads.

Call sets: the seed-0 7fc7396 runs at the main threshold (validation recall 0.9) of E4b / Table R3.2:
RUNS/<S>_<P>_<K>_{linear,d46}_7fc7396/vcf_chr1/eval_bed/vcfeval/{raw_calls,pon_calls}/fp.vcf.gz (chr1, benchmark BED).
Alleles normalized (calls_lib.normalize; MNPs also as single bases). Germline: the set's germline truth VCF (normal's
dipcall), calls_lib.Germline keys. HPRC v2.1: carriers among the 463 haplotypes (E7 v2.1 hprc_alleles chr1); >= 46 =
the d46 filter threshold (a carrier-frequency proxy for the graph's content). On v2.1: the allele has a chr1 tensor
(labels' grch38 projection), and is a PASS false positive before / after the PoN of the d46 run.
Writes r32_prepon.tsv and r32_prepon.md next to this script."""
import gzip, json, sys, collections
from pathlib import Path
sys.path.insert(0, "/scratch/jshen/Github/Pansoma/analysis/reviewer_response_20261011/experiments/call_level_common")
import calls_lib as L

HERE = Path(__file__).resolve().parent
V21 = Path("/scratch/jshen/data/reviewer_experiments_20261013/E7_graph_composition/v2.1/hprc_alleles/chr1.tsv.gz")
SETS = [("HG008", "HG008T", p) for p in ("Illumina", "PacBio", "ONT")] + [("COLO829", "COLO829T", p) for p in ("Illumina", "fiberseq", "ONT")]


def fp_keys(run, which, seq):
    f = L.RUNS / run / "vcf_chr1/eval_bed/vcfeval" / which / "fp.vcf.gz"
    keys = set()
    for r in L.read_records(f, with_format=False):
        for a in r.alts:
            n = L.normalize(seq, r.pos, r.ref, a)
            if n:
                keys.add(n)
    return keys


def tensor_keys(set_name, kind, seq):
    keys = set()
    for line in open(L.TENSORS / set_name / "tensors" / kind / "chr1_labels.ndjson"):
        g = json.loads(line).get("grch38")
        if not g or g["chrom"] != "chr1":
            continue
        p0, ref, alt = g["pos0"], g["ref"], g["alt"]
        if len(ref) == len(alt):
            n = L.normalize(seq, p0 + 1, ref, alt)
        else:  # INDEL: pad with the base before pos0
            b = seq[p0 - 1]
            n = L.normalize(seq, p0, b + ref, b + alt)
        if n:
            keys.add(n)
    return keys


def main():
    rows, need = [], set()
    jobs = []
    for s, t, p in SETS:
        for kind in ("SNV", "INDEL"):
            lin, d46 = f"{s}_{p}_{kind}_linear_7fc7396", f"{s}_{p}_{kind}_d46_7fc7396"
            m = L.manifest(f"{t}_{p}_linear", kind)
            seq = L.chrom_seq(m["fasta"])
            g = L.germline(m["germline_vcf"], m["fasta"])
            raw, pon = fp_keys(lin, "raw_calls", seq), fp_keys(lin, "pon_calls", seq)
            draw, dpon = fp_keys(d46, "raw_calls", seq), fp_keys(d46, "pon_calls", seq)
            tk = tensor_keys(f"{t}_{p}_d46", kind, seq)
            L.log(t, p, kind, "linear raw FP", len(raw), "d46 raw FP", len(draw), "d46 tensors", len(tk))
            jobs.append((t, p, kind, raw, pon, draw, dpon, tk, g))
            need |= raw
    car = {}
    with gzip.open(V21, "rt") as fh:
        fh.readline()
        for line in fh:
            k, rest = line.split("\t", 1)
            parts = k.split(":")
            if len(parts) != 4:  # symbolic / breakend ALTs
                continue
            key = (int(parts[1]), parts[2], parts[3])
            if key in need:
                car[key] = int(rest.split("\t")[5])
    L.log("v2.1 carrier keys found", len(car), "of", len(need))
    for t, p, kind, raw, pon, draw, dpon, tk, g in jobs:
        c = collections.Counter()
        for k in raw:
            germ = "germline" if k in g.keys else "other"
            n = car.get(k, 0)
            grp = "ge46" if n >= 46 else "1to45" if n > 0 else "0"
            for cls in (f"{germ}|{grp}", f"{germ}|all", "all|all"):
                c[cls, "n"] += 1
                c[cls, "pon_removes"] += k not in pon
                c[cls, "d46_tensor"] += k in tk
                c[cls, "d46_raw_fp"] += k in draw
                c[cls, "d46_pon_fp"] += k in dpon
        for cls in sorted({x for x, _ in c}):
            rows.append(dict(tumor=t, platform=p, kind=kind, cls=cls, **{f: c[cls, f] for f in ("n", "pon_removes", "d46_tensor", "d46_raw_fp", "d46_pon_fp")},
                             linear_raw_fp=len(raw), linear_pon_fp=len(pon), d46_raw_fp_all=len(draw), d46_pon_fp_all=len(dpon)))
    with open(HERE / "r32_prepon.tsv", "w") as f:
        f.write("\t".join(rows[0]) + "\n")
        for r in rows:
            f.write("\t".join(str(v) for v in r.values()) + "\n")
    # compact markdown
    md = ["| dataset | kind | GRCh38-only FP before PoN (after) | HPRC v2.1 FP before PoN (after) | germline allele of the normal | "
          "germline, >= 46 v2.1 haplotypes: n; tensor on v2.1; FP on v2.1 | germline, < 46: n; tensor on v2.1; FP on v2.1 | "
          "not germline, >= 46: n; tensor on v2.1 | not germline, < 46: n | removed by the PoN (GRCh38-only) |", "|---|---|---|---|---|---|---|---|---|---|"]
    by = collections.defaultdict(dict)
    for r in rows:
        by[(r["tumor"], r["platform"], r["kind"])][r["cls"]] = r
    pct = lambda a, b: f"{100 * a / b:.0f}%" if b else "-"
    for (t, p, kind), d in by.items():
        a = d["all|all"]; n = a["n"]
        z = lambda k: d.get(k, dict(n=0, d46_tensor=0, d46_raw_fp=0))
        ge, lt1, lt0 = z("germline|ge46"), z("germline|1to45"), z("germline|0")
        lt = dict(n=lt1["n"] + lt0["n"], d46_tensor=lt1["d46_tensor"] + lt0["d46_tensor"], d46_raw_fp=lt1["d46_raw_fp"] + lt0["d46_raw_fp"])
        oge, o1, o0 = z("other|ge46"), z("other|1to45"), z("other|0")
        md.append(f"| {t} {p} | {kind} | {a['linear_raw_fp']:,} ({a['linear_pon_fp']:,}) | {a['d46_raw_fp_all']:,} ({a['d46_pon_fp_all']:,}) | "
                  f"{z('germline|all')['n']:,} ({pct(z('germline|all')['n'], n)}) | "
                  f"{ge['n']:,} ({pct(ge['n'], n)}); {ge['d46_tensor']:,}; {ge['d46_raw_fp']:,} | "
                  f"{lt['n']:,} ({pct(lt['n'], n)}); {lt['d46_tensor']:,}; {lt['d46_raw_fp']:,} | "
                  f"{oge['n']:,} ({pct(oge['n'], n)}); {oge['d46_tensor']:,} | {o1['n'] + o0['n']:,} ({pct(o1['n'] + o0['n'], n)}) | "
                  f"{a['pon_removes']:,} ({pct(a['pon_removes'], n)}) |")
    (HERE / "r32_prepon.md").write_text("\n".join(md) + "\n")
    print("\n".join(md))


if __name__ == "__main__":
    main()
