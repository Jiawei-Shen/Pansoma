#!/usr/bin/env python3
"""Diagnostic for the ver4 INDEL table (not a letter item): why the graph-represented germline INDEL alleles keep a
candidate on HPRC v2.1 far less often on HCC1395T long reads than on HG008T / COLO829T long reads. The class
"matched-normal germline, >= 46 of 463 haplotypes" of the GRCh38-only FP before the filter, split by the repeat class
of the allele in the E7 carrier table (rep_class) and by the genotype of the germline record (GT of the first allele
match: hom = both alleles ALT). Uses the functions of fp_decomp.py. Writes diag_indel_longread.md."""
import collections, gzip, re, sys
from pathlib import Path
HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
import fp_decomp as F
L = F.L

RC = ("HP>=7", "HP2-6", "STR", "other", "SV", "none")  # rep_class values of the E7 table ("" = none)
cells = [c for c in F.ver4_cells() if c["kind"] == "INDEL" and c["platform"] in ("PacBio", "ONT")]
need = set()
for c in cells:
    c["lin_set"], c["d46_set"] = F.set_of(c["lin_dir"], "INDEL"), F.set_of(c["d46_dir"], "INDEL")
    m = L.manifest(c["lin_set"], "INDEL")
    seq = L.chrom_seq(m["fasta"])
    c["g"] = L.germline(m["germline_vcf"], m["fasta"])
    c["raw"], _ = F.fp_keys(c["lin_dir"], "raw_calls", seq)
    c["tk"], _ = F.tensor_keys(c["d46_set"], "INDEL", seq, keep=c["raw"])
    # zygosity of the germline alleles (normalized key -> hom / het)
    if "zyg" not in c["g"].__dict__:
        z = {}
        for r in L.read_records(m["germline_vcf"]):
            gt = [g for g in re.split(r"[|/]", r.fmt.get("GT", ".")) if g.isdigit()]
            for i in {int(g) for g in gt if int(g) > 0}:
                if i - 1 < len(r.alts):
                    for k in L.keys_of(seq, r.pos, r.ref, r.alts[i - 1]):
                        z.setdefault(k, "hom" if len(gt) == 2 and gt[0] == gt[1] else "het")
        c["g"].zyg = z
    need |= c["raw"]
car, rep = {}, {}
with gzip.open(F.V21, "rt") as fh:
    fh.readline()
    for line in fh:
        f = line.rstrip("\n").split("\t")
        parts = f[0].split(":")
        if len(parts) != 4:
            continue
        key = (int(parts[1]), parts[2], parts[3])
        if key in need:
            car[key], rep[key] = int(f[6]), f[2]
md = ["# Diagnostic: graph-represented germline INDEL alleles (GRCh38-only FP before the filter) that keep a candidate on HPRC v2.1, long reads", "",
      "Cells: n; of them with a v2.1 candidate (%).", "",
      "| Dataset | all | " + " | ".join(f"rep_class {x}" for x in RC) + " | germline GT hom | germline GT het |", "|---|---|" + "---|" * (len(RC) + 2)]
for c in cells:
    cnt = collections.Counter()
    for k in c["raw"]:
        if k in c["g"].keys and car.get(k, 0) >= 46:
            rc = rep[k] or "none"
            assert rc in RC, rc
            for cls in ("all", rc, c["g"].zyg.get(k, "?")):
                cnt[cls, "n"] += 1
                cnt[cls, "t"] += k in c["tk"]
    f = lambda x: f"{cnt[x, 'n']:,}; {cnt[x, 't']:,} ({100 * cnt[x, 't'] / cnt[x, 'n']:.1f}%)" if cnt[x, "n"] else "0"
    md.append(f"| {c['tumor']} {F.pname(c['tumor'], c['platform'])} | " + " | ".join(f(x) for x in ("all",) + RC + ("hom", "het")) + " |")
    print(sorted({k[0] for k in cnt}), file=sys.stderr)
(HERE / "diag_indel_longread.md").write_text("\n".join(md) + "\n")
print("\n".join(md))
