"""Are the multi-cancer models' HP>=7 / STR false INDEL calls on HCC1395T real variants the SEQC2 truth lacks, or
platform artifacts? A real variant is in the tumor's DNA and so in every platform's reads; an ONT or PacBio error is not.
For each platform's multicancer4 false calls in HP>=7 / STR (PASS after the PoN, chr1), the share also present in the
other two platforms' call sets: Pansoma multicancer4 PASS (TP + FP after the PoN) and DeepSomatic-TO / ClairS-TO records
of any FILTER (calls.all), exactly (same POS / REF / ALT, anchor removed) or nearby (same kind, within 5 bp).
Standard library only. Writes crossplat.md next to this script.
"""
import importlib.util
from collections import defaultdict
from pathlib import Path

HERE = Path(__file__).resolve().parent
spec = importlib.util.spec_from_file_location("indel_fp", HERE / "indel_fp.py")
F = importlib.util.module_from_spec(spec)
spec.loader.exec_module(F)
spec2 = importlib.util.spec_from_file_location("hcc", HERE / "hcc1395_indel.py")
H = importlib.util.module_from_spec(spec2)
spec2.loader.exec_module(H)
PLATFORMS = ("Illumina", "PacBio", "ONT")


def kind(ref, alt):
    return "DEL" if ref and not alt else "INS" if alt and not ref else "CPX"


def index(keys):
    near = defaultdict(set)
    for c, p, r, a in keys:
        near[(c, kind(r, a))].add(p)
    return set(keys), near


def hit(k, idx, bp=5):
    exact, near = idx
    if k in exact:
        return "exact"
    c, p, r, a = k
    ps = near.get((c, kind(r, a)), ())
    return "near" if any(abs(p - q) <= bp for q in ps) else None


def main():
    genome = F.Fasta(F.FASTA)
    pans, comp = {}, {}
    for p in PLATFORMS:
        d = H.RUNS / f"multicancer4_{p}_INDEL" / f"vcf_chr1_HCC1395T_{p}_d46_normfix" / "eval_bed" / "vcfeval" / "pon_calls"
        fp = F.vcf_records(d / "fp.vcf.gz")
        tp = F.vcf_records(d / "tp.vcf.gz")
        pans[p] = dict(fp=fp, all=index([H.key(r) for r in fp + tp]))
        comp[p] = {c: index(list(H.raw_filters(H.SCORE / c / f"HCC1395T_{p}" / "INDEL" / "calls.all.vcf.gz")))
                   for c in ("deepsomatic", "clairs_to")}
    L = ["# HCC1395T chr1: multicancer4's HP>=7 / STR false INDEL calls on the other platforms\n",
         "(`crossplat.py`. Exact = same POS / REF / ALT; near = same kind (DEL / INS) within 5 bp. Pansoma = the other "
         "platform's multicancer4 PASS calls after the PoN (TP + FP); DS / CT = DeepSomatic-TO / ClairS-TO records of any "
         "FILTER, RefCall excluded.)\n",
         "| platform | HP>=7 / STR false calls | in Pansoma of the other platforms: exact / near | in DS of the other platforms: "
         "exact / near | in CT of the other platforms: exact / near | in any of them, near |", "|---|---:|---|---|---|---|"]
    for p in PLATFORMS:
        others = [q for q in PLATFORMS if q != p]
        rep = [H.key(r) for r in pans[p]["fp"] if F.repeat_class(genome, r[0], r[1], r[2], r[3]) in H.REPEAT]
        cells, anyhit = [], set()
        for src in ("pansoma", "deepsomatic", "clairs_to"):
            ex = nr = 0
            for k in rep:
                hs = [hit(k, pans[q]["all"] if src == "pansoma" else comp[q][src]) for q in others]
                if "exact" in hs:
                    ex += 1
                if any(hs):
                    nr += 1
                    anyhit.add(k)
            cells.append(f"{ex:,} ({ex / max(1, len(rep)):.1%}) / {nr:,} ({nr / max(1, len(rep)):.1%})")
        L.append(f"| {p} | {len(rep):,} | " + " | ".join(cells) + f" | {len(anyhit):,} ({len(anyhit) / max(1, len(rep)):.1%}) |")
        print(p, "done", flush=True)
    (HERE / "crossplat.md").write_text("\n".join(L) + "\n")
    print("\n".join(L))


if __name__ == "__main__":
    main()
