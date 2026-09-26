"""Why the INS-span fix changes 37,065 v5 INDEL labels but only 7,142 v6 ones.

Builds the truth sets once (current code), records each allele's placements to also get the old-rule spans
(INS span = [lo, lo+1)), then for every tensor labelled near_truth_allele_mismatch (new rule) checks whether the old
rule also calls it near. The ones that are not are exactly the 0 -> -1 changes; they are characterised by type, length
and distance from the leftmost boundary of the insertion truth that now covers them.
"""
import bisect, collections, glob, json, sys
import numpy as np
sys.path.insert(0, "/scratch/jshen/Github/Pansoma")
from indexed_gam_pipeline_v2.tensor_postprocessing import truth_labels as tl
from indexed_gam_pipeline_v2.tensor_postprocessing.reference_path import ReferencePath

record = []
orig = tl.placements
def spy(genome, pos, ref, alt, kind):
    places = orig(genome, pos, ref, alt, kind)
    record.append((pos, len(ref), kind, places))
    return places
tl.placements = spy

T = "/scratch/jshen/data/HG008_GIAB"
path = ReferencePath(f"/scratch/jshen/data/pansoma_v2_tensors/graph_index/hprc-v1.1-mc-grch38.d9.grch38_path")
locator = tl.Locator(path)
fasta = "/scratch/jshen/data/HapMap/GCA_000001405.15_GRCh38_no_alt_analysis_set.fasta"
sets = {}
for name, vcf, bed in (("somatic", f"{T}/draft_v02_benchmark/HG008-T_somatic_smvar_benchmark_v0.2_tumorvariants.vcf.gz",
                        f"{T}/draft_v02_benchmark/HG008-T_somatic_smvar_benchmark_v0.2_all.bed"),
                       ("germline", f"{T}/dipcall_HG008N_GRCh38/HG008N_GRCh38_dipcall.dip.vcf.gz",
                        f"{T}/dipcall_HG008N_GRCh38/HG008N_GRCh38_dipcall.dip.bed")):
    record.clear()
    ts = tl.TruthSet(name, vcf, bed, fasta, locator)
    assert len(record) == len(ts.alleles)
    old, ins = collections.defaultdict(list), collections.defaultdict(list)
    for a, (pos, rlen, kind, places) in zip(ts.alleles, record):
        lo = min((p[0] for p in places), default=pos)
        hi = max((p[0] + len(p[1]) for p in places), default=pos + rlen)
        old[a["chrom"]].append((lo, lo + 1 if kind == "INS" else max(hi, lo + 1)))
        if kind == "INS":
            ins[a["chrom"]].append((lo, max(hi, lo + 1)))
    for c in old:
        old[c].sort(); ins[c].sort()
    sets[name] = (ts, {c: (np.array([x[0] for x in v]), np.maximum.accumulate(np.array([x[1] for x in v]))) for c, v in old.items()},
                  {c: v for c, v in ins.items()})
    print(name, "loaded", len(ts.alleles), flush=True)

def near(spans, chrom, s, e, d=tl.NEAR_BP):
    if chrom not in spans: return False
    lo, maxhi = spans[chrom]
    k = int(np.searchsorted(lo, e + d, side="left"))
    return bool(k > 0 and maxhi[k - 1] > s - d)

def ins_distance(chrom, s, e):
    """Distance of the tensor start from the leftmost boundary of the nearest covering insertion truth."""
    best = None
    for name in ("somatic", "germline"):
        items = sets[name][2].get(chrom, [])
        k = bisect.bisect_left(items, (e + tl.NEAR_BP + 1,))
        for lo, hi in items[max(0, k - 200):k]:
            if hi > s - tl.NEAR_BP:
                d = s - lo
                best = d if best is None or abs(d) < abs(best) else best
    return best

def bucket(d):
    if d is None: return "none"
    for b in (10, 20, 50, 100):
        if d <= b: return f"<={b}"
    return ">100"

ROOT = "/scratch/jshen/data/pansoma_v2_tensors/Liss_lab_PacBio_Revio_20240125"
for version, directory in (("v5", f"{ROOT}/v5_tensors"), ("v6", f"{ROOT}/v6_tensors")):
    for kind in ("SNV", "INDEL"):
        changed, near_total, types, dist = 0, 0, collections.Counter(), collections.Counter()
        for f in glob.glob(f"{directory}/{kind}/chr*_labels.ndjson"):
            for line in open(f):
                if "near_truth_allele_mismatch" not in line: continue
                l = json.loads(line); g = l["grch38"]; near_total += 1
                ev = l["candidate_id"].split(":")[2]
                s, e = tl.linear_interval(g, ev)
                if near(sets["somatic"][1], g["chrom"], s, e) or near(sets["germline"][1], g["chrom"], s, e):
                    continue
                changed += 1
                ref, alt = l["candidate_id"].split(":", 3)[3].split("@")[0].split(">")
                types[ev if ev == "SNP" else f"{ev}{'1' if abs(len(alt) - len(ref)) == 1 else '>1'}"] += 1
                dist[bucket(ins_distance(g["chrom"], s, e))] += 1
        print(version, kind, "near_truth", near_total, "changed_by_fix", changed, dict(types),
              "distance_from_INS_leftmost", dict(dist), flush=True)
