"""For every partial-somatic tensor (truth-labels-v6), recompute the raw read evidence behind the haplotype overlap:
median over A1 reads of d_ref - d_truth (no background subtraction). 'no evidence' = median <= 0 (the reads are not
closer to the truth haplotype than to GRCh38); such tensors get overlap 0.5 whenever the REF-read error b >= d_ref."""
import csv, glob, json, sys, collections, statistics
sys.path.insert(0, "/scratch/jshen/data/pansoma_v2_tensors/pipeline_code")
from indexed_gam_pipeline_v4.tensor_postprocessing import truth_labels as tl
from indexed_gam_pipeline_v4.tensor_postprocessing.reference_path import rc
csv.field_size_limit(sys.maxsize)
R = "/scratch/jshen/data/pansoma_v2_tensors"; FA = "/scratch/jshen/data/HapMap/GCA_000001405.15_GRCh38_no_alt_analysis_set.fasta"
SETS = {"HG008 Illumina": "Liss_lab_BCM_Illumina-WGS_20240313", "HG008 PacBio": "Liss_lab_PacBio_Revio_20240125",
        "HG008 ONT": "Liss_lab_Northeastern-ONT-UL-20241216", "COLO829T Illumina": "COLO829T_Illumina",
        "COLO829T ONT": "COLO829T_ONT", "COLO829T fiberseq": "COLO829T_fiberseq"}
for name, d in SETS.items():
    S = f"{R}/{d}"; truth = {}
    for r in csv.DictReader(open(f"{S}/truth/somatic.graph.tsv"), delimiter="\t"):
        truth[r["truth_id"]] = r
    out = collections.Counter()
    for kind in ("SNV", "INDEL"):
        ev = tl.Evidence(f"{S}/v3_tensors/{kind}", FA)
        for f in sorted(glob.glob(f"{S}/v3_tensors/{kind}/chr*_labels.ndjson")):
            for a, b in zip(open(f), open(f.replace("_labels.ndjson", "_variant_summary.ndjson"))):
                if "partial_somatic_truth" not in a: continue
                l = json.loads(a); rec = json.loads(b); t = truth[str(l["partial_truth"]["truth_id"])]
                ov = l["overlap"]; ob = "0.45-0.6" if ov <= 0.6 else ">0.6"
                pos0 = int(t["pos0"]); lin = l["grch38"]
                # allele overlap alone decides when it is > 0.45 (INDEL vs INDEL of one kind): not the haplotype path
                if lin and t["kind"] in ("INS", "DEL") and rec["event_type"] == t["kind"]:
                    a_ov = tl.allele_overlap(rec["event_type"], lin["pos0"], lin["ref"], lin["alt"],
                                             dict(kind=t["kind"], ref=t["ref"], lo=int(t["pos0"]), hi=int(t["pos0"]) + len(t["ref"]), left_alt=t["alt"]))
                else:
                    a_ov = 0.0
                if a_ov > tl.MIN_OVERLAP:
                    out[(kind, ob, "allele overlap")] += 1; continue
                alt_rows, ref_rows = ev.rows(rec)
                if lin and lin["node_reverse"]:
                    alt_rows, ref_rows = [rc(q) for q in alt_rows], [rc(q) for q in ref_rows]
                start = max(0, pos0 - tl.HAPLOTYPE_WINDOW)
                ref = ev.window(t["chrom"], start, pos0 + len(t["ref"]) + tl.HAPLOTYPE_WINDOW)
                o = pos0 - start; hap = ref[:o] + t["alt"] + ref[o + len(t["ref"]):]
                def fits(q):
                    pairs = [(tl.fit_distance(x, ref), tl.fit_distance(x, hap)) for x in ([q] if lin else [q, rc(q)])]
                    return min(pairs, key=min)
                diffs = [dr - dt for dr, dt in (fits(q) for q in alt_rows)]
                if not diffs:
                    out[(kind, ob, "no A1 reads?")] += 1; continue
                m = statistics.median(diffs)
                event = 1 if t["kind"] == "SNP" else max(len(t["ref"]), len(t["alt"]))
                cls = "no evidence (median d_ref-d_truth <= 0)" if m <= 0 else ("reads carry >= half the truth" if m >= event / 2 else "weak evidence")
                out[(kind, ob, cls)] += 1
    print(f"\n=== {name}")
    for k in sorted(out): print(f"   {k[0]:5s} overlap {k[1]:8s} {k[2]:42s} {out[k]:6d}")
