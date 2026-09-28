"""Replay the haplotype overlap of the sampled HG008 Illumina partial-somatic SNV tensors, read by read (pipeline_code, v4)."""
import csv, json, sys
sys.path.insert(0, "/scratch/jshen/data/pansoma_v2_tensors/pipeline_code")  # v4 snapshot (same rules as truth-labels-v6)
from indexed_gam_pipeline_v4.tensor_postprocessing import truth_labels as tl
from indexed_gam_pipeline_v4.tensor_postprocessing.reference_path import rc
csv.field_size_limit(sys.maxsize)
S = "/scratch/jshen/data/pansoma_v2_tensors/Liss_lab_BCM_Illumina-WGS_20240313"
FA = "/scratch/jshen/data/HapMap/GCA_000001405.15_GRCh38_no_alt_analysis_set.fasta"
truth = {}
for r in csv.DictReader(open(f"{S}/truth/somatic.graph.tsv"), delimiter="\t"):
    truth.setdefault((r["chrom"], r["vcf_pos"]), []).append(r)
ev = tl.Evidence(f"{S}/v3_tensors/SNV", FA)
for e in json.load(open("partial_snv_examples.json")):
    chrom, rest = e["truth"].split(":"); vpos = rest.split()[0]
    t = truth[(chrom, vpos)][0]; pos0 = int(t["pos0"])
    rec = next(json.loads(l) for l in open(f"{S}/v3_tensors/SNV/{e['chrom']}_variant_summary.ndjson") if f'"candidate_id": "{e["cid"]}"' in l)
    alt_rows, ref_rows = ev.rows(rec)
    if rec["grch38"] and rec["grch38"]["node_reverse"]:
        alt_rows, ref_rows = [rc(q) for q in alt_rows], [rc(q) for q in ref_rows]
    start = max(0, pos0 - tl.HAPLOTYPE_WINDOW)
    reference = ev.window(chrom, start, pos0 + len(t["ref"]) + tl.HAPLOTYPE_WINDOW)
    o = pos0 - start; hap = reference[:o] + t["alt"] + reference[o + len(t["ref"]):]
    event = 1 if t["kind"] == "SNP" else max(len(t["ref"]), len(t["alt"]))
    fits = lambda q: (tl.fit_distance(q, reference), tl.fit_distance(q, hap))  # noqa: E731
    base = sorted(min(fits(q)) for q in ref_rows); b = base[len(base) // 2] if base else 0
    d = [fits(q) for q in alt_rows]
    print(f"\n{e['cid']} ({e['chrom']}:{e['pos']}, AF {e['af']:.3f}, ALT {e['alt']}/{e['cov']}, strands fwd {e['fwd']} rev {e['rev']})"
          f"  truth {chrom}:{vpos} {t['vcf_ref'][:15]}>{t['vcf_alt'][:15]} ({t['kind']}, event {event}), overlap {e['ov']}, truth has own tensor: {e['own']}")
    snv_off = (e["pos"] - 1 - start) if e["pos"] else None
    print(f"   GRCh38 around truth: {reference[max(0,o-20):o]}|{reference[o:o+len(t['ref'])]}|{reference[o+len(t['ref']):o+len(t['ref'])+20]}"
          + (f"   SNV at {e['pos'] - int(vpos):+d} bp from truth" if e["pos"] else ""))
    print(f"   REF-read median error b={b}; A1 reads (d_ref, d_truth): {d[:10]}")
