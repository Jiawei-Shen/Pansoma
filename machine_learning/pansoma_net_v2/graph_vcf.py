"""pansoma_net_v2 predictions -> a VCF in graph coordinates, one record per tensor on the node its event starts on.

    cd machine_learning
    python -m pansoma_net_v2.graph_vcf --predictions <predict output>/<name>.SNV.predictions.ndjson.gz \
        --tensors <set> --kind SNV --graph-index <graph>.graph_index.sqlite --output <out>/<name>.SNV.graph.vcf.gz

Inputs: the predictions of one kind (predict.py) and the set's summaries, which give each tensor's event and read
support. Record k of chromosome C in the predictions is line k of <set>/<KIND>/C_variant_summary.ndjson; the
candidate IDs must agree. Labels and truth are not read.

A record is the site's representative allele A1, the allele the model scores. The other alleles of the site
(INFO NALLELES > 1) get no probability of their own and no record.
- CHROM is the node ID, POS a 1-based offset on the node's forward strand, ID the candidate ID.
- SNP: POS = start + 1, and REF/ALT are the two bases.
- INS and DEL: VCF padding with the node base before the event (POS = start). At node offset 0 the padding is
  the base after the event, with POS = 1 (the VCF rule for an event at position 1): the node's first base for
  an INS, and the base after the deleted bases for a DEL. A deletion that ends at the end of its last node has
  no known base after it; it is padded with N and flagged NOPAD. A deletion over several nodes (INFO PATH)
  has a REF that runs on over the further nodes.
- INFO keeps the event exactly (START, KIND, EVREF, EVALT, PATH), and that is what linear_vcf reads. It also
  keeps the class probabilities and the site. FORMAT has the read support of A1: GT is always 0/1 (the
  model does not genotype), plus DP, AD and AF. On a node over the read cap (max_node_reads) these count the
  sampled records only, and INFO NODE_RECORDS gives the node's records before the sample.
- QUAL = -10 log10(1 - p_somatic), which is 50 at p = 1 (the probabilities have 5 decimals). rtg vcfeval
  rounds scores to 3 decimals, and on the phred scale the probabilities near 1 stay apart.
- FILTER:
  - PASS: called somatic and AF >= the AF floor.
  - LowQual: not called somatic. Without --threshold, the call is predict's own (`pred`, made with the
    checkpoint's threshold before the probabilities were rounded), so PASS is exactly predict's somatic call
    set; the threshold is read from the .metrics.json beside the predictions for the header. With
    --threshold, the call is p_somatic >= threshold.
  - LowAF: AF < --min-af. The default is the label floor of the kind that predict's test used (the labels
    entry of its .metrics.json: snv_min_af / indel_min_af; e.g. SNV 0.07 on the Illumina sets, none on the
    long-read sets), else the set's labels.manifest.json, else none: the labels ignore those tensors
    (below_snv_min_af) and the test leaves them out. A relabel of the set after predict does not change it.
- Tensors with p_somatic < --min-score (default 0.01) that are not called are not written. The LowQual records
  are kept so that vcfeval can draw the precision/recall curve below the threshold.
"""
import argparse
import gzip
import json
import math
import re
from pathlib import Path

from .vcfio import format_info, import_pipeline, stats_path, write_vcf

KIND_EVENTS = {"SNV": {"SNP"}, "INDEL": {"INS", "DEL"}}
FIRST_ID = re.compile(r'\{"candidate_id": ?"([^"\\]*)"')
PREDICTION_FIELDS = ("candidate_id", "p_somatic", "p_germline", "p_non", "pred")
SUMMARY_FIELDS = ("candidate_id", "node_id", "start", "ref", "alt", "event_type", "path", "chrom", "site_id",
                  "allele_count", "coverage", "ref_count", "alt_count", "af", "downsampled_from")  # the rest is kilobytes
LABEL_FLOOR = {"SNV": "snv_min_af", "INDEL": "indel_min_af"}
MAX_QUAL = 50.0

HEADER = """##FILTER=<ID=PASS,Description="Called somatic (##pansoma_call) and AF at or above the AF floor">
##FILTER=<ID=LowQual,Description="Not called somatic: p_somatic below the somatic threshold">
##FILTER=<ID=LowAF,Description="A1 allele fraction below the AF floor of the training labels">
##INFO=<ID=BLOCK,Number=1,Type=String,Description="Chromosome block of the node (tensor_postprocessing chr_index)">
##INFO=<ID=START,Number=1,Type=Integer,Description="0-based start of the event on the node's forward strand">
##INFO=<ID=KIND,Number=1,Type=String,Description="Event type: SNP, INS or DEL">
##INFO=<ID=EVREF,Number=1,Type=String,Description="Replaced or deleted bases, node-forward, without the padding base (absent for an INS)">
##INFO=<ID=EVALT,Number=1,Type=String,Description="Substituted or inserted bases, node-forward, without the padding base (absent for a DEL)">
##INFO=<ID=PATH,Number=.,Type=String,Description="Further nodes of a deletion over several nodes: node:start:end, 0-based half-open, forward strand">
##INFO=<ID=NOPAD,Number=0,Type=Flag,Description="Deletion from node offset 0 to the end of its last node: the base after it is unknown, padded with N">
##INFO=<ID=P_SOMATIC,Number=1,Type=Float,Description="Model probability of somatic">
##INFO=<ID=P_GERMLINE,Number=1,Type=Float,Description="Model probability of germline">
##INFO=<ID=P_NON,Number=1,Type=Float,Description="Model probability of non-variant">
##INFO=<ID=SITE,Number=1,Type=String,Description="Tensor site: node:start:SNV|INDEL">
##INFO=<ID=NALLELES,Number=1,Type=Integer,Description="Passing alleles at the site; the model scores A1, this record, only">
##INFO=<ID=NODE_RECORDS,Number=1,Type=Integer,Description="Records on the node before the max_node_reads sample; FORMAT DP and AD count reads within the sample">
##FORMAT=<ID=GT,Number=1,Type=String,Description="Always 0/1: the model does not genotype">
##FORMAT=<ID=DP,Number=1,Type=Integer,Description="Reads at the site: A1 + REF + other alleles (within the sample on a node with NODE_RECORDS)">
##FORMAT=<ID=AD,Number=R,Type=Integer,Description="Reads supporting REF and A1">
##FORMAT=<ID=AF,Number=A,Type=Float,Description="A1 allele fraction: A1 reads / DP">"""


def phred(p):
    return min(MAX_QUAL, -10.0 * math.log10(max(1.0 - p, 1e-12)))


def read_predictions(path):
    """Yield (chrom, [records]) per chromosome block of a predictions file (predict.py order)."""
    seen, chrom, block = set(), None, []
    with gzip.open(path, "rt") as f:
        for line in f:
            p = json.loads(line)
            if p["chrom"] != chrom:
                if block:
                    yield chrom, block
                chrom, block = p["chrom"], []
                if chrom in seen:
                    raise ValueError(f"{path}: the records of {chrom} are not contiguous")
                seen.add(chrom)
            block.append(p)
    if block:
        yield chrom, block


def predict_report(predictions):
    """The .metrics.json predict wrote beside the predictions, or {}."""
    name = Path(predictions).name
    if not name.endswith(".predictions.ndjson.gz"):
        return {}
    report = Path(predictions).with_name(name[:-len(".predictions.ndjson.gz")] + ".metrics.json")
    return json.loads(report.read_text()) if report.exists() else {}


def default_threshold(predictions):
    """The threshold predict applied."""
    return predict_report(predictions).get("threshold")


def default_min_af(predictions, tensors, kind):
    """(AF floor, source): the labels of predict's test, else the set's labels, else none."""
    key = LABEL_FLOOR[kind]
    labels = predict_report(predictions).get("labels") or {}
    if key in labels:
        return labels[key], f".metrics.json labels {key}"
    manifest = Path(tensors) / kind / "labels.manifest.json"
    if manifest.exists():
        return json.loads(manifest.read_text()).get(key), f"labels.manifest.json {key}"
    return None, "none (no labels)"


def select(predictions, summaries_dir, kind, called, min_score, threshold, stats):
    """The kept (prediction, summary) pairs, and the predicted chromosomes in file order."""
    kept, chroms = [], []
    for chrom, block in read_predictions(predictions):
        chroms.append(chrom)
        summary_path = Path(summaries_dir) / f"{chrom}_variant_summary.ndjson"
        with open(summary_path) as f:
            for k, p in enumerate(block):
                line = f.readline()
                if not line:
                    raise ValueError(f"{summary_path} has {k} records, the predictions {len(block)}")
                if p["p_somatic"] < min_score and not called(p):
                    m = FIRST_ID.match(line)  # the first key; the full parse is only for kept records
                    if (m.group(1) if m else json.loads(line)["candidate_id"]) != p["candidate_id"]:
                        raise ValueError(f"{summary_path} line {k + 1} is not {p['candidate_id']}")
                    continue
                s = json.loads(line)
                if s["candidate_id"] != p["candidate_id"]:
                    raise ValueError(f"{summary_path} line {k + 1} is {s['candidate_id']}, the prediction {p['candidate_id']}")
                if s["event_type"] not in KIND_EVENTS[kind]:
                    raise ValueError(f"{p['candidate_id']}: event {s['event_type']} is not {kind}")
                kept.append(({k: p[k] for k in PREDICTION_FIELDS}, {k: s.get(k) for k in SUMMARY_FIELDS}))
            if f.readline():
                raise ValueError(f"{summary_path} has more records than the {len(block)} predictions of {chrom}")
        stats["predictions"] += len(block)
        stats["called_by_predict"] += sum(p["pred"] == "somatic" for p in block)
        stats["called"] += sum(called(p) for p in block)
        stats["pred_differs_from_threshold"] += sum((p["pred"] == "somatic") != (p["p_somatic"] >= threshold)
                                                    for p in block)
    return kept, chroms


def event_record(s, sequences):
    """(POS, REF, ALT, INFO event fields) of a summary record in node coordinates; checks the bases against the
    graph."""
    node, start, ref, alt, kind = int(s["node_id"]), int(s["start"]), s["ref"], s["alt"], s["event_type"]
    path = [tuple(int(v) for v in step) for step in s.get("path") or []]
    seq = sequences[node]
    first = len(ref) - sum(e - b for _, b, e in path)
    graph_ref = seq[start:start + first] + "".join(sequences[n][b:e] for n, b, e in path)
    if (kind == "INS" and (ref or not alt or start > len(seq))) or (kind != "INS" and (graph_ref != ref or first < 0)):
        raise ValueError(f"{s['candidate_id']}: event does not fit the graph (node bases {graph_ref!r})")
    if kind == "SNP" and (len(ref) != 1 or len(alt) != 1):
        raise ValueError(f"{s['candidate_id']}: an SNP of {len(ref)} > {len(alt)} bases")
    info = dict(START=start, KIND=kind, EVREF=ref or None, EVALT=alt or None,
                PATH=",".join(f"{n}:{b}:{e}" for n, b, e in path) or None, NOPAD=None)
    if kind == "SNP":
        return start + 1, ref, alt, info
    if start > 0:
        pad = seq[start - 1]
        return start, pad + ref, pad + alt, info
    if kind == "INS":
        return 1, seq[0], alt + seq[0], info
    last, end = (path[-1][0], path[-1][2]) if path else (node, start + first)
    if end < len(sequences[last]):
        pad = sequences[last][end]
    else:
        pad, info["NOPAD"] = "N", True
    return 1, ref + pad, alt + pad, info


def build(args):
    import_pipeline()
    from indexed_gam_pipeline_v4.graph_index import GraphIndex

    threshold = args.threshold if args.threshold is not None else default_threshold(args.predictions)
    if threshold is None:
        raise SystemExit("no --threshold, and no .metrics.json with the threshold beside the predictions")
    if args.threshold is None:
        call_rule = "predict's call (pred)"
        called = lambda p: p["pred"] == "somatic"  # noqa: E731
    else:
        call_rule = f"p_somatic >= {threshold}"
        called = lambda p: p["p_somatic"] >= threshold  # noqa: E731
    if args.min_af is not None:
        min_af, min_af_source = args.min_af, "--min-af"
    else:
        min_af, min_af_source = default_min_af(args.predictions, args.tensors, args.kind)
    sample = args.sample or Path(args.tensors).resolve().parent.name
    stats = dict(predictions=0, called_by_predict=0, called=0, pred_differs_from_threshold=0)
    kept, chroms = select(args.predictions, Path(args.tensors) / args.kind, args.kind, called, args.min_score,
                          threshold, stats)

    nodes = {int(s["node_id"]) for _, s in kept} | {int(n) for _, s in kept for n, _, _ in s.get("path") or []}
    with GraphIndex(args.graph_index) as graph:
        sequences = {n: v["sequence"] for n, v in graph.get_nodes(nodes).items()}
    records, filters = [], {"PASS": 0, "LowQual": 0, "LowAF": 0}
    for p, s in kept:
        pos, ref, alt, event = event_record(s, sequences)
        flags = [name for name, bad in (("LowQual", not called(p)),
                                        ("LowAF", min_af is not None and s["af"] < min_af)) if bad]
        for name in flags or ["PASS"]:
            filters[name] += 1
        info = dict(BLOCK=s["chrom"], **event, P_SOMATIC=f"{p['p_somatic']:.5f}", P_GERMLINE=f"{p['p_germline']:.5f}",
                    P_NON=f"{p['p_non']:.5f}", SITE=s["site_id"], NALLELES=s["allele_count"],
                    NODE_RECORDS=s["downsampled_from"])
        sample_field = f"0/1:{s['coverage']}:{s['ref_count']},{s['alt_count']}:{s['af']:.6g}"
        records.append((int(s["node_id"]), pos, [str(s["node_id"]), str(pos), s["candidate_id"], ref, alt,
                                                 f"{phred(p['p_somatic']):.3f}", ";".join(flags) or "PASS",
                                                 format_info(info), "GT:DP:AD:AF", sample_field]))
    records.sort(key=lambda r: (r[0], r[1], r[2][3], r[2][4]))
    used = sorted({r[0] for r in records})
    header = ["##fileformat=VCFv4.2", "##source=pansoma_net_v2.graph_vcf",
              "##pansoma_coordinates=graph: CHROM node ID, POS 1-based on the node's forward strand",
              f"##pansoma_kind={args.kind}", f"##pansoma_chromosomes={','.join(chroms)}",
              f"##pansoma_predictions={Path(args.predictions).resolve()}",
              f"##pansoma_tensors={Path(args.tensors).resolve()}",
              f"##pansoma_graph_index={Path(args.graph_index).resolve()}",
              f"##pansoma_threshold={threshold}", f"##pansoma_call={call_rule}",
              f"##pansoma_min_score={args.min_score}", f"##pansoma_min_af={min_af}",
              f"##pansoma_min_af_source={min_af_source}", *HEADER.split("\n"),
              *(f"##contig=<ID={n},length={len(sequences[n])}>" for n in used)]
    columns = "\t".join(["#CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER", "INFO", "FORMAT", sample])
    write_vcf(args.output, header, columns, (r[2] for r in records))
    stats.update(kind=args.kind, chromosomes=chroms, threshold=threshold, call=call_rule, min_score=args.min_score,
                 min_af=min_af, min_af_source=min_af_source,
                 records=len(records), filters=filters, nopad=sum(";NOPAD" in r[2][7] for r in records),
                 multi_allele_sites=sum(s["allele_count"] > 1 for _, s in kept), output=str(Path(args.output).resolve()))
    stats_path(args.output).write_text(json.dumps(stats, indent=2) + "\n")
    print(f"{args.output}: {stats['predictions']:,} predictions, {len(records):,} records "
          f"(PASS {filters['PASS']:,}, LowQual {filters['LowQual']:,}, LowAF {filters['LowAF']:,}); "
          f"threshold {threshold:.5f}, min score {args.min_score}, min AF {min_af}", flush=True)
    return stats


def parse_args(argv=None):
    p = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    p.add_argument("--predictions", required=True, help="<name>.<KIND>.predictions.ndjson.gz of predict")
    p.add_argument("--tensors", required=True, help="the tensor set predicted (directory with SNV/ and INDEL/)")
    p.add_argument("--kind", required=True, choices=sorted(KIND_EVENTS))
    p.add_argument("--graph-index", required=True, help="graph SQLite with the node sequences (indexed_gam_pipeline_v4)")
    p.add_argument("--output", required=True, help="<name>.graph.vcf.gz")
    p.add_argument("--threshold", type=float, help="somatic threshold (default: predict's calls)")
    p.add_argument("--min-score", type=float, default=0.01,
                   help="leave out tensors with a lower p_somatic that are not called")
    p.add_argument("--min-af", type=float,
                   help="AF floor (default: the labels' snv_min_af / indel_min_af of predict's .metrics.json, "
                        "else of the set's labels.manifest.json)")
    p.add_argument("--sample", help="sample column name (default: the name of the set's parent directory)")
    return p.parse_args(argv)


def main(argv=None):
    build(parse_args(argv))


if __name__ == "__main__":
    main()
