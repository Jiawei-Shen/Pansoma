"""A graph VCF of graph_vcf -> a GRCh38 VCF, plus the records that have no GRCh38 position.

    cd machine_learning
    python -m pansoma_net_v2.linear_vcf --graph-vcf <name>.SNV.graph.vcf.gz \
        --reference-path <graph_index>/<graph>.grch38_path --fasta GRCh38.fasta --output <name>.SNV.linear.vcf.gz

Every record's event (INFO START, KIND, EVREF, EVALT, PATH on node CHROM) is projected with
tensor_postprocessing's ReferencePath.linear, the same projection that merge stored as the summaries' grch38.
- A node with no unique GRCh38 visit, or a deletion whose nodes are not neighbours on GRCh38, has no position.
  Its record goes unchanged to <output>.unplaced.vcf.gz, in graph coordinates. Off-reference events need the
  graph's edges to be placed, and the reference path has none.
- The projected REF must equal the FASTA; any difference stops the run, because it means that the reference
  path and the FASTA disagree.
- INS and DEL are left-aligned on the FASTA (INFO LEFTSHIFT = bases moved) and padded with the FASTA base
  before them (at contig position 1, with the base after them).
- Records that become the same (CHROM, POS, REF, ALT) are merged, for example one insertion placed at the end
  of one node by some reads and at the start of the next node by others. The record that stays is a PASS one if
  any, else the most probable; the others' candidate IDs are added to its ID, and INFO MERGED counts the
  records (the kept one included). Its FORMAT (DP, AD, AF) stays its own: the other records' A1 reads are
  inside its REF and DP there, so AF understates a merged allele. INFO MERGED_ALT sums the A1 reads.
- INFO NODE is the graph node. The graph INFO, FILTER, QUAL and sample fields are kept.
Records are sorted in FASTA contig order, and ##contig lines come from the FASTA index.
"""
import argparse
import json
from pathlib import Path

import pysam

from .vcfio import format_info, import_pipeline, meta, parse_info, read_vcf, stats_path, write_vcf

EXTRA_INFO = """##INFO=<ID=NODE,Number=1,Type=Integer,Description="Graph node of the event (the graph VCF's CHROM)">
##INFO=<ID=LEFTSHIFT,Number=1,Type=Integer,Description="Bases the event moved left when left-aligned on the reference">
##INFO=<ID=MERGED,Number=1,Type=Integer,Description="Graph records of this GRCh38 allele, this one included; ID lists their candidate IDs">
##INFO=<ID=MERGED_ALT,Number=1,Type=Integer,Description="A1 reads summed over the merged records (FORMAT DP/AD/AF are the kept record's)">"""
WINDOW = 256


class Genome:
    """Upper-case FASTA bases with a cached window per contig for the left-alignment walks."""

    def __init__(self, fasta):
        self.fasta = pysam.FastaFile(str(fasta))
        self.lengths = dict(zip(self.fasta.references, self.fasta.lengths))

    def get(self, chrom, start, end):
        return self.fasta.fetch(chrom, max(0, start), min(end, self.lengths[chrom])).upper()

    def left_align(self, chrom, pos0, ref, alt, kind):
        """Leftmost equivalent (pos0, ref, alt) of an INS (alt) or DEL (ref) on the contig; bases moved."""
        moved = 0
        while pos0 > 0:
            window = self.get(chrom, pos0 - WINDOW, pos0)
            for base in reversed(window):
                if kind == "DEL" and base == ref[-1]:
                    ref = base + ref[:-1]
                elif kind == "INS" and base == alt[-1]:
                    alt = alt[-1] + alt[:-1]
                else:
                    return pos0, ref, alt, moved
                pos0, moved = pos0 - 1, moved + 1
        return pos0, ref, alt, moved

    def vcf_allele(self, chrom, pos0, ref, alt, kind):
        """(POS, REF, ALT) with VCF padding."""
        if kind == "SNP":
            return pos0 + 1, ref, alt
        if pos0 > 0:
            pad = self.get(chrom, pos0 - 1, pos0)
            return pos0, pad + ref, pad + alt
        pad = self.get(chrom, len(ref), len(ref) + 1)
        return 1, ref + pad, alt + pad


def a1_reads(record):
    """A1 reads of a graph record (FORMAT GT:DP:AD:AF, AD = REF,A1)."""
    fields = dict(zip(record[8].split(":"), record[9].split(":")))
    return int(fields["AD"].split(",")[1])


def graph_path(text):
    return [tuple(int(v) for v in step.split(":")) for step in text.split(",")] if text else []


def convert(args):
    import_pipeline()
    from indexed_gam_pipeline_v4.tensor_postprocessing.reference_path import ReferencePath

    path, genome = ReferencePath(args.reference_path), Genome(args.fasta)
    header, columns, records = read_vcf(args.graph_vcf)
    if meta(header, "source") != "pansoma_net_v2.graph_vcf":
        raise SystemExit(f"{args.graph_vcf} is not a graph_vcf output")
    unplaced_out = args.unplaced or str(args.output)[:-len(".vcf.gz")] + ".unplaced.vcf.gz"
    placed, unplaced = {}, []
    stats = dict(graph_records=0, unplaced=0, unplaced_pass=0, left_shifted=0, merged=0)
    for record in records:
        stats["graph_records"] += 1
        info = parse_info(record[7])
        kind = info["KIND"]
        ref, alt = (info.get("EVREF") or "").upper(), (info.get("EVALT") or "").upper()
        lin = path.linear(int(record[0]), int(info["START"]), ref, alt, kind, graph_path(info.get("PATH")))
        if lin is None:
            unplaced.append(record)
            stats["unplaced"] += 1
            stats["unplaced_pass"] += record[6] == "PASS"
            continue
        chrom, pos0, ref, alt = lin["chrom"], lin["pos0"], lin["ref"], lin["alt"]
        if chrom not in genome.lengths:
            raise SystemExit(f"{record[2]}: contig {chrom} is not in {args.fasta}")
        bases = genome.get(chrom, pos0, pos0 + len(ref))
        if bases != ref:
            raise SystemExit(f"{record[2]}: GRCh38 {chrom}:{pos0 + 1} REF {ref} but the FASTA has {bases}; "
                             f"the reference path and the FASTA disagree")
        moved = 0
        if kind in ("INS", "DEL"):
            pos0, ref, alt, moved = genome.left_align(chrom, pos0, ref, alt, kind)
            stats["left_shifted"] += moved > 0
        pos, vref, valt = genome.vcf_allele(chrom, pos0, ref, alt, kind)
        new_info = dict(NODE=record[0], **info, LEFTSHIFT=moved or None)
        line = [chrom, str(pos), record[2], vref, valt, record[5], record[6], None, *record[8:]]
        key, rank, alt_reads = (chrom, pos, vref, valt), (record[6] == "PASS", float(info["P_SOMATIC"])), a1_reads(record)
        if key in placed:
            stats["merged"] += 1
            kept, ids, kept_rank, reads = placed[key]
            if rank > kept_rank:
                placed[key] = ((line, new_info), [record[2]] + ids, rank, reads + alt_reads)
            else:
                ids.append(record[2])
                placed[key] = (kept, ids, kept_rank, reads + alt_reads)
        else:
            placed[key] = ((line, new_info), [record[2]], rank, alt_reads)
    order = {name: k for k, name in enumerate(genome.fasta.references)}
    out = []
    for key in sorted(placed, key=lambda k: (order[k[0]], k[1], k[2], k[3])):
        (line, info), ids, _, reads = placed[key]
        if len(ids) > 1:
            info["MERGED"], info["MERGED_ALT"] = len(ids), reads
        line[2], line[7] = ";".join(ids), format_info(info)
        out.append(line)

    kept_header = [h for h in header if not h.startswith(("##contig=", "##fileformat", "##source", "##pansoma_coordinates"))]
    info_end = max(k for k, h in enumerate(kept_header) if h.startswith(("##INFO", "##FILTER", "##FORMAT")))
    new_header = (["##fileformat=VCFv4.2", "##source=pansoma_net_v2.linear_vcf",
                   "##pansoma_coordinates=GRCh38", f"##pansoma_graph_vcf={Path(args.graph_vcf).resolve()}",
                   f"##pansoma_reference_path={Path(args.reference_path).resolve()}",
                   f"##reference=file://{Path(args.fasta).resolve()}"]
                  + kept_header[:info_end + 1] + EXTRA_INFO.split("\n") + kept_header[info_end + 1:]
                  + [f"##contig=<ID={n},length={genome.lengths[n]}>" for n in genome.fasta.references])
    write_vcf(args.output, new_header, columns, out)
    write_vcf(unplaced_out, header, columns, unplaced)
    filters = {}
    for line in out:
        for name in line[6].split(";"):
            filters[name] = filters.get(name, 0) + 1
    stats.update(records=len(out), filters=filters, output=str(Path(args.output).resolve()),
                 unplaced_output=str(Path(unplaced_out).resolve()))
    stats_path(args.output).write_text(json.dumps(stats, indent=2) + "\n")
    print(f"{args.output}: {stats['graph_records']:,} graph records -> {len(out):,} GRCh38 records "
          f"({stats['merged']:,} merged, {stats['left_shifted']:,} left-aligned); {stats['unplaced']:,} without "
          f"a GRCh38 position ({stats['unplaced_pass']:,} PASS) in {unplaced_out}", flush=True)
    return stats


def parse_args(argv=None):
    p = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    p.add_argument("--graph-vcf", required=True)
    p.add_argument("--reference-path", required=True, help="reference path directory of tools.graph_prep ref-path-scan")
    p.add_argument("--fasta", required=True, help="GRCh38 FASTA with .fai (the one the reference path was checked with)")
    p.add_argument("--output", required=True, help="<name>.linear.vcf.gz")
    p.add_argument("--unplaced", help="default: <output without .vcf.gz>.unplaced.vcf.gz")
    return p.parse_args(argv)


def main(argv=None):
    convert(parse_args(argv))


if __name__ == "__main__":
    main()
