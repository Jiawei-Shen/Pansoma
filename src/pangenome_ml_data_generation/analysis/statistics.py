#!/usr/bin/env python3

import argparse
import json
import math
import subprocess
import statistics
from collections import Counter


def run_cmd(cmd):
    return subprocess.Popen(cmd, stdout=subprocess.PIPE, text=True)


def summarize_gam(gam_path):
    """
    Requires: vg
    Uses: vg view -aj input.gam
    """
    stats = Counter()
    mapqs = []

    p = run_cmd(["vg", "view", "-aj", gam_path])

    for line in p.stdout:
        if not line.strip():
            continue

        aln = json.loads(line)

        stats["total_alignments"] += 1

        mapq = aln.get("mapping_quality")
        if mapq is not None:
            mapqs.append(mapq)

        path = aln.get("path", {})
        mappings = path.get("mapping", [])

        if mappings:
            stats["total_aligned"] += 1

        perfect = True
        gapless_softclips_allowed = True

        for m in mappings:
            for e in m.get("edit", []):
                from_len = e.get("from_length", 0)
                to_len = e.get("to_length", 0)
                seq = e.get("sequence", "")

                # Match or mismatch
                if from_len > 0 and to_len > 0:
                    if seq:
                        # substitution block
                        stats["substitutions_bp"] += to_len
                        perfect = False
                        gapless_softclips_allowed = False
                    else:
                        # exact match
                        stats["perfect_match_bp"] += to_len

                # Insertion relative to graph/reference
                elif from_len == 0 and to_len > 0:
                    stats["insertions_bp"] += to_len
                    stats["softclips_bp"] += to_len
                    perfect = False
                    # insertion/softclip allowed for gapless-softclip metric

                # Deletion relative to graph/reference
                elif from_len > 0 and to_len == 0:
                    stats["deletions_bp"] += from_len
                    perfect = False
                    gapless_softclips_allowed = False

        if perfect:
            stats["total_perfect"] += 1

        if mappings and gapless_softclips_allowed:
            stats["total_gapless_softclips_allowed"] += 1

    p.wait()

    add_mapq_summary(stats, mapqs)
    return stats


def parse_cigar(cigar):
    """
    Parse BAM CIGAR string.
    Example: 10M2I5M3D
    """
    ops = []
    num = ""

    for c in cigar:
        if c.isdigit():
            num += c
        else:
            if num:
                ops.append((int(num), c))
                num = ""

    return ops


def summarize_bam(bam_path):
    """
    Requires: samtools
    Uses: samtools view input.bam

    Notes:
    - Perfect alignment uses NM:i:0 if available.
    - Substitutions require NM plus indel counts.
    - If NM is missing, substitution counts are approximate/incomplete.
    """
    stats = Counter()
    mapqs = []

    p = run_cmd(["samtools", "view", bam_path])

    for line in p.stdout:
        fields = line.rstrip("\n").split("\t")
        if len(fields) < 11:
            continue

        flag = int(fields[1])
        cigar = fields[5]
        mapq = int(fields[4])

        # Skip unmapped reads
        if flag & 0x4:
            continue

        # Optional: skip secondary/supplementary
        if flag & 0x100 or flag & 0x800:
            continue

        stats["total_alignments"] += 1
        stats["total_aligned"] += 1
        mapqs.append(mapq)

        tags = fields[11:]
        nm = None

        for tag in tags:
            if tag.startswith("NM:i:"):
                nm = int(tag.split(":")[-1])
                break

        insertion_bp = 0
        deletion_bp = 0
        softclip_bp = 0
        match_like_bp = 0

        has_gap = False

        for length, op in parse_cigar(cigar):
            if op in {"M", "=", "X"}:
                match_like_bp += length
                if op == "X":
                    stats["substitutions_bp"] += length
            elif op == "I":
                insertion_bp += length
                has_gap = True
            elif op == "D":
                deletion_bp += length
                has_gap = True
            elif op == "S":
                softclip_bp += length

        stats["insertions_bp"] += insertion_bp
        stats["deletions_bp"] += deletion_bp
        stats["softclips_bp"] += softclip_bp

        if nm is not None:
            substitutions = max(nm - insertion_bp - deletion_bp, 0)
            stats["substitutions_bp"] += substitutions

            if nm == 0 and softclip_bp == 0:
                stats["total_perfect"] += 1

        if not has_gap:
            stats["total_gapless_softclips_allowed"] += 1

    p.wait()

    add_mapq_summary(stats, mapqs)
    return stats


def add_mapq_summary(stats, mapqs):
    if mapqs:
        stats["mapq_mean"] = statistics.mean(mapqs)
        stats["mapq_sd"] = statistics.stdev(mapqs) if len(mapqs) > 1 else 0
        stats["mapq_median"] = statistics.median(mapqs)
        stats["mapq_n"] = len(mapqs)
    else:
        stats["mapq_mean"] = "NA"
        stats["mapq_sd"] = "NA"
        stats["mapq_median"] = "NA"
        stats["mapq_n"] = 0


def write_row(out, sample, platform, source, stats):
    keys = [
        "total_alignments",
        "total_aligned",
        "total_perfect",
        "total_gapless_softclips_allowed",
        "softclips_bp",
        "substitutions_bp",
        "deletions_bp",
        "insertions_bp",
        "mapq_mean",
        "mapq_sd",
        "mapq_median",
        "mapq_n",
    ]

    values = [sample, platform, source] + [stats.get(k, 0) for k in keys]
    out.write("\t".join(map(str, values)) + "\n")


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--sample", default="HG008T")
    parser.add_argument("--pacbio-gam")
    parser.add_argument("--pacbio-bam")
    parser.add_argument("--ont-gam")
    parser.add_argument("--ont-bam")
    parser.add_argument("-o", "--output", required=True)

    args = parser.parse_args()

    with open(args.output, "w") as out:
        header = [
            "sample",
            "platform",
            "source",
            "total_alignments",
            "total_aligned",
            "total_perfect",
            "total_gapless_softclips_allowed",
            "softclips_bp",
            "substitutions_bp",
            "deletions_bp",
            "insertions_bp",
            "mapq_mean",
            "mapq_sd",
            "mapq_median",
            "mapq_n",
        ]
        out.write("\t".join(header) + "\n")

        if args.pacbio_gam:
            write_row(out, args.sample, "PacBio", "Graph", summarize_gam(args.pacbio_gam))

        if args.pacbio_bam:
            write_row(out, args.sample, "PacBio", "Linear", summarize_bam(args.pacbio_bam))

        if args.ont_gam:
            write_row(out, args.sample, "ONT", "Graph", summarize_gam(args.ont_gam))

        if args.ont_bam:
            write_row(out, args.sample, "ONT", "Linear", summarize_bam(args.ont_bam))


if __name__ == "__main__":
    main()