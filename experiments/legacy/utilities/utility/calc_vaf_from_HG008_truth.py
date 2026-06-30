#!/usr/bin/env python3
import pysam
import sys

BAM_PATH = "/scratch/jshen/data/HG008_GIAB/raw_sequencing_data/Liss_lab_PacBio_Revio_20240125/HG008-T_PacBio-HiFi-Revio_20240125_116x_GRCh38-GIABv3.bam"
VCF_PATH = "/scratch/jshen/data/HG008_GIAB/draft_v02_benchmark/HG008-T_somatic_smvar_benchmark_v0.2_tumorvariants.vcf.gz"
OUT_PATH = "/scratch/jshen/tmp/HG008-T_somatic_smvar_benchmark_v0.2_tumorvariants.vaf.tsv"

MIN_MAPQ = 0
MIN_BASEQ = 0
PASS_ONLY = False


def infer_variant_type(ref: str, alt: str) -> str:
    if len(ref) == 1 and len(alt) == 1:
        return "SNV"
    if len(ref) < len(alt) and alt.startswith(ref):
        return "INS"
    if len(ref) > len(alt) and ref.startswith(alt):
        return "DEL"
    return "COMPLEX"


def count_snv(bam, chrom, pos1, ref, alt, min_mapq=0, min_baseq=0):
    depth = 0
    ref_count = 0
    alt_count = 0
    pos0 = pos1 - 1

    for pileupcolumn in bam.pileup(
        chrom,
        pos0,
        pos0 + 1,
        truncate=True,
        stepper="all",
        min_base_quality=min_baseq,
        ignore_overlaps=False,
        ignore_orphans=False,
    ):
        if pileupcolumn.reference_pos != pos0:
            continue

        for pileupread in pileupcolumn.pileups:
            aln = pileupread.alignment

            if aln.is_unmapped or aln.is_duplicate or aln.is_qcfail:
                continue
            if aln.mapping_quality < min_mapq:
                continue
            if pileupread.is_del or pileupread.is_refskip:
                continue
            if pileupread.query_position is None:
                continue

            qpos = pileupread.query_position
            base = aln.query_sequence[qpos]
            bq = aln.query_qualities[qpos] if aln.query_qualities is not None else 0
            if bq < min_baseq:
                continue

            depth += 1
            if base.upper() == ref.upper():
                ref_count += 1
            elif base.upper() == alt.upper():
                alt_count += 1

    return depth, ref_count, alt_count


def count_insertion(bam, chrom, pos1, ref, alt, min_mapq=0, min_baseq=0):
    ins_seq = alt[len(ref):]
    ins_len = len(ins_seq)

    depth = 0
    ref_count = 0
    alt_count = 0
    pos0 = pos1 - 1

    for pileupcolumn in bam.pileup(
        chrom,
        pos0,
        pos0 + 1,
        truncate=True,
        stepper="all",
        min_base_quality=min_baseq,
        ignore_overlaps=False,
        ignore_orphans=False,
    ):
        if pileupcolumn.reference_pos != pos0:
            continue

        for pileupread in pileupcolumn.pileups:
            aln = pileupread.alignment

            if aln.is_unmapped or aln.is_duplicate or aln.is_qcfail:
                continue
            if aln.mapping_quality < min_mapq:
                continue
            if pileupread.is_del or pileupread.is_refskip:
                continue
            if pileupread.query_position is None:
                continue

            qpos = pileupread.query_position
            bq = aln.query_qualities[qpos] if aln.query_qualities is not None else 0
            if bq < min_baseq:
                continue

            depth += 1

            if pileupread.indel == ins_len:
                qseq = aln.query_sequence
                inserted = qseq[qpos + 1:qpos + 1 + ins_len]
                if inserted.upper() == ins_seq.upper():
                    alt_count += 1
            elif pileupread.indel == 0:
                ref_count += 1

    return depth, ref_count, alt_count


def count_deletion(bam, chrom, pos1, ref, alt, min_mapq=0, min_baseq=0):
    del_len = len(ref) - len(alt)

    depth = 0
    ref_count = 0
    alt_count = 0
    pos0 = pos1 - 1

    for pileupcolumn in bam.pileup(
        chrom,
        pos0,
        pos0 + 1,
        truncate=True,
        stepper="all",
        min_base_quality=min_baseq,
        ignore_overlaps=False,
        ignore_orphans=False,
    ):
        if pileupcolumn.reference_pos != pos0:
            continue

        for pileupread in pileupcolumn.pileups:
            aln = pileupread.alignment

            if aln.is_unmapped or aln.is_duplicate or aln.is_qcfail:
                continue
            if aln.mapping_quality < min_mapq:
                continue
            if pileupread.is_del or pileupread.is_refskip:
                continue
            if pileupread.query_position is None:
                continue

            qpos = pileupread.query_position
            bq = aln.query_qualities[qpos] if aln.query_qualities is not None else 0
            if bq < min_baseq:
                continue

            depth += 1

            if pileupread.indel == -del_len:
                alt_count += 1
            elif pileupread.indel == 0:
                ref_count += 1

    return depth, ref_count, alt_count


def calc_variant_counts(bam, chrom, pos1, ref, alt, min_mapq, min_baseq):
    vtype = infer_variant_type(ref, alt)

    if vtype == "SNV":
        depth, ref_count, alt_count = count_snv(
            bam, chrom, pos1, ref, alt, min_mapq=min_mapq, min_baseq=min_baseq
        )
    elif vtype == "INS":
        depth, ref_count, alt_count = count_insertion(
            bam, chrom, pos1, ref, alt, min_mapq=min_mapq, min_baseq=min_baseq
        )
    elif vtype == "DEL":
        depth, ref_count, alt_count = count_deletion(
            bam, chrom, pos1, ref, alt, min_mapq=min_mapq, min_baseq=min_baseq
        )
    else:
        depth, ref_count, alt_count = 0, 0, 0

    return vtype, depth, ref_count, alt_count


def main():
    print(f"[INFO] BAM: {BAM_PATH}", file=sys.stderr)
    print(f"[INFO] VCF: {VCF_PATH}", file=sys.stderr)
    print(f"[INFO] OUT: {OUT_PATH}", file=sys.stderr)

    bam = pysam.AlignmentFile(BAM_PATH, "rb")
    vcf = pysam.VariantFile(VCF_PATH)

    with open(OUT_PATH, "w") as fout:
        fout.write(
            "\t".join([
                "chrom",
                "pos",
                "ref",
                "alt",
                "type",
                "depth",
                "ref_count",
                "alt_count",
                "other_count",
                "vaf",
                "filter"
            ]) + "\n"
        )

        n = 0
        for rec in vcf:
            if PASS_ONLY:
                filt_keys = list(rec.filter.keys())
                if len(filt_keys) > 0 and "PASS" not in filt_keys:
                    continue

            if rec.alts is None or len(rec.alts) == 0:
                continue

            chrom = rec.chrom
            pos = rec.pos
            ref = rec.ref
            alt = rec.alts[0]

            vtype, depth, ref_count, alt_count = calc_variant_counts(
                bam, chrom, pos, ref, alt, MIN_MAPQ, MIN_BASEQ
            )

            other_count = depth - ref_count - alt_count

            if depth > 0:
                vaf = alt_count / depth
                vaf_str = f"{vaf:.6f}"
            else:
                vaf_str = "NA"

            filt_str = ";".join(rec.filter.keys()) if len(rec.filter.keys()) > 0 else "."

            if vtype == "COMPLEX":
                vtype = "UNSUPPORTED_COMPLEX"

            fout.write(
                "\t".join([
                    chrom,
                    str(pos),
                    ref,
                    alt,
                    vtype,
                    str(depth),
                    str(ref_count),
                    str(alt_count),
                    str(other_count),
                    vaf_str,
                    filt_str
                ]) + "\n"
            )

            n += 1
            if n % 100 == 0:
                print(f"[INFO] processed {n} variants", file=sys.stderr)

    bam.close()
    vcf.close()
    print(f"[INFO] done, wrote: {OUT_PATH}", file=sys.stderr)


if __name__ == "__main__":
    main()