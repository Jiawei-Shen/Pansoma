#!/usr/bin/env python3

import argparse
import time
from collections import Counter
import pysam


def log(msg):
    print(f"[{time.strftime('%Y-%m-%d %H:%M:%S')}] {msg}", flush=True)


def write_header(out_txt, args):
    with open(out_txt, "w") as f:
        f.write("DeepVariant-like lightweight candidate-site summary\n")
        f.write("=" * 80 + "\n")
        f.write(f"Input BAM: {args.bam}\n")
        f.write(f"Input FASTA: {args.fasta}\n")
        f.write(f"Minimum AF: {args.min_af}\n")
        f.write(f"Minimum count: {args.min_count}\n")
        f.write(f"Minimum BQ: {args.min_bq}\n")
        f.write(f"Minimum MQ: {args.min_mq}\n")
        f.write(f"Autosomes only: {args.autosomes_only}\n")
        f.write("=" * 80 + "\n\n")
        f.write(
            "chrom\tchrom_length\tpileup_sites\tSNV_candidate_sites\t"
            "INS_candidate_sites\tDEL_candidate_sites\tINDEL_candidate_sites\t"
            "ANY_candidate_sites\telapsed_minutes\n"
        )


def append_chrom_summary(
    out_txt,
    chrom,
    chrom_len,
    pileup_sites,
    snv_sites,
    ins_sites,
    del_sites,
    any_sites,
    elapsed_min,
):
    indel_sites = ins_sites + del_sites

    with open(out_txt, "a") as f:
        f.write(
            f"{chrom}\t{chrom_len}\t{pileup_sites}\t"
            f"{snv_sites}\t{ins_sites}\t{del_sites}\t"
            f"{indel_sites}\t{any_sites}\t{elapsed_min:.2f}\n"
        )
        f.flush()


def append_final_summary(
    out_txt,
    total_pileup,
    total_snv,
    total_ins,
    total_del,
    total_any,
    elapsed_min,
):
    with open(out_txt, "a") as f:
        f.write("\n")
        f.write("=" * 80 + "\n")
        f.write("Final total summary\n")
        f.write("=" * 80 + "\n")
        f.write(f"Total pileup sites: {total_pileup}\n")
        f.write(f"Total SNV candidate sites: {total_snv}\n")
        f.write(f"Total INS candidate sites: {total_ins}\n")
        f.write(f"Total DEL candidate sites: {total_del}\n")
        f.write(f"Total INDEL candidate sites: {total_ins + total_del}\n")
        f.write(f"Total ANY candidate sites: {total_any}\n")
        f.write(f"Total elapsed minutes: {elapsed_min:.2f}\n")


def get_chrom_list(bam, autosomes_only):
    bam_refs = set(bam.references)

    if autosomes_only:
        chroms = [f"chr{i}" for i in range(1, 23)]
        return [c for c in chroms if c in bam_refs]

    return list(bam.references)


def detect_candidate_sites(
    bam_path,
    fasta_path,
    out_txt,
    min_af,
    min_count,
    min_bq,
    min_mq,
    log_interval,
    autosomes_only,
):
    log("Starting lightweight candidate-site detection")
    log(f"BAM: {bam_path}")
    log(f"FASTA: {fasta_path}")
    log(f"Output summary: {out_txt}")

    bam = pysam.AlignmentFile(bam_path, "rb")
    fasta = pysam.FastaFile(fasta_path)

    chrom_list = get_chrom_list(bam, autosomes_only)

    total_pileup = 0
    total_snv = 0
    total_ins = 0
    total_del = 0
    total_any = 0

    start_time = time.time()

    for chrom in chrom_list:
        chrom_start = time.time()
        chrom_len = bam.get_reference_length(chrom)

        log("=" * 80)
        log(f"Processing {chrom}, length={chrom_len:,}")

        log(f"Loading FASTA sequence for {chrom}")
        chrom_seq = fasta.fetch(chrom).upper()

        pileup_sites = 0
        snv_sites = 0
        ins_sites = 0
        del_sites = 0
        any_sites = 0

        last_logged_pos = 0

        for col in bam.pileup(
            chrom,
            0,
            chrom_len,
            truncate=True,
            stepper="samtools",
            min_base_quality=min_bq,
            min_mapping_quality=min_mq,
        ):
            pos0 = col.reference_pos
            pos1 = pos0 + 1

            pileup_sites += 1
            total_pileup += 1

            if pos1 - last_logged_pos >= log_interval:
                pct = pos1 / chrom_len * 100
                log(
                    f"{chrom}: {pos1:,}/{chrom_len:,} "
                    f"({pct:.2f}%), scanned={pileup_sites:,}, "
                    f"SNV={snv_sites:,}, INS={ins_sites:,}, DEL={del_sites:,}, ANY={any_sites:,}"
                )
                last_logged_pos = pos1

            ref_base = chrom_seq[pos0]
            if ref_base not in {"A", "C", "G", "T"}:
                continue

            depth = 0
            mismatch_count = 0
            ins_count = 0
            del_count = 0

            for pr in col.pileups:
                aln = pr.alignment

                if (
                    aln.is_unmapped
                    or aln.is_duplicate
                    or aln.is_secondary
                    or aln.is_supplementary
                ):
                    continue

                if aln.mapping_quality < min_mq:
                    continue

                # Count SNV-like evidence
                if not pr.is_del and not pr.is_refskip and pr.query_position is not None:
                    qpos = pr.query_position

                    if aln.query_qualities and aln.query_qualities[qpos] < min_bq:
                        continue

                    base = aln.query_sequence[qpos].upper()

                    if base in {"A", "C", "G", "T"}:
                        depth += 1
                        if base != ref_base:
                            mismatch_count += 1

                # Count INDEL-like evidence
                if pr.indel > 0:
                    ins_count += 1
                elif pr.indel < 0:
                    del_count += 1

            if depth == 0:
                continue

            snv_af = mismatch_count / depth
            ins_af = ins_count / depth
            del_af = del_count / depth

            is_snv = mismatch_count > min_count and snv_af > min_af
            is_ins = ins_count > min_count and ins_af > min_af
            is_del = del_count > min_count and del_af > min_af

            if is_snv:
                snv_sites += 1
                total_snv += 1

            if is_ins:
                ins_sites += 1
                total_ins += 1

            if is_del:
                del_sites += 1
                total_del += 1

            if is_snv or is_ins or is_del:
                any_sites += 1
                total_any += 1

        elapsed_min = (time.time() - chrom_start) / 60

        append_chrom_summary(
            out_txt=out_txt,
            chrom=chrom,
            chrom_len=chrom_len,
            pileup_sites=pileup_sites,
            snv_sites=snv_sites,
            ins_sites=ins_sites,
            del_sites=del_sites,
            any_sites=any_sites,
            elapsed_min=elapsed_min,
        )

        log(
            f"Finished {chrom}: scanned={pileup_sites:,}, "
            f"SNV={snv_sites:,}, INS={ins_sites:,}, DEL={del_sites:,}, "
            f"ANY={any_sites:,}, elapsed={elapsed_min:.2f} min"
        )
        log(f"Updated summary file: {out_txt}")

    bam.close()
    fasta.close()

    total_elapsed = (time.time() - start_time) / 60

    append_final_summary(
        out_txt=out_txt,
        total_pileup=total_pileup,
        total_snv=total_snv,
        total_ins=total_ins,
        total_del=total_del,
        total_any=total_any,
        elapsed_min=total_elapsed,
    )

    log("=" * 80)
    log("Finished all chromosomes")
    log(f"Total SNV candidate sites: {total_snv:,}")
    log(f"Total INS candidate sites: {total_ins:,}")
    log(f"Total DEL candidate sites: {total_del:,}")
    log(f"Total ANY candidate sites: {total_any:,}")
    log(f"Final summary written to: {out_txt}")


def main():
    parser = argparse.ArgumentParser(
        description="Fast DeepVariant-like candidate-site counting from BAM."
    )

    parser.add_argument("-b", "--bam", required=True)
    parser.add_argument("-f", "--fasta", required=True)
    parser.add_argument("-o", "--out", required=True)

    parser.add_argument("--min-af", type=float, default=0.1)
    parser.add_argument("--min-count", type=int, default=2)
    parser.add_argument("--min-bq", type=int, default=13)
    parser.add_argument("--min-mq", type=int, default=20)

    parser.add_argument(
        "--log-interval",
        type=int,
        default=1_000_000,
    )

    parser.add_argument(
        "--autosomes-only",
        action="store_true",
        help="Only process chr1-chr22",
    )

    args = parser.parse_args()

    write_header(args.out, args)

    detect_candidate_sites(
        bam_path=args.bam,
        fasta_path=args.fasta,
        out_txt=args.out,
        min_af=args.min_af,
        min_count=args.min_count,
        min_bq=args.min_bq,
        min_mq=args.min_mq,
        log_interval=args.log_interval,
        autosomes_only=args.autosomes_only,
    )


if __name__ == "__main__":
    main()