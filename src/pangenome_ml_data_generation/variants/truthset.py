import pysam
import json
import argparse


def extract_pileup_data(bam_file, vcf_file, output_json, flanking_bp=200):
    # Open BAM file
    bam = pysam.AlignmentFile(bam_file, "rb")

    # Read variants from VCF file
    variants = []
    with open(vcf_file, "r") as vcf:
        for line in vcf:
            if line.startswith("#"):
                continue  # Skip headers
            cols = line.strip().split("\t")
            chrom, pos, ref, alt = cols[0], int(cols[1]), cols[3], cols[4]
            variants.append((chrom, pos, ref, alt))

    # Dictionary to store pileup data
    pileup_data = {}

    for chrom, pos, ref, alt in variants:
        key = f"{chrom}:{pos}"
        pileup_data[key] = {
            "variant": {
                "chromosome": chrom,
                "position": pos,
                "reference": ref,
                "alternate": alt
            },
            "pileup_reads": []
        }

        # Fetch reads within the ±200 bp region
        for read in bam.fetch(chrom, max(0, pos - flanking_bp), pos + flanking_bp):
            read_data = {
                "query_name": read.query_name,
                "reference_start": read.reference_start,
                "cigar": read.cigarstring,
                "sequence": read.query_sequence,
                "mapping_quality": read.mapping_quality,
                "is_reverse": read.is_reverse,
            }
            pileup_data[key]["pileup_reads"].append(read_data)

    # Close BAM file
    bam.close()

    # Write JSON output
    with open(output_json, "w") as json_file:
        json.dump(pileup_data, json_file, indent=4)

    print(f"Pileup data saved to {output_json}")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Extract pileup reads within 200bp of variants and save to JSON.")
    parser.add_argument("-b", "--bam", required=True, help="Input BAM file")
    parser.add_argument("-v", "--vcf", required=True, help="Input VCF file")
    parser.add_argument("-o", "--output", required=True, help="Output JSON file")
    args = parser.parse_args()

    extract_pileup_data(args.bam, args.vcf, args.output)
