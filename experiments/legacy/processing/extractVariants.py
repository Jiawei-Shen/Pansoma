import pysam
import json
import time
import argparse


def read_vcf_pysam(file_path, chromosome=None):
    vcf_file = pysam.VariantFile(file_path)
    variants = []

    for record in vcf_file:
        if chromosome and record.chrom != chromosome:
            continue  # Skip records not on the specified chromosome

        variant_data = {
            "chromosome": record.chrom,
            "position": record.pos,
            "reference": record.ref,
            "alternative": record.alts,
            "id": record.id,
            "quality": record.qual,
            "samples": {sample: record.samples[sample].alleles for sample in record.samples}
        }
        variants.append(variant_data)

    return variants


def format_chromosome(chromosome):
    """Ensure chromosome is formatted with 'chr' prefix."""
    if not chromosome.startswith('chr'):
        return f'chr{chromosome}'
    return chromosome


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Read VCF file and filter by chromosome.")
    parser.add_argument("vcf_path", help="Path to the VCF file.")
    parser.add_argument("--chromosome", help="Chromosome to filter by (e.g., 'chr1', '1').", default=None)
    args = parser.parse_args()

    # Format chromosome if provided
    chromosome = format_chromosome(args.chromosome) if args.chromosome else None

    start_time = time.time()  # Start time

    vcf_data = read_vcf_pysam(args.vcf_path, chromosome)

    # Generate output filename
    output_filename = f"{args.vcf_path.split('/')[-1].split('.')[0]}_variants"
    if chromosome:
        output_filename += f"_{chromosome}"
    output_filename += ".json"

    # Save to JSON file
    with open(output_filename, "w") as f:
        json.dump(vcf_data, f, indent=4)

    # Print first 5 entries
    # print(vcf_data[:5])

    end_time = time.time()  # End time
    elapsed_time = end_time - start_time
    print(f"Elapsed time: {elapsed_time:.6f} seconds")
