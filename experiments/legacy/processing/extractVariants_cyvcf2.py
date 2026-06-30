import gzip
import vcfpy
import json
import time


def read_vcf_vcfpy(file_path):
    with gzip.open(file_path, "rt") as f:
        reader = vcfpy.Reader.from_stream(f)
        variants = []

        for record in reader:
            variant_data = {
                "chromosome": record.CHROM,
                "position": record.POS,
                "reference": record.REF,
                "alternative": record.ALT,
                "samples": {call.sample: call.data for call in record.calls}  # Sample genotype data
            }
            variants.append(variant_data)

    return variants


if __name__ == "__main__":
    start_time = time.time()  # Start time

    vcf_path = "./Genome_Bottle_VCF/HG005_GRCh38_1_22_v4.2.1_benchmark.vcf.gz"
    vcf_data = read_vcf_vcfpy(vcf_path)

    end_time = time.time()  # End time

    elapsed_time = end_time - start_time
    print(f"Elapsed time: {elapsed_time:.6f} seconds")

    # Save to JSON file
    with open("HG005_variants_vcfpy.json", "w") as f:
        json.dump(vcf_data, f, indent=4)

    # Print first 5 entries
    print(vcf_data[:5])

