"""scripts/filter_panel_of_normals.py with another INDEL AF floor (ris experiment; the script is not changed).

    python ris/tools/pon_tag.py FLOOR INPUT.vcf.gz OUTPUT.vcf.gz GNOMAD DBSNP 1000G COLORSDB

The INDEL rule becomes gnomAD / CoLoRSdb exact allele at population AF >= FLOOR for every INDEL, one flat floor (the
script's own INDEL rule is 0.01 in homopolymers >= 7 bp and 0.001 otherwise); the SNV rule is unchanged. The output is
the tagged VCF that vcfeval --pon-vcf takes."""
import sys

import filter_panel_of_normals as F

floor = float(sys.argv[1])
F.INDEL_MIN_AF = floor  # only the --help text: the dict rule below replaces the script's homopolymer split
F.PON_RULES["INDEL"] = (dict(min_af=floor), None, None, dict(min_af=floor))
sys.argv = [F.__file__, sys.argv[2], sys.argv[3], "--pon", *sys.argv[4:8]]
sys.exit(F.main())
