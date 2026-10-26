#!/bin/bash
# The alleles of the mgi_non_HG005.d9 graph on GRCh38 (vg deconstruct VCF of the d9 GBZ, every snarl level): contigs
# GRCh38#0#chrN -> chrN, multi-allelic records split, left-normalized against the GRCh38 no-alt FASTA, symbolic / '*'
# alleles dropped; written as CHROM POS REF ALT (gzip). A truth allele found here is a path of the graph, so reads that
# carry it can align without an edit (absorbed): used by tensor_vaf.py --graph-alleles.
# ATOMIZE=1: also split multi-base alleles (MNPs, complex) of <= 50 bp into their SNVs and indels (bcftools norm -a), so
# an SNV carried inside a larger allele of a bubble counts as a graph allele (an upper bound on absorption: a read with
# only that SNV still has an edit against such a path). Longer alleles are kept whole (atomizing the 1.6 Mb chr5 SMN
# allele alone exceeds 64 GB).
# usage: [ATOMIZE=1] sbatch -c 4 --mem=16G graph_alleles.sh OUT.tsv.gz
set -euo pipefail
OUT=$1
ml load bcftools
VCF=/scratch/jshen/data/non_HG005_VG_Indexes/mgi_non_HG005.d9.vcf.gz
FASTA=/scratch/jshen/data/HapMap/GCA_000001405.15_GRCh38_no_alt_analysis_set.fasta
TMP=$(mktemp -d)
for c in $(seq 1 22) X Y M; do printf 'GRCh38#0#chr%s\tchr%s\n' $c $c >> $TMP/map; done
bcftools view -G -r $(seq -s, 1 22 | sed 's/\([0-9]*\)/GRCh38#0#chr\1/g') $VCF -Ou \
  | bcftools annotate --rename-chrs $TMP/map -Ou \
  | bcftools norm -m -any -f $FASTA -c s -Ob -o $TMP/split.bcf
KEEP='ALT!="*" && ALT!="." && ALT!~"<"'
FMT='%CHROM\t%POS\t%REF\t%ALT\n'
if [ -z "${ATOMIZE:-}" ]; then
  bcftools query -i "$KEEP" -f "$FMT" $TMP/split.bcf | gzip -1 > $OUT
else
  SHORT='strlen(REF)<=50 && strlen(ALT)<=50'
  { bcftools query -i "$KEEP" -f "$FMT" $TMP/split.bcf
    bcftools view -i "$SHORT" $TMP/split.bcf -Ou | bcftools norm -a --atom-overlaps . -f $FASTA -c s -Ou \
      | bcftools query -i "$KEEP" -f "$FMT"; } | LC_ALL=C sort -u -S 4G -T $TMP | gzip -1 > $OUT
fi
rm -f $TMP/map $TMP/split.bcf; rmdir $TMP
zcat $OUT | awk '{n++; if (length($3)==1 && length($4)==1) s++} END {print n " graph alleles, " s " SNVs"}'
