#!/bin/bash
# Every INDEL call set on RIS re-scored with the INDEL PoN rule of 71cb51a (gnomAD / CoLoRSdb exact allele at AF >= 0.01
# when the ALT changes a homopolymer of >= 8 bp, else AF >= 1e-4; was >= 7 bp and 1e-3; SNV rule unchanged): the PoN +
# rtg vcfeval steps only (ris/vcfeval_only.sh, this checkout's scripts/filter_panel_of_normals.py), on each run's
# existing GRCh38 calls <run>/<SRC>/*.INDEL.linear.vcf.gz, into <run>/<SRC>_pon71cb51a/eval_bed and eval_nobed. The
# calls and their PASS threshold are unchanged. Run from the repository root on the login node: bash <this file>
# The list: every run with args.json kinds == [INDEL] and an evaluated <SRC> = vcf_chr1* (not an earlier *_pon* re-score),
# plus the two COLO829T ONT INDEL no-ch6 runs still training on 2026-10-16 (after their own vcf.sh).
set -euo pipefail
cd "$(git rev-parse --show-toplevel)"
RUNS=/storage2/fs1/epigenome/Active/jiawei.s/pansoma_net_v2_runs
list=()
for d in "$RUNS"/*/; do
  r=$(basename "$d")
  [ -f "$d/args.json" ] || continue
  python3 -c "import json, sys; sys.exit(0 if json.load(open(sys.argv[1])).get('kinds') == ['INDEL'] else 1)" "$d/args.json" || continue
  for v in "$d"vcf_chr1*/; do
    s=$(basename "$v")
    case $s in *_pon*) continue;; esac
    [ -f "$d/${s}_pon71cb51a/eval_bed/report.json" ] && continue  # done (HG008_Illumina_INDEL_7fc7396: the test job)
    [ -f "$v/eval_bed/report.json" ] && ls "$v"/*.INDEL.linear.vcf.gz >/dev/null 2>&1 && list+=("$r:$s")
  done
done
for rs in "${list[@]}"; do
  r=${rs%%:*} s=${rs#*:}
  ris/submit.sh -g 0 -n 4 -m 32 "pnv2_rescore_${r}_${s}" vcfeval_only.sh "$r" INDEL "$s" "${s}_pon71cb51a"
done
for r in COLO829_ONT_INDEL_nopc_7fc7396 COLO829_ONT_INDEL_nopc_d46_7fc7396; do
  if [[ " ${list[*]} " != *" $r:vcf_chr1 "* ]]; then
    ris/submit.sh -g 0 -n 4 -m 32 -w "done(pnv2_vcf_$r)" "pnv2_rescore_${r}_vcf_chr1" vcfeval_only.sh "$r" INDEL vcf_chr1 \
        vcf_chr1_pon71cb51a
  fi
done
echo "${#list[@]} re-scores submitted now"
