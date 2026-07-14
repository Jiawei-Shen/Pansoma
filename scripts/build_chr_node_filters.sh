#!/usr/bin/env bash
set -euo pipefail

export LC_ALL=C

GBZ="${GBZ:-/scratch/jshen/data/AF-Filtered_VG_Indexes/hprc-v1.1-mc-grch38.d9.gbz}"
GFA="${GFA:-/scratch/jshen/data/AF-Filtered_VG_Indexes/hprc-v1.1-mc-grch38.d9.gfa}"

OUTDIR="${OUTDIR:-/scratch/jshen/tmp/chr_component_vs_GRCh38_summary}"
CHROMOSOMES="${CHROMOSOMES:-$(seq 1 22)}"
mkdir -p "${OUTDIR}"

SUMMARY="${OUTDIR}/summary.tsv"

echo -e "chr\tcomponent_total_nodes\tGRCh38_path_total_nodes\tshared_nodes\tcomponent_not_on_GRCh38\tGRCh38_not_in_component\tpct_component_on_GRCh38\tpct_component_not_on_GRCh38\tpct_extra_vs_GRCh38" > "${SUMMARY}"

for c in ${CHROMOSOMES}; do
    if [[ "${c}" == chr* ]]; then
        chr="${c}"
    else
        chr="chr${c}"
    fi
    echo "Processing ${chr} ..."

    CHR_DIR="${OUTDIR}/${chr}"
    mkdir -p "${CHR_DIR}"

    COMPONENT_VG="${CHR_DIR}/${chr}.component.vg"
    COMPONENT_NODES_RAW="${CHR_DIR}/${chr}.component.nodes.raw.txt"
    COMPONENT_NODES="${CHR_DIR}/${chr}.component.nodes.txt"

    GRCH38_RAW="${CHR_DIR}/${chr}.GRCh38_path.nodes.raw.txt"
    GRCH38_NODES="${CHR_DIR}/${chr}.GRCh38_path.nodes.txt"

    SHARED="${CHR_DIR}/${chr}.shared.txt"
    COMPONENT_NOT_GRCH38="${CHR_DIR}/${chr}.component_not_on_GRCh38.txt"
    GRCH38_NOT_COMPONENT="${CHR_DIR}/${chr}.GRCh38_not_in_component.txt"

    # --------------------------------------------------
    # 1) Extract full connected component from GBZ
    # --------------------------------------------------
    vg chunk \
      -x "${GBZ}" \
      -p "GRCh38#0#${chr}" \
      -C > "${COMPONENT_VG}"

    vg view -j "${COMPONENT_VG}" \
      | jq -r '.node[].id' > "${COMPONENT_NODES_RAW}"

    awk 'NF{print $1}' "${COMPONENT_NODES_RAW}" | sort -u > "${COMPONENT_NODES}"

    # --------------------------------------------------
    # 2) Extract nodes on GRCh38#0#chr only from GFA
    #    Prefer W-lines, fallback to P-lines
    # --------------------------------------------------
    if ! awk -v target="${chr}" '
    $1=="W" && $2=="GRCh38" && $3=="0" && $4==target {
        walk=$NF
        gsub(/[<>]/, " ", walk)
        n=split(walk, a, /[[:space:]]+/)
        for(i=1;i<=n;i++) if(a[i]!="") print a[i]
        found=1
    }
    END{ exit(found ? 0 : 1) }
    ' "${GFA}" > "${GRCH38_RAW}"; then

        awk -v target="GRCh38#0#${chr}" '
        $1=="P" && $2==target {
            path=$3
            gsub(/,/, " ", path)
            n=split(path, a, /[[:space:]]+/)
            for(i=1;i<=n;i++){
                x=a[i]
                sub(/[+-]$/, "", x)
                if(x!="") print x
            }
            found=1
        }
        END{
            if(!found){
                print "ERROR: could not find W-line or P-line for " target > "/dev/stderr"
                exit 1
            }
        }
        ' "${GFA}" > "${GRCH38_RAW}"
    fi

    awk 'NF{print $1}' "${GRCH38_RAW}" | sort -u > "${GRCH38_NODES}"

    # --------------------------------------------------
    # 3) Check sort order and compare sets
    # --------------------------------------------------
    comm --check-order "${COMPONENT_NODES}" "${GRCH38_NODES}" >/dev/null

    comm -12 "${COMPONENT_NODES}" "${GRCH38_NODES}" > "${SHARED}"
    comm -23 "${COMPONENT_NODES}" "${GRCH38_NODES}" > "${COMPONENT_NOT_GRCH38}"
    comm -13 "${COMPONENT_NODES}" "${GRCH38_NODES}" > "${GRCH38_NOT_COMPONENT}"

    component_total=$(wc -l < "${COMPONENT_NODES}")
    grch38_total=$(wc -l < "${GRCH38_NODES}")
    shared_total=$(wc -l < "${SHARED}")
    component_not_total=$(wc -l < "${COMPONENT_NOT_GRCH38}")
    grch38_not_total=$(wc -l < "${GRCH38_NOT_COMPONENT}")

    pct_component_on_grch38=$(awk -v a="${shared_total}" -v b="${component_total}" 'BEGIN{ if(b==0) print "0.000000"; else printf "%.6f", 100*a/b }')
    pct_component_not_on_grch38=$(awk -v a="${component_not_total}" -v b="${component_total}" 'BEGIN{ if(b==0) print "0.000000"; else printf "%.6f", 100*a/b }')
    pct_extra_vs_grch38=$(awk -v a="${component_not_total}" -v b="${grch38_total}" 'BEGIN{ if(b==0) print "0.000000"; else printf "%.6f", 100*a/b }')

    echo -e "${chr}\t${component_total}\t${grch38_total}\t${shared_total}\t${component_not_total}\t${grch38_not_total}\t${pct_component_on_grch38}\t${pct_component_not_on_grch38}\t${pct_extra_vs_grch38}" >> "${SUMMARY}"
done

echo "Done."
echo "Summary: ${SUMMARY}"
