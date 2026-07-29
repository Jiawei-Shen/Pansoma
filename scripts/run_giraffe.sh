#!/usr/bin/env bash
set -euo pipefail

# FASTQ -> GAM with vg giraffe.
#
# Required variables:
#   GBZ, MIN_INDEX, DIST_INDEX, FASTQ1, OUT_GAM
#
# Optional variables:
#   FASTQ2      paired-end FASTQ
#   PLATFORM    illumina, pacbio-hifi, or ont-r10
#   MAPPER_PRESET
#               illumina: default, chaining-sr, fast, or srold
#               pacbio-hifi: hifi
#               ont-r10: r10
#   THREADS     default: 12

: "${GBZ:?Set GBZ to the graph .gbz path}"
: "${MIN_INDEX:?Set MIN_INDEX to the .min index path}"
: "${DIST_INDEX:?Set DIST_INDEX to the .dist index path}"
: "${FASTQ1:?Set FASTQ1 to read 1 FASTQ path}"
: "${OUT_GAM:?Set OUT_GAM to output .gam path}"

THREADS="${THREADS:-12}"
PLATFORM="${PLATFORM:-}"
MAPPER_PRESET="${MAPPER_PRESET:-}"

# Backward compatibility for the unambiguous legacy values only.
if [[ -z "${PLATFORM}" && -n "${READ_TYPE:-}" ]]; then
  case "${READ_TYPE}" in
    illumina) PLATFORM="illumina" ;;
    hifi) PLATFORM="pacbio-hifi" ;;
    r10) PLATFORM="ont-r10" ;;
    *)
      echo "ERROR: legacy READ_TYPE supports only illumina, hifi, or r10; use PLATFORM" >&2
      exit 2
      ;;
  esac
fi

if [[ -z "${PLATFORM}" ]]; then
  echo "ERROR: Set PLATFORM to illumina, pacbio-hifi, or ont-r10" >&2
  exit 2
fi

mkdir -p "$(dirname "${OUT_GAM}")"

args=(
  vg giraffe
  -Z "${GBZ}"
  -m "${MIN_INDEX}"
  -d "${DIST_INDEX}"
  -f "${FASTQ1}"
  -t "${THREADS}"
  -o gam
)

case "${PLATFORM}" in
  illumina)
    preset="${MAPPER_PRESET:-default}"
    case "${preset}" in
      default|chaining-sr|fast|srold) ;;
      *)
        echo "ERROR: Illumina MAPPER_PRESET must be default, chaining-sr, fast, or srold" >&2
        exit 2
        ;;
    esac
    ;;
  pacbio-hifi)
    preset="${MAPPER_PRESET:-hifi}"
    if [[ "${preset}" != "hifi" ]]; then
      echo "ERROR: pacbio-hifi requires MAPPER_PRESET=hifi" >&2
      exit 2
    fi
    if [[ -n "${FASTQ2:-}" ]]; then
      echo "ERROR: pacbio-hifi accepts one FASTQ; do not set FASTQ2" >&2
      exit 2
    fi
    ;;
  ont-r10)
    preset="${MAPPER_PRESET:-r10}"
    if [[ "${preset}" != "r10" ]]; then
      echo "ERROR: ont-r10 requires MAPPER_PRESET=r10" >&2
      exit 2
    fi
    if [[ -n "${FASTQ2:-}" ]]; then
      echo "ERROR: ont-r10 accepts one FASTQ; do not set FASTQ2" >&2
      exit 2
    fi
    ;;
  *)
    echo "ERROR: PLATFORM must be illumina, pacbio-hifi, or ont-r10" >&2
    exit 2
    ;;
esac

if [[ -n "${FASTQ2:-}" ]]; then
  args+=( -f "${FASTQ2}" )
fi

args+=( -b "${preset}" )

echo "Pansoma alignment: platform=${PLATFORM} vg_preset=${preset}" >&2

"${args[@]}" > "${OUT_GAM}"
