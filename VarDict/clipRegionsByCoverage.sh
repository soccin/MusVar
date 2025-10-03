#!/bin/bash
#
# Clip BED regions based on coverage thresholds.
# Computes coverage for each BED region and filters out high-coverage regions.

set -euo pipefail

readonly SCRIPT_DIR="$(dirname "$(readlink -f "$0")")"

usage() {
  cat << EOF
Usage: $(basename "$0") BAM BED

Compute per-region coverage from BAM file and filter BED regions.
Regions with coverage >= MAX_COVERAGE threshold are written to .BAD file.

Arguments:
  BAM    Input BAM file
  BED    Input BED file (regions to analyze)

Output:
  BED.clip.bed    Regions with normal coverage
  BED.BAD         Regions with excessive coverage (if any)
  cov/BED.cov     Raw coverage data
EOF
  exit 1
}

main() {
  # Check arguments
  if [[ $# -ne 2 ]]; then
    usage
  fi

  local bam_file="$1"
  local bed_file="$2"

  # Validate input files
  if [[ ! -f "${bam_file}" ]]; then
    echo "ERROR: BAM file not found: ${bam_file}" >&2
    exit 1
  fi

  if [[ ! -f "${bed_file}" ]]; then
    echo "ERROR: BED file not found: ${bed_file}" >&2
    exit 1
  fi

  module load samtools

  # Store coverage data in subdirectory to avoid clutter
  local cov_dir="$(dirname "${bed_file}")/cov"
  mkdir -p "${cov_dir}"

  local cov_file="${cov_dir}/$(basename "${bed_file/.bed/.cov}")"

  # Compute per-region coverage
  samtools bedcov "${bed_file}" "${bam_file}" > "${cov_file}"

  # Filter regions by coverage threshold
  Rscript "${SCRIPT_DIR}/clipBedFile.R" "${bed_file}" "${cov_file}"
}

main "$@"
