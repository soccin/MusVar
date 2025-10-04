#!/bin/bash
#
# Download and process blacklist regions for GRCm38
#
# This script:
#   1. Downloads mm10 blacklist regions
#   2. Converts chromosome names (removes 'chr' prefix)
#   3. Subtracts blacklist regions from main bed file
#   4. Filters out mitochondrial regions

set -euo pipefail

readonly BLACKLIST_URL="https://github.com/Boyle-Lab/Blacklist/raw/refs/heads/master/lists/mm10-blacklist.v2.bed.gz"
readonly MM10_BLACKLIST="mm10-blacklist.v2.bed.gz"
readonly GRCM38_BLACKLIST="GRCm38-blacklist.v2.bed"
readonly MAIN_BED="GRCm38.main.bed"
readonly OUTPUT_BED="GRCm38.subBlack.bed"

main() {
  module load bedtools

  # Download blacklist regions
  wget $BLACKLIST_URL

  # Convert chromosome names from mm10 (chr1) to GRCm38 (1)
  zcat $MM10_BLACKLIST | sed 's/^chr//' > $GRCM38_BLACKLIST

  # Subtract blacklist regions and filter out mitochondrial DNA
  bedtools subtract -a $MAIN_BED -b $GRCM38_BLACKLIST \
    | grep -E -v "^MT" > $OUTPUT_BED
}

main "$@"
