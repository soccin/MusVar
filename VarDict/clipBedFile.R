argv <- commandArgs(trailing = TRUE)

bed_file <- argv[1]
cov_file <- argv[2]

MAX_COVERAGE <- 5000

suppressPackageStartupMessages(require(tidyverse))

# Read coverage file and calculate coverage per region
cov <- read_tsv(cov_file, col_names = FALSE) %>%
  mutate(coverage = .[[ncol(.)]] / (X3 - X2))

# Write normal coverage regions (< MAX_COVERAGE)
cov |>
  filter(coverage < MAX_COVERAGE) |>
  select(X1:X3) |>
  write_tsv(gsub(".bed$", ".clip.bed", bed_file), col_names = FALSE)

# Write abnormal high-coverage regions (>= MAX_COVERAGE) if any exist
high_cov <- cov |>
  filter(coverage >= MAX_COVERAGE)

if (nrow(high_cov) > 0) {
  write_tsv(high_cov, gsub(".bed$", ".BAD", bed_file), col_names = FALSE)
}
