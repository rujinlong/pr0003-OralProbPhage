#!/usr/bin/env Rscript
# Compare our MetaPhlAn4 species table for a bioproject against Geng's
# published Figure 01.Rdata. Pass / fail gate: median |delta| < 0.05.
suppressPackageStartupMessages({
  library(readr); library(dplyr); library(tidyr)
})

args <- commandArgs(trailingOnly = TRUE)
long_tsv <- args[1]
bp       <- args[2]

geng_rdata <- Sys.getenv(
  "GENG_RDATA",
  unset = "../../../analyses/data/00-raw/d020-data-source-registry/wei2024-imeta-repo/Rdata/Figure 01.Rdata"
)
stopifnot(file.exists(geng_rdata))
load(geng_rdata)  # creates `my_data`

ours <- read_tsv(long_tsv, show_col_types = FALSE) |>
  transmute(run_id, species, our_relabund = rel_abund / 100)  # MetaPhlAn rel_ab is %

mat <- my_data$feat_list[[bp]]
stopifnot(!is.null(mat))
geng <- as.data.frame(mat) |>
  tibble::rownames_to_column("run_id") |>
  pivot_longer(-run_id, names_to = "species", values_to = "geng_relabund") |>
  filter(geng_relabund > 0)

joined <- full_join(geng, ours, by = c("run_id","species")) |>
  mutate(across(c(geng_relabund, our_relabund), ~tidyr::replace_na(.x, 0))) |>
  mutate(delta = abs(geng_relabund - our_relabund))

med <- median(joined$delta)
n_geng_only <- sum(joined$geng_relabund > 0 & joined$our_relabund == 0)
frac_missing <- n_geng_only / sum(joined$geng_relabund > 0)

cat(sprintf("[validate] %s: median |delta| = %.4f (gate: <0.05)\n", bp, med))
cat(sprintf("[validate] %s: %d Geng species missing in ours (%.1f%%)\n",
            bp, n_geng_only, 100 * frac_missing))

out <- sub("_long\\.tsv$", "_validation.tsv", long_tsv)
write_tsv(joined, out)
cat(sprintf("[validate] wrote per-(run,species) deltas: %s\n", out))

if (med >= 0.05 || frac_missing > 0.10) {
  cat("[validate] FAIL — investigate DB version / filtering thresholds\n")
  quit(status = 1)
} else {
  cat("[validate] PASS\n")
}
