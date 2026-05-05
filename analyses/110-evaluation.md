# 110-evaluation
Your Name
2026-05-05

- [<span class="toc-section-number">1</span> Pull contract
  tables](#pull-contract-tables)
- [<span class="toc-section-number">2</span> Summary kable: row counts
  per contract table](#summary-kable-row-counts-per-contract-table)
- [<span class="toc-section-number">3</span> Figure 1: Sample × source
  breakdown](#figure-1-sample--source-breakdown)
- [<span class="toc-section-number">4</span> Figure 2: Candidate-microbe
  volcano](#figure-2-candidate-microbe-volcano)
- [<span class="toc-section-number">5</span> Figure 3: Phage-host
  coverage of candidates](#figure-3-phage-host-coverage-of-candidates)
- [<span class="toc-section-number">6</span> Health
  checks](#health-checks)
- [<span class="toc-section-number">7</span> 110-r1: Bootstrap
  stability + leave-one-cohort-out
  (LOOCV)](#110-r1-bootstrap-stability--leave-one-cohort-out-loocv)
  - [<span class="toc-section-number">7.1</span> Section 1: bootstrap
    stability
    (within-cohort)](#section-1-bootstrap-stability-within-cohort)
  - [<span class="toc-section-number">7.2</span> Section 2:
    leave-one-cohort-out (LOOCV) on the
    meta-analysis](#section-2-leave-one-cohort-out-loocv-on-the-meta-analysis)
  - [<span class="toc-section-number">7.3</span> Section 3: 110-r1
    summary](#section-3-110-r1-summary)
- [<span class="toc-section-number">8</span> Files
  written](#files-written)

**Updated: 2026-05-05 22:43:07 CET.**

End-to-end QC of the MVP knowledge layer. Reads the six contract tables
and emits three diagnostic figures plus a summary kable. **No DB
writes** — this module is an evaluator, not a producer.

**Ripple ladder**:

- `110-r1`: bootstrap stability — resample samples and re-derive
  `candidate_microbe`; report fraction of candidates surviving 100
  boots.
- `110-r2`: leave-one-cohort-out validation across all Geng sub-cohorts
  (after 020-r1 + 040-r1 expand the sample set).
- `110-r3`: literature triangulation — for each top candidate, hit
  PubMed via the bioRxiv MCP server and the project’s Zotero library for
  prior evidence.

<details class="code-fold">
<summary>Code</summary>

``` r
suppressPackageStartupMessages({
  library(here)
  library(conflicted)
  library(tidyverse)
  library(DBI)
  devtools::load_all()
})
conflicted::conflicts_prefer(
  dplyr::filter,
  dplyr::lag,
  dplyr::first,
  dplyr::last,
  dplyr::between,
  .quiet = TRUE
)
```

</details>

## Pull contract tables

<details class="code-fold">
<summary>Code</summary>

``` r
con <- load_db()
df_source    <- read_table_db(con, "data_source")
df_sample    <- read_table_db(con, "sample")
df_profile   <- read_table_db(con, "taxon_profile")
df_candidate <- read_table_db(con, "candidate_microbe")
df_link      <- read_table_db(con, "phage_host_link")
df_packet    <- read_table_db(con, "evidence_packet")
close_db(con)
```

</details>

## Summary kable: row counts per contract table

<details class="code-fold">
<summary>Code</summary>

``` r
df_counts <- tibble::tibble(
  table = c("data_source","sample","taxon_profile",
            "candidate_microbe","phage_host_link","evidence_packet"),
  rows  = c(nrow(df_source), nrow(df_sample), nrow(df_profile),
            nrow(df_candidate), nrow(df_link), nrow(df_packet))
)
knitr::kable(df_counts, caption = "Contract table populations after MVP")
```

</details>

| table             |  rows |
|:------------------|------:|
| data_source       |     7 |
| sample            |   227 |
| taxon_profile     | 30063 |
| candidate_microbe |   764 |
| phage_host_link   |   770 |
| evidence_packet   |   798 |

Contract table populations after MVP

## Figure 1: Sample × source breakdown

<details class="code-fold">
<summary>Code</summary>

``` r
df_sample_summary <- df_sample |>
  dplyr::count(source_id, disease_status, body_site, seq_type)
p1 <- ggplot(df_sample_summary,
             aes(x = source_id, y = n, fill = disease_status)) +
  geom_col(position = "stack") +
  geom_text(aes(label = paste0(disease_status, ": ", n)),
            position = position_stack(vjust = 0.5), size = 3) +
  labs(x = NULL, y = "samples",
       title = "Samples per source × disease status",
       subtitle = sprintf("Body sites: %s | seq_types: %s",
                          paste(unique(df_sample$body_site), collapse = ", "),
                          paste(unique(df_sample$seq_type), collapse = ", "))) +
  theme_minimal(base_size = 11) +
  theme(legend.position = "none")
ggsave(path_target("fig1_samples.png"), p1,
       width = 7, height = 3.5, dpi = 150)
print(p1)
```

</details>

<div id="fig-1">

<img src="110-evaluation_files/figure-commonmark/fig-1-1.png"
id="fig-1" />

Figure 1

</div>

## Figure 2: Candidate-microbe volcano

<details class="code-fold">
<summary>Code</summary>

``` r
df_volcano <- df_candidate |>
  dplyr::mutate(
    nlog10_p = -log10(p_value),
    label = ifelse(rank(p_value) <= 10, sub("^STR_", "", taxon_id), NA_character_)
  )
p2 <- ggplot(df_volcano,
             aes(x = effect_size, y = nlog10_p, color = direction)) +
  geom_point(alpha = 0.85, size = 2.4) +
  geom_text(aes(label = label), hjust = -0.15, size = 3, na.rm = TRUE) +
  scale_color_manual(values = c(disease_enriched = "#c84a4a",
                                health_enriched = "#3f86b8")) +
  labs(x = "log2 fold change (peri / healthy)",
       y = expression(-log[10](p)),
       title = "M1a candidate microbes — Wilcoxon volcano (PRJDB11203)") +
  theme_minimal(base_size = 11)
ggsave(path_target("fig2_volcano.png"), p2,
       width = 7, height = 5, dpi = 150)
print(p2)
```

</details>

<div id="fig-2">

<img src="110-evaluation_files/figure-commonmark/fig-2-1.png"
id="fig-2" />

Figure 2

</div>

## Figure 3: Phage-host coverage of candidates

<details class="code-fold">
<summary>Code</summary>

``` r
df_cov <- df_candidate |>
  dplyr::left_join(
    df_packet |>
      dplyr::group_by(candidate_id) |>
      dplyr::summarise(n_phages = dplyr::n_distinct(link_id),
                       .groups = "drop"),
    by = "candidate_id"
  ) |>
  dplyr::mutate(n_phages = tidyr::replace_na(n_phages, 0L),
                covered  = n_phages > 0)

cov_pct <- mean(df_cov$covered) * 100
p3 <- ggplot(df_cov, aes(x = forcats::fct_reorder(taxon_name, n_phages),
                         y = n_phages, fill = direction)) +
  geom_col() +
  coord_flip() +
  scale_fill_manual(values = c(disease_enriched = "#c84a4a",
                               health_enriched = "#3f86b8")) +
  labs(x = NULL, y = "phages predicted to infect this species' genus",
       title = "Phage-host coverage of M1a candidates",
       subtitle = sprintf("%.0f%% of candidates have ≥1 predicted phage",
                          cov_pct)) +
  theme_minimal(base_size = 10)
ggsave(path_target("fig3_phage_coverage.png"), p3,
       width = 7, height = 4, dpi = 150)
print(p3)
```

</details>

<div id="fig-3">

<img src="110-evaluation_files/figure-commonmark/fig-3-1.png"
id="fig-3" />

Figure 3

</div>

## Health checks

<details class="code-fold">
<summary>Code</summary>

``` r
checks <- list(
  six_tables_nonempty = all(df_counts$rows > 0),
  evidence_packets_have_rows = nrow(df_packet) >= 1,
  three_figures_written = length(
    list.files(path_target(), pattern = "^fig.*\\.png$")
  ) >= 3,
  no_orphan_packets = all(df_packet$candidate_id %in% df_candidate$candidate_id) &&
                      all(df_packet$link_id %in% df_link$link_id)
)
df_checks <- tibble::tibble(check = names(checks),
                            pass  = unlist(checks))
knitr::kable(df_checks, caption = "MVP DoD health checks")
```

</details>

| check                      | pass |
|:---------------------------|:-----|
| six_tables_nonempty        | TRUE |
| evidence_packets_have_rows | TRUE |
| three_figures_written      | TRUE |
| no_orphan_packets          | TRUE |

MVP DoD health checks

<details class="code-fold">
<summary>Code</summary>

``` r
stopifnot(all(df_checks$pass))
```

</details>

## 110-r1: Bootstrap stability + leave-one-cohort-out (LOOCV)

This section quantifies the robustness of `050-r1` (per-cohort Wilcoxon)
and `050-r2` (Stouffer’s Z meta-analysis) candidates. Two complementary
checks:

- **Within-cohort bootstrap** (`B = 200`, stratified by
  `disease_status`): for each Geng 2024 cohort, resample with
  replacement preserving group sizes, re-run two-sided Wilcoxon,
  BH-correct within bootstrap, and record the fraction of bootstraps in
  which each species crosses `q < 0.05`.
- **Leave-one-cohort-out (LOOCV)**: drop each of the 6 Geng 2024 cohorts
  in turn and re-run the meta-analysis (Stouffer’s Z, sqrt-n weighted)
  on the remaining 5. A meta_stouffer candidate is “LOOCV-robust” if it
  remains significant (q \< 0.05) in `>= ceiling(0.83 * eligible_folds)`
  LOO folds where the species was present in `>= 2` of the 5 retained
  cohorts.

<details class="code-fold">
<summary>Code</summary>

``` r
con <- load_db()
df_sample_geng <- read_table_db(con, "sample") |>
  dplyr::filter(grepl("^geng2024_", source_id),
                disease_status %in% c("periodontitis", "healthy"))
df_profile_geng <- read_table_db(con, "taxon_profile") |>
  dplyr::filter(sample_id %in% df_sample_geng$sample_id,
                taxon_kind == "bacterium")
df_meta_stouffer_full <- read_table_db(con, "candidate_microbe") |>
  dplyr::filter(method == "meta_stouffer", phenotype == "periodontitis")
close_db(con)

cohorts_geng <- sort(unique(df_sample_geng$source_id))
cat("Cohorts:", length(cohorts_geng),
    "| samples:", nrow(df_sample_geng),
    "| meta_stouffer rows in DB:", nrow(df_meta_stouffer_full), "\n")
```

</details>

    Cohorts: 6 | samples: 223 | meta_stouffer rows in DB: 334 

### Section 1: bootstrap stability (within-cohort)

<details class="code-fold">
<summary>Code</summary>

``` r
B_BOOT     <- 200L
SEED_BOOT  <- 20260505L
PREV_CUT   <- 0.10
Q_CUT      <- 0.05

build_cohort_mat <- function(source_id_x) {
  df_s <- df_sample_geng |> dplyr::filter(source_id == source_id_x)
  df_p <- df_profile_geng |> dplyr::filter(sample_id %in% df_s$sample_id)
  df_w <- df_p |>
    dplyr::select(sample_id, taxon_id, taxon_name, abundance) |>
    tidyr::pivot_wider(id_cols = c(taxon_id, taxon_name),
                       names_from = sample_id, values_from = abundance,
                       values_fill = 0)
  ids_case <- df_s$sample_id[df_s$disease_status == "periodontitis"]
  ids_ctrl <- df_s$sample_id[df_s$disease_status == "healthy"]
  mat <- as.matrix(df_w[, c(ids_case, ids_ctrl)])
  rownames(mat) <- df_w$taxon_id
  prev <- rowMeans(mat > 0)
  keep <- prev >= PREV_CUT
  list(
    source_id  = source_id_x,
    mat        = mat[keep, , drop = FALSE],
    taxon_name = setNames(df_w$taxon_name, df_w$taxon_id)[rownames(mat)[keep]],
    ids_case   = ids_case,
    ids_ctrl   = ids_ctrl
  )
}

cohort_data <- purrr::map(cohorts_geng, build_cohort_mat) |>
  setNames(cohorts_geng)

bootstrap_one_cohort <- function(cd, B, seed) {
  set.seed(seed)
  taxa <- rownames(cd$mat)
  n1 <- length(cd$ids_case); n2 <- length(cd$ids_ctrl)
  hits <- integer(length(taxa)); names(hits) <- taxa
  if (n1 < 5 || n2 < 5 || length(taxa) == 0) {
    return(tibble::tibble(taxon_id = taxa,
                          n_bootstraps_significant = hits,
                          stability_score = 0))
  }
  for (b in seq_len(B)) {
    idx_case <- sample(cd$ids_case, n1, replace = TRUE)
    idx_ctrl <- sample(cd$ids_ctrl, n2, replace = TRUE)
    sub_case <- cd$mat[, idx_case, drop = FALSE]
    sub_ctrl <- cd$mat[, idx_ctrl, drop = FALSE]
    pv <- vapply(seq_along(taxa), function(i) {
      a <- sub_case[i, ]; b2 <- sub_ctrl[i, ]
      if (length(unique(c(a, b2))) < 2) return(NA_real_)
      suppressWarnings(stats::wilcox.test(a, b2, exact = FALSE)$p.value)
    }, numeric(1))
    qv <- stats::p.adjust(pv, method = "BH")
    hits <- hits + as.integer(!is.na(qv) & qv < Q_CUT)
  }
  tibble::tibble(
    taxon_id                 = taxa,
    n_bootstraps_significant = hits,
    stability_score          = hits / B
  )
}

t0 <- Sys.time()
df_boot <- purrr::imap_dfr(cohort_data, function(cd, src) {
  bootstrap_one_cohort(cd, B_BOOT, SEED_BOOT + which(cohorts_geng == src)) |>
    dplyr::mutate(cohort = src,
                  taxon_name = unname(cd$taxon_name[taxon_id]))
})
t1 <- Sys.time()
cat(sprintf("Bootstrap done: %d cohorts x %d boots = %d total bootstraps in %.1f min\n",
            length(cohort_data), B_BOOT, length(cohort_data) * B_BOOT,
            as.numeric(difftime(t1, t0, units = "mins"))))
```

</details>

    Bootstrap done: 6 cohorts x 200 boots = 1200 total bootstraps in 0.4 min

<details class="code-fold">
<summary>Code</summary>

``` r
df_boot_out <- df_boot |>
  dplyr::select(cohort, taxon_id, taxon_name,
                n_bootstraps_significant, stability_score) |>
  dplyr::arrange(cohort, dplyr::desc(stability_score))

write.table(df_boot_out, path_target("bootstrap_stability.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)
cat("Wrote bootstrap_stability.tsv: ", nrow(df_boot_out), " rows\n")
```

</details>

    Wrote bootstrap_stability.tsv:  1535  rows

<details class="code-fold">
<summary>Code</summary>

``` r
p_boot <- ggplot(df_boot_out, aes(x = stability_score)) +
  geom_histogram(binwidth = 0.05, boundary = 0,
                 fill = "#3f86b8", color = "white") +
  facet_wrap(~ cohort, ncol = 3, scales = "free_y") +
  labs(x = "stability score (fraction of B=200 bootstraps with q < 0.05)",
       y = "species",
       title = "Bootstrap stability per cohort",
       subtitle = sprintf("B = %d, stratified by disease_status, BH within bootstrap",
                          B_BOOT)) +
  theme_minimal(base_size = 10)
ggsave(path_target("fig4_bootstrap_stability_hist.png"), p_boot,
       width = 9, height = 5.5, dpi = 150)
print(p_boot)
```

</details>

![](110-evaluation_files/figure-commonmark/r1-bootstrap-plot-1.png)

<details class="code-fold">
<summary>Code</summary>

``` r
df_boot_cohort <- df_boot_out |>
  dplyr::group_by(cohort) |>
  dplyr::summarise(
    n_species          = dplyr::n(),
    median_stability   = stats::median(stability_score),
    mean_stability     = mean(stability_score),
    n_high             = sum(stability_score >= 0.8),
    .groups = "drop"
  )
knitr::kable(df_boot_cohort, digits = 3,
             caption = "Per-cohort bootstrap stability summary")
```

</details>

| cohort               | n_species | median_stability | mean_stability | n_high |
|:---------------------|----------:|-----------------:|---------------:|-------:|
| geng2024_PRJDB11203  |       136 |            0.035 |          0.112 |      0 |
| geng2024_PRJNA230363 |       280 |            0.180 |          0.272 |     18 |
| geng2024_PRJNA396840 |       210 |            0.015 |          0.056 |      0 |
| geng2024_PRJNA678453 |       293 |            0.150 |          0.263 |     27 |
| geng2024_PRJNA717815 |       283 |            0.210 |          0.352 |     50 |
| geng2024_PRJNA932553 |       333 |            0.085 |          0.146 |      2 |

Per-cohort bootstrap stability summary

### Section 2: leave-one-cohort-out (LOOCV) on the meta-analysis

<details class="code-fold">
<summary>Code</summary>

``` r
per_cohort_wilcox_loo <- function(source_id_x) {
  cd <- cohort_data[[source_id_x]]
  if (is.null(cd)) return(NULL)
  if (length(cd$ids_case) < 5 || length(cd$ids_ctrl) < 5) return(NULL)
  taxa <- rownames(cd$mat)
  if (length(taxa) == 0) return(NULL)
  ids_case <- cd$ids_case; ids_ctrl <- cd$ids_ctrl
  n1 <- length(ids_case); n2 <- length(ids_ctrl); n_total <- n1 + n2
  mat_case <- cd$mat[, ids_case, drop = FALSE]
  mat_ctrl <- cd$mat[, ids_ctrl, drop = FALSE]
  out <- purrr::map_dfr(seq_along(taxa), function(i) {
    a <- mat_case[i, ]; b <- mat_ctrl[i, ]
    if (length(unique(c(a, b))) < 2) {
      return(tibble::tibble(taxon_id = taxa[i],
                            p_value = NA_real_,
                            log2fc = 0, n_total = n_total))
    }
    wt <- suppressWarnings(stats::wilcox.test(a, b, exact = FALSE))
    tibble::tibble(taxon_id = taxa[i],
                   p_value  = wt$p.value,
                   log2fc   = log2((mean(a) + 1e-6) / (mean(b) + 1e-6)),
                   n_total  = n_total)
  })
  out |>
    dplyr::mutate(source_id  = source_id_x,
                  taxon_name = unname(cd$taxon_name[taxon_id])) |>
    dplyr::filter(!is.na(p_value), is.finite(p_value),
                  p_value > 0, p_value < 1)
}

df_per_cohort_stats <- purrr::map_dfr(cohorts_geng, per_cohort_wilcox_loo)
cat("Per-cohort Wilcoxon table:", nrow(df_per_cohort_stats),
    "rows across", dplyr::n_distinct(df_per_cohort_stats$source_id), "cohorts\n")
```

</details>

    Per-cohort Wilcoxon table: 1513 rows across 6 cohorts

<details class="code-fold">
<summary>Code</summary>

``` r
stouffer_pool <- function(df_one) {
  pv <- df_one$p_value
  fc <- df_one$log2fc
  n  <- df_one$n_total
  z_two <- stats::qnorm(1 - pv / 2)
  z_two[!is.finite(z_two)] <- 0
  z_signed <- z_two * ifelse(fc >= 0, 1, -1)
  w <- sqrt(n)
  z_pool <- sum(w * z_signed) / sqrt(sum(w^2))
  p_pool <- 2 * stats::pnorm(-abs(z_pool))
  tibble::tibble(z_pool = z_pool, p_value = p_pool,
                 effect_size = mean(fc),
                 direction   = ifelse(z_pool >= 0, "disease_enriched",
                                                   "health_enriched"))
}
```

</details>

<details class="code-fold">
<summary>Code</summary>

``` r
run_loo_meta <- function(dropped_cohort) {
  df_kept <- df_per_cohort_stats |>
    dplyr::filter(source_id != dropped_cohort)
  df_grp <- df_kept |>
    dplyr::group_by(taxon_id) |>
    dplyr::summarise(n_cohorts = dplyr::n_distinct(source_id),
                     taxon_name = dplyr::first(stats::na.omit(taxon_name)),
                     .groups = "drop") |>
    dplyr::filter(n_cohorts >= 2)
  if (nrow(df_grp) == 0) return(tibble::tibble())
  df_meta <- df_kept |>
    dplyr::semi_join(df_grp, by = "taxon_id") |>
    dplyr::group_by(taxon_id) |>
    dplyr::group_modify(~ stouffer_pool(.x)) |>
    dplyr::ungroup() |>
    dplyr::left_join(df_grp, by = "taxon_id") |>
    dplyr::mutate(q_value = stats::p.adjust(p_value, method = "BH"),
                  cohort_dropped = dropped_cohort)
  df_meta |>
    dplyr::select(cohort_dropped, taxon_id, taxon_name,
                  n_cohorts, direction, effect_size,
                  z_pool, p_value, q_value)
}

df_loo <- purrr::map_dfr(cohorts_geng, run_loo_meta)
cat("LOOCV: 6 folds x species testable (>=2 of 5 cohorts) = ",
    nrow(df_loo), " rows\n", sep = "")
```

</details>

    LOOCV: 6 folds x species testable (>=2 of 5 cohorts) = 1920 rows

<details class="code-fold">
<summary>Code</summary>

``` r
df_loo_out <- df_loo |>
  dplyr::select(cohort_dropped, taxon_id, taxon_name,
                n_cohorts, direction, effect_size, p_value, q_value)
write.table(df_loo_out, path_target("loocv.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)
cat("Wrote loocv.tsv\n")
```

</details>

    Wrote loocv.tsv

<details class="code-fold">
<summary>Code</summary>

``` r
df_meta_full_q05 <- df_meta_stouffer_full |>
  dplyr::filter(q_value < 0.05)

df_loo_robust <- df_loo |>
  dplyr::filter(taxon_id %in% df_meta_full_q05$taxon_id) |>
  dplyr::group_by(taxon_id) |>
  dplyr::summarise(
    taxon_name      = dplyr::first(stats::na.omit(taxon_name)),
    eligible_folds  = dplyr::n(),
    folds_q05       = sum(q_value < 0.05, na.rm = TRUE),
    folds_q01       = sum(q_value < 0.01, na.rm = TRUE),
    direction_consistent =
      length(unique(stats::na.omit(direction))) == 1,
    .groups = "drop"
  ) |>
  dplyr::mutate(
    cutoff_n_robust = pmax(1L, as.integer(ceiling(0.83 * eligible_folds))),
    is_loo_robust   = folds_q05 >= cutoff_n_robust & direction_consistent
  )

n_loo_robust <- sum(df_loo_robust$is_loo_robust)
cat("LOOCV-robust species (among meta_stouffer q<0.05, n=",
    nrow(df_meta_full_q05), "): ", n_loo_robust, "\n", sep = "")
```

</details>

    LOOCV-robust species (among meta_stouffer q<0.05, n=111): 83

<details class="code-fold">
<summary>Code</summary>

``` r
top30_taxa <- df_meta_stouffer_full |>
  dplyr::arrange(q_value) |>
  head(30) |>
  dplyr::pull(taxon_id)

df_heat <- df_loo |>
  dplyr::filter(taxon_id %in% top30_taxa) |>
  dplyr::mutate(
    nlogq      = -log10(pmax(q_value, 1e-300)),
    cohort_lab = sub("^geng2024_", "", cohort_dropped),
    sp_lab     = sub("^STR_", "", taxon_id)
  )

taxa_ordered <- df_meta_stouffer_full |>
  dplyr::filter(taxon_id %in% top30_taxa) |>
  dplyr::arrange(dplyr::desc(q_value)) |>
  dplyr::mutate(sp_lab = sub("^STR_", "", taxon_id)) |>
  dplyr::pull(sp_lab)

df_heat$sp_lab <- factor(df_heat$sp_lab, levels = taxa_ordered)

p_heat <- ggplot(df_heat,
                 aes(x = cohort_lab, y = sp_lab, fill = nlogq)) +
  geom_tile(color = "white", linewidth = 0.3) +
  geom_text(aes(label = ifelse(q_value < 0.001, "***",
                       ifelse(q_value < 0.01,  "**",
                       ifelse(q_value < 0.05,  "*", "")))),
            size = 3, color = "white") +
  scale_fill_gradient(low = "#f1eef6", high = "#67000d",
                      name = expression(-log[10](q))) +
  labs(x = "cohort dropped", y = NULL,
       title = "LOOCV: top 30 meta_stouffer candidates across 6 LOO folds",
       subtitle = "* q<0.05, ** q<0.01, *** q<0.001 in the 5-cohort meta after dropping x") +
  theme_minimal(base_size = 9) +
  theme(axis.text.x = element_text(angle = 30, hjust = 1),
        panel.grid = element_blank())
ggsave(path_target("fig5_loocv_top30_heatmap.png"), p_heat,
       width = 8.5, height = 9, dpi = 150)
print(p_heat)
```

</details>

![](110-evaluation_files/figure-commonmark/r1-loocv-heatmap-1.png)

### Section 3: 110-r1 summary

<details class="code-fold">
<summary>Code</summary>

``` r
n_high_stable_any <- df_boot_out |>
  dplyr::filter(stability_score >= 0.8) |>
  dplyr::pull(taxon_id) |>
  unique() |>
  length()

n_meta_q05 <- nrow(df_meta_full_q05)
cross_table <- tibble::tibble(
  bucket = c("meta_stouffer q<0.05",
             "of which LOOCV-robust (>=ceil(0.83*eligible))",
             "of which NOT LOOCV-robust"),
  n      = c(n_meta_q05, n_loo_robust, n_meta_q05 - n_loo_robust)
)

df_r1_summary <- tibble::tibble(
  metric = c(
    "B (bootstraps per cohort)",
    "n cohorts",
    "n bootstraps total (B x cohorts)",
    "n highly-stable species (>=0.8 in any cohort)",
    "n meta_stouffer candidates (q<0.05, full 6-cohort)",
    "n LOOCV-robust species among those candidates"
  ),
  value = c(
    as.character(B_BOOT),
    as.character(length(cohorts_geng)),
    as.character(B_BOOT * length(cohorts_geng)),
    as.character(n_high_stable_any),
    as.character(n_meta_q05),
    as.character(n_loo_robust)
  )
)
knitr::kable(df_r1_summary, caption = "110-r1 headline numbers")
```

</details>

| metric                                              | value |
|:----------------------------------------------------|:------|
| B (bootstraps per cohort)                           | 200   |
| n cohorts                                           | 6     |
| n bootstraps total (B x cohorts)                    | 1200  |
| n highly-stable species (\>=0.8 in any cohort)      | 74    |
| n meta_stouffer candidates (q\<0.05, full 6-cohort) | 111   |
| n LOOCV-robust species among those candidates       | 83    |

110-r1 headline numbers

<details class="code-fold">
<summary>Code</summary>

``` r
knitr::kable(df_boot_cohort, digits = 3,
             caption = "Per-cohort bootstrap (median stability + n>=0.8)")
```

</details>

| cohort               | n_species | median_stability | mean_stability | n_high |
|:---------------------|----------:|-----------------:|---------------:|-------:|
| geng2024_PRJDB11203  |       136 |            0.035 |          0.112 |      0 |
| geng2024_PRJNA230363 |       280 |            0.180 |          0.272 |     18 |
| geng2024_PRJNA396840 |       210 |            0.015 |          0.056 |      0 |
| geng2024_PRJNA678453 |       293 |            0.150 |          0.263 |     27 |
| geng2024_PRJNA717815 |       283 |            0.210 |          0.352 |     50 |
| geng2024_PRJNA932553 |       333 |            0.085 |          0.146 |      2 |

Per-cohort bootstrap (median stability + n\>=0.8)

<details class="code-fold">
<summary>Code</summary>

``` r
knitr::kable(cross_table,
             caption = "Cross-table: meta_stouffer candidates x LOOCV robustness")
```

</details>

| bucket                                          |   n |
|:------------------------------------------------|----:|
| meta_stouffer q\<0.05                           | 111 |
| of which LOOCV-robust (\>=ceil(0.83\*eligible)) |  83 |
| of which NOT LOOCV-robust                       |  28 |

Cross-table: meta_stouffer candidates x LOOCV robustness

<details class="code-fold">
<summary>Code</summary>

``` r
known <- c("STR_Filifactor_alocis",
           "STR_Porphyromonas_gingivalis",
           "STR_Eubacterium_nodatum",
           "STR_Treponema_denticola",
           "STR_Prevotella_intermedia")

df_known_boot <- df_boot_out |>
  dplyr::filter(taxon_id %in% known) |>
  dplyr::group_by(taxon_id, taxon_name) |>
  dplyr::summarise(
    n_cohorts_high  = sum(stability_score >= 0.8),
    median_score    = stats::median(stability_score),
    .groups = "drop"
  ) |>
  dplyr::arrange(dplyr::desc(median_score))
knitr::kable(df_known_boot, digits = 3,
             caption = "Sanity check: known periodontitis species in bootstrap")
```

</details>

| taxon_id | taxon_name | n_cohorts_high | median_score |
|:---|:---|---:|---:|
| STR_Eubacterium_nodatum | Eubacterium_nodatum | 2 | 0.720 |
| STR_Treponema_denticola | Treponema_denticola | 3 | 0.712 |
| STR_Porphyromonas_gingivalis | Porphyromonas_gingivalis | 2 | 0.705 |
| STR_Filifactor_alocis | Filifactor_alocis | 3 | 0.675 |
| STR_Prevotella_intermedia | Prevotella_intermedia | 2 | 0.515 |

Sanity check: known periodontitis species in bootstrap

<details class="code-fold">
<summary>Code</summary>

``` r
df_known_loo <- df_loo_robust |>
  dplyr::filter(taxon_id %in% known) |>
  dplyr::select(taxon_id, taxon_name, eligible_folds, folds_q05, is_loo_robust)
knitr::kable(df_known_loo,
             caption = "Sanity check: same species in LOOCV")
```

</details>

| taxon_id | taxon_name | eligible_folds | folds_q05 | is_loo_robust |
|:---|:---|---:|---:|:---|
| STR_Eubacterium_nodatum | Eubacterium_nodatum | 6 | 6 | TRUE |
| STR_Filifactor_alocis | Filifactor_alocis | 6 | 6 | TRUE |
| STR_Porphyromonas_gingivalis | Porphyromonas_gingivalis | 6 | 6 | TRUE |
| STR_Prevotella_intermedia | Prevotella_intermedia | 6 | 6 | TRUE |
| STR_Treponema_denticola | Treponema_denticola | 6 | 6 | TRUE |

Sanity check: same species in LOOCV

## Files written

<details class="code-fold">
<summary>Code</summary>

``` r
knitr::kable(qproj::proj_dir_info(path_target(), tz = "CET"))
```

</details>

| path                              | type |   size | modification_time   |
|:----------------------------------|:-----|-------:|:--------------------|
| bootstrap_stability.tsv           | file | 119.4K | 2026-05-05 22:43:31 |
| fig1_samples.png                  | file |  44.5K | 2026-05-05 22:43:08 |
| fig2_volcano.png                  | file |  90.5K | 2026-05-05 22:43:08 |
| fig3_phage_coverage.png           | file | 223.1K | 2026-05-05 22:43:09 |
| fig4_bootstrap_stability_hist.png | file |  57.5K | 2026-05-05 22:43:31 |
| fig5_loocv_top30_heatmap.png      | file | 130.1K | 2026-05-05 22:43:33 |
| loocv.tsv                         | file | 272.3K | 2026-05-05 22:43:33 |
