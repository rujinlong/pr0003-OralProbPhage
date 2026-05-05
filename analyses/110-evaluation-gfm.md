Your Name

- [<span class="toc-section-number">1</span>
  110-evaluation](#110-evaluation)
  - [<span class="toc-section-number">1.1</span> Pull contract
    tables](#pull-contract-tables)
  - [<span class="toc-section-number">1.2</span> Summary kable: row
    counts per contract
    table](#summary-kable-row-counts-per-contract-table)
  - [<span class="toc-section-number">1.3</span> Figure 1: Sample ×
    source breakdown](#figure-1-sample--source-breakdown)
  - [<span class="toc-section-number">1.4</span> Figure 2:
    Candidate-microbe volcano](#figure-2-candidate-microbe-volcano)
  - [<span class="toc-section-number">1.5</span> Figure 3: Phage-host
    coverage of candidates](#figure-3-phage-host-coverage-of-candidates)
  - [<span class="toc-section-number">1.6</span> Health
    checks](#health-checks)
  - [<span class="toc-section-number">1.7</span> Files
    written](#files-written)

# 110-evaluation

Your Name 2026-05-05

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
- [<span class="toc-section-number">7</span> Files
  written](#files-written)

**Updated: 2026-05-05 21:13:03 CET.**

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

<summary>

Code
</summary>

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

<summary>

Code
</summary>

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

<summary>

Code
</summary>

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
| candidate_microbe |    96 |
| phage_host_link   |   770 |
| evidence_packet   |   344 |

Contract table populations after MVP

## Figure 1: Sample × source breakdown

<details class="code-fold">

<summary>

Code
</summary>

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

Figure 1: Figure 1

</div>

## Figure 2: Candidate-microbe volcano

<details class="code-fold">

<summary>

Code
</summary>

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

Figure 2: Figure 2

</div>

## Figure 3: Phage-host coverage of candidates

<details class="code-fold">

<summary>

Code
</summary>

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

Figure 3: Figure 3

</div>

## Health checks

<details class="code-fold">

<summary>

Code
</summary>

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

<summary>

Code
</summary>

``` r
stopifnot(all(df_checks$pass))
```

</details>

## Files written

<details class="code-fold">

<summary>

Code
</summary>

``` r
knitr::kable(qproj::proj_dir_info(path_target(), tz = "CET"))
```

</details>

| path                    | type |   size | modification_time   |
|:------------------------|:-----|-------:|:--------------------|
| fig1_samples.png        | file |  44.5K | 2026-05-05 21:13:04 |
| fig2_volcano.png        | file |    59K | 2026-05-05 21:13:04 |
| fig3_phage_coverage.png | file | 109.6K | 2026-05-05 21:13:04 |
