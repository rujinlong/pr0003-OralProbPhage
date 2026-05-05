Your Name

- [<span class="toc-section-number">1</span>
  090-m2a-mapping](#090-m2a-mapping)
  - [<span class="toc-section-number">1.1</span> Pull contract
    tables](#pull-contract-tables)
  - [<span class="toc-section-number">1.2</span> Strategy 1: strict
    equality](#strategy-1-strict-equality)
  - [<span class="toc-section-number">1.3</span> Strategy 2:
    genus-prefix join (MVP
    primary)](#strategy-2-genus-prefix-join-mvp-primary)
  - [<span class="toc-section-number">1.4</span> Per-genus
    summary](#per-genus-summary)
  - [<span class="toc-section-number">1.5</span> Persist mapping
    artefacts](#persist-mapping-artefacts)
  - [<span class="toc-section-number">1.6</span> Quick
    visualization](#quick-visualization)
  - [<span class="toc-section-number">1.7</span> Files
    written](#files-written)

# 090-m2a-mapping

Your Name 2026-05-05

- [<span class="toc-section-number">1</span> Pull contract
  tables](#pull-contract-tables)
- [<span class="toc-section-number">2</span> Strategy 1: strict
  equality](#strategy-1-strict-equality)
- [<span class="toc-section-number">3</span> Strategy 2: genus-prefix
  join (MVP primary)](#strategy-2-genus-prefix-join-mvp-primary)
- [<span class="toc-section-number">4</span> Per-genus
  summary](#per-genus-summary)
- [<span class="toc-section-number">5</span> Persist mapping
  artefacts](#persist-mapping-artefacts)
- [<span class="toc-section-number">6</span> Quick
  visualization](#quick-visualization)
- [<span class="toc-section-number">7</span> Files
  written](#files-written)

**Updated: 2026-05-05 21:13:35 CET.**

Joins `candidate_microbe` (M1a output) with `phage_host_link` (M2a
output) to surface the phages whose predicted hosts overlap with
periodontitis candidates. The MVP uses two join strategies:

1.  **Strict** —
    `candidate_microbe.taxon_id == phage_host_link.host_taxon_id`.
    Likely zero hits in this MVP because Geng’s bacteriome is at species
    level (e.g. `STR_Streptococcus_salivarius`) while Yahara’s CAT
    predictions are at genus level (e.g. `STR_Streptococcus`).
2.  **Genus prefix** — candidate species belongs to the predicted host
    genus (`candidate.taxon_id LIKE host_taxon_id || '_%'`). This is the
    primary mapping for MVP downstream use.

`100-knowledge-integration` consumes the genus-prefix mapping to
assemble `evidence_packet`. This module’s outputs are exploratory
artefacts (CSVs + plots), not contract-layer writes.

**Ripple ladder**:

- `090-r1`: replace genus prefix with eHOMD-aware fuzzy taxonomy (handle
  synonyms, basonyms, misspellings, Candidatus naming).
- `090-r2`: bipartite-graph visualizations (heatmap, network);
  centrality metrics; partitioning of phage–host modules.

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
df_candidate <- read_table_db(con, "candidate_microbe")
df_link      <- read_table_db(con, "phage_host_link")
close_db(con)
cat("candidate_microbe:", nrow(df_candidate),
    " | phage_host_link:", nrow(df_link), "\n")
```

</details>

    candidate_microbe: 96  | phage_host_link: 770 

## Strategy 1: strict equality

<details class="code-fold">

<summary>

Code
</summary>

``` r
df_strict <- df_candidate |>
  dplyr::inner_join(df_link, by = c("taxon_id" = "host_taxon_id")) |>
  dplyr::select(candidate_id, taxon_id, taxon_name, phenotype, direction,
                effect_size, q_value,
                link_id, vOTU_id, evidence_kind, confidence)
cat("Strict-equality hits:", nrow(df_strict), "\n")
```

</details>

    Strict-equality hits: 8 

## Strategy 2: genus-prefix join (MVP primary)

<details class="code-fold">

<summary>

Code
</summary>

``` r
df_prefix <- df_candidate |>
  dplyr::inner_join(df_link,
                    by = character(),  # cross join
                    relationship = "many-to-many") |>
  dplyr::filter(stringr::str_starts(taxon_id,
                                    paste0(host_taxon_id, "_"))) |>
  dplyr::select(candidate_id, taxon_id, taxon_name, phenotype, direction,
                effect_size, q_value,
                link_id, vOTU_id, host_taxon_id, evidence_kind, confidence)
cat("Genus-prefix hits:",
    nrow(df_prefix),
    " | unique candidates with phage hits:",
    dplyr::n_distinct(df_prefix$candidate_id),
    " | unique phages:",
    dplyr::n_distinct(df_prefix$vOTU_id), "\n")
```

</details>

    Genus-prefix hits: 798  | unique candidates with phage hits: 48  | unique phages: 480 

## Per-genus summary

<details class="code-fold">

<summary>

Code
</summary>

``` r
df_genus_summary <- df_prefix |>
  dplyr::group_by(host_taxon_id) |>
  dplyr::summarise(
    n_candidates = dplyr::n_distinct(candidate_id),
    n_phages     = dplyr::n_distinct(vOTU_id),
    direction_breakdown = paste(
      paste0(unique(direction), "(",
             vapply(unique(direction), function(d) dplyr::n_distinct(candidate_id[direction == d]), integer(1)),
             ")"),
      collapse = ", "
    ),
    mean_confidence = mean(confidence, na.rm = TRUE),
    .groups = "drop"
  ) |>
  dplyr::arrange(dplyr::desc(n_candidates), dplyr::desc(n_phages))
knitr::kable(df_genus_summary, caption = "Phage–host genus mapping summary")
```

</details>

| host_taxon_id | n_candidates | n_phages | direction_breakdown | mean_confidence |
|:---|---:|---:|:---|---:|
| STR_Actinomyces | 9 | 6 | health_enriched(9) | 0.7216667 |
| STR_Treponema | 9 | 1 | disease_enriched(9) | 0.9300000 |
| STR_Corynebacterium | 5 | 2 | health_enriched(5) | 0.7400000 |
| STR_Porphyromonas | 5 | 2 | disease_enriched(5) | 0.5900000 |
| STR_Prevotella | 4 | 41 | disease_enriched(4) | 0.8482927 |
| STR_Tannerella | 4 | 11 | disease_enriched(4) | 0.8990909 |
| STR_Leptotrichia | 2 | 65 | health_enriched(1), disease_enriched(1) | 0.8761538 |
| STR_Capnocytophaga | 2 | 24 | health_enriched(2) | 0.8095833 |
| STR_Rothia | 2 | 1 | health_enriched(2) | 0.9200000 |
| STR_Streptococcus | 1 | 251 | health_enriched(1) | 0.7090040 |
| STR_Selenomonas | 1 | 27 | disease_enriched(1) | 0.6822222 |
| STR_Mogibacterium | 1 | 14 | disease_enriched(1) | 0.8785714 |
| STR_Neisseria | 1 | 14 | health_enriched(1) | 0.6921429 |
| STR_Peptostreptococcus | 1 | 11 | disease_enriched(1) | 0.8409091 |
| STR_Campylobacter | 1 | 10 | disease_enriched(1) | 0.9440000 |

Phage–host genus mapping summary

## Persist mapping artefacts

<details class="code-fold">

<summary>

Code
</summary>

``` r
write_csv(df_strict, path_target("mapping_strict.csv"))
write_csv(df_prefix, path_target("mapping_prefix.csv"))
write_csv(df_genus_summary, path_target("mapping_genus_summary.csv"))
```

</details>

## Quick visualization

<details class="code-fold">

<summary>

Code
</summary>

``` r
df_plot <- df_genus_summary |>
  dplyr::slice_max(n_phages, n = 12) |>
  dplyr::mutate(host = sub("^STR_", "", host_taxon_id))
p <- ggplot(df_plot,
            aes(x = reorder(host, n_phages), y = n_phages)) +
  geom_col(aes(fill = n_candidates)) +
  coord_flip() +
  scale_fill_viridis_c(name = "Candidate species") +
  labs(x = NULL, y = "Phages predicted to infect this genus",
       title = "M2a phage hits per M1a-implicated host genus") +
  theme_minimal(base_size = 11)
ggsave(path_target("fig_genus_mapping.png"),
       p, width = 7, height = 4, dpi = 150)
print(p)
```

</details>

<div id="fig-mapping-genus">

<img
src="090-m2a-mapping_files/figure-commonmark/fig-mapping-genus-1.png"
id="fig-mapping-genus" />

Figure 1: Figure 1

</div>

## Files written

<details class="code-fold">

<summary>

Code
</summary>

``` r
knitr::kable(qproj::proj_dir_info(path_target(), tz = "CET"))
```

</details>

| path                      | type |    size | modification_time   |
|:--------------------------|:-----|--------:|:--------------------|
| fig_genus_mapping.png     | file |  48.77K | 2026-05-05 21:13:31 |
| mapping_genus_summary.csv | file |     925 | 2026-05-05 21:13:31 |
| mapping_prefix.csv        | file | 173.59K | 2026-05-05 21:13:31 |
| mapping_strict.csv        | file |   1.73K | 2026-05-05 21:13:31 |
