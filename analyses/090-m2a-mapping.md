# 090-m2a-mapping
Your Name
2026-05-05

- [<span class="toc-section-number">1</span> Pull contract
  tables](#pull-contract-tables)
- [<span class="toc-section-number">2</span> Strategy 1: strict
  equality](#strategy-1-strict-equality)
- [<span class="toc-section-number">3</span> Strategy 2: genus-prefix
  join (MVP primary)](#strategy-2-genus-prefix-join-mvp-primary)
- [<span class="toc-section-number">4</span> Genus-level fallback
  (090-r1)](#genus-level-fallback-090-r1)
- [<span class="toc-section-number">5</span> Per-genus
  summary](#per-genus-summary)
- [<span class="toc-section-number">6</span> Persist mapping
  artefacts](#persist-mapping-artefacts)
- [<span class="toc-section-number">7</span> Quick
  visualization](#quick-visualization)
- [<span class="toc-section-number">8</span> Files
  written](#files-written)

**Updated: 2026-05-05 22:45:40 CET.**

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

- `090-r1` (this version): add a **genus-token equality fallback** layer
  on top of the species/prefix join. Any candidate whose
  `STR_<Genus>_<...>` shares a genus with a phage host’s
  `STR_<Genus>_<...>` produces a (lower-confidence) mapping, annotated
  with `match_level = 'genus'`. Rows already covered by the
  species/prefix join are kept as `match_level = 'species'`; the genus
  layer only fills gaps. The output feeds 100-knowledge-integration via
  `mapping_with_fallback.tsv`.
- `090-r2`: replace string-based genus extraction with eHOMD-aware fuzzy
  taxonomy (handle synonyms, basonyms, misspellings, Candidatus naming);
  bipartite-graph visualizations (heatmap, network); centrality metrics;
  partitioning of phage–host modules.

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
df_candidate <- read_table_db(con, "candidate_microbe")
df_link      <- read_table_db(con, "phage_host_link")
close_db(con)
cat("candidate_microbe:", nrow(df_candidate),
    " | phage_host_link:", nrow(df_link), "\n")
```

</details>

    candidate_microbe: 764  | phage_host_link: 770 

## Strategy 1: strict equality

<details class="code-fold">
<summary>Code</summary>

``` r
df_strict <- df_candidate |>
  dplyr::inner_join(df_link, by = c("taxon_id" = "host_taxon_id")) |>
  dplyr::select(candidate_id, taxon_id, taxon_name, phenotype, direction,
                effect_size, q_value,
                link_id, vOTU_id, evidence_kind, confidence)
cat("Strict-equality hits:", nrow(df_strict), "\n")
```

</details>

    Strict-equality hits: 100 

## Strategy 2: genus-prefix join (MVP primary)

<details class="code-fold">
<summary>Code</summary>

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

    Genus-prefix hits: 16760  | unique candidates with phage hits: 430  | unique phages: 649 

## Genus-level fallback (090-r1)

The species/prefix join above (`df_prefix`) catches the easy case:
candidate `STR_Streptococcus_salivarius` × host `STR_Streptococcus`. But
it misses pairs where the host token tree is more specific than just the
genus (e.g. host `STR_Solobacterium_moorei` × candidate
`STR_Solobacterium_*`), and pairs where both sides are at the same depth
but neither prefixes the other (e.g. host `STR_Prevotella_shahii` ×
candidate `STR_Prevotella_denticola`).

The fix is a genus-token equality layer: extract the first non-`STR_`
token on each side and join on it. Wherever the species/prefix join
already produced a row, we keep it as `match_level = 'species'`; the
genus layer only fills the remaining `(candidate_id, link_id)` slots
with `match_level = 'genus'`. Downstream scoring (in 100) discounts
`match_level = 'genus'` rows by a 0.6× factor to reflect the lower
taxonomic resolution.

Edge cases tolerated by the simple regex `^STR_([^_]+).*$`:

- `STR_Candidatus_Saccharimonas_sp` → genus token “Candidatus” (the same
  token is extracted on both sides, so Candidatus-named candidates still
  match Candidatus-named hosts — apples-to-apples).
- `STR_Lachnospiraceae_bacterium` → “Lachnospiraceae” (a family, not
  strictly a genus, but again symmetric across both tables).
- `STR_GGB1611_SGB2208` → “GGB1611” (an unnamed-genome bin); these only
  match other GGB-prefixed tokens, which is the correct behaviour.

<details class="code-fold">
<summary>Code</summary>

``` r
extract_genus <- function(x) sub("^STR_([^_]+).*$", "\\1", x)

df_candidate_g <- df_candidate |>
  dplyr::mutate(candidate_genus = extract_genus(taxon_id))
df_link_g <- df_link |>
  dplyr::mutate(host_genus = extract_genus(host_taxon_id))

df_species_layer <- df_prefix |>
  dplyr::mutate(
    candidate_genus = extract_genus(taxon_id),
    host_genus      = extract_genus(host_taxon_id),
    match_level     = "species"
  ) |>
  dplyr::rename(candidate_taxon_id = taxon_id)

df_genus_all <- df_candidate_g |>
  dplyr::inner_join(df_link_g, by = c("candidate_genus" = "host_genus"),
                    relationship = "many-to-many") |>
  dplyr::mutate(host_genus = candidate_genus)

species_keys <- paste(df_species_layer$candidate_id,
                      df_species_layer$link_id, sep = "|")
df_genus_layer <- df_genus_all |>
  dplyr::filter(!paste(candidate_id, link_id, sep = "|") %in% species_keys) |>
  dplyr::transmute(
    candidate_id, link_id,
    candidate_taxon_id = taxon_id,
    taxon_name, phenotype, direction, effect_size, q_value,
    vOTU_id, host_taxon_id, evidence_kind, confidence,
    candidate_genus, host_genus,
    match_level = "genus"
  )

df_mapping_with_fallback <- dplyr::bind_rows(
  df_species_layer |>
    dplyr::select(candidate_id, link_id, candidate_taxon_id, host_taxon_id,
                  candidate_genus, host_genus, match_level,
                  taxon_name, phenotype, direction, effect_size, q_value,
                  vOTU_id, evidence_kind, confidence),
  df_genus_layer |>
    dplyr::select(candidate_id, link_id, candidate_taxon_id, host_taxon_id,
                  candidate_genus, host_genus, match_level,
                  taxon_name, phenotype, direction, effect_size, q_value,
                  vOTU_id, evidence_kind, confidence)
)

cat("Mapping after fallback union:",
    nrow(df_mapping_with_fallback),
    "rows |  species:",
    sum(df_mapping_with_fallback$match_level == "species"),
    "  genus:",
    sum(df_mapping_with_fallback$match_level == "genus"),
    "\n")
```

</details>

    Mapping after fallback union: 18461 rows |  species: 16760   genus: 1701 

<details class="code-fold">
<summary>Code</summary>

``` r
cat("Distinct candidates with mapping  | species-only:",
    dplyr::n_distinct(df_species_layer$candidate_id),
    "  with-fallback:",
    dplyr::n_distinct(df_mapping_with_fallback$candidate_id), "\n")
```

</details>

    Distinct candidates with mapping  | species-only: 430   with-fallback: 490 

<details class="code-fold">
<summary>Code</summary>

``` r
readr::write_tsv(df_mapping_with_fallback,
                 path_target("mapping_with_fallback.tsv"))
```

</details>

<details class="code-fold">
<summary>Code</summary>

``` r
df_fallback_genus_summary <- df_genus_layer |>
  dplyr::group_by(host_genus) |>
  dplyr::summarise(
    n_new_pairs       = dplyr::n(),
    n_new_candidates  = dplyr::n_distinct(candidate_id),
    n_new_phages      = dplyr::n_distinct(vOTU_id),
    .groups = "drop"
  ) |>
  dplyr::arrange(dplyr::desc(n_new_pairs))
knitr::kable(
  utils::head(df_fallback_genus_summary, 15),
  caption = "Top genera gaining mapping rows from the 090-r1 genus fallback"
)
```

</details>

| host_genus     | n_new_pairs | n_new_candidates | n_new_phages |
|:---------------|------------:|-----------------:|-------------:|
| Candidatus     |         569 |               23 |           25 |
| Prevotella     |         360 |               72 |            5 |
| Streptococcus  |         256 |               35 |            8 |
| Veillonella    |          96 |               12 |            8 |
| Oribacterium   |          70 |               10 |            7 |
| Selenomonas    |          54 |               27 |            2 |
| Solobacterium  |          44 |                2 |           22 |
| Actinomyces    |          43 |               43 |            1 |
| Porphyromonas  |          43 |               15 |            3 |
| Alloprevotella |          26 |                8 |            4 |
| Eubacterium    |          17 |               17 |            1 |
| Catonella      |          14 |                2 |            7 |
| Neisseria      |          13 |               13 |            1 |
| Gemella        |          12 |                6 |            2 |
| Stomatobaculum |          12 |                4 |            3 |

Top genera gaining mapping rows from the 090-r1 genus fallback

<details class="code-fold">
<summary>Code</summary>

``` r
readr::write_tsv(df_fallback_genus_summary,
                 path_target("fallback_genus_summary.tsv"))
```

</details>

## Per-genus summary

<details class="code-fold">
<summary>Code</summary>

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
| STR_Prevotella | 72 | 41 | disease_enriched(50), health_enriched(22) | 0.8482927 |
| STR_Actinomyces | 43 | 6 | health_enriched(32), disease_enriched(11) | 0.7216667 |
| STR_Streptococcus | 35 | 251 | health_enriched(25), disease_enriched(10) | 0.7090040 |
| STR_Treponema | 33 | 1 | disease_enriched(31), health_enriched(2) | 0.9300000 |
| STR_Selenomonas | 27 | 27 | disease_enriched(16), health_enriched(11) | 0.6822222 |
| STR_Capnocytophaga | 24 | 24 | health_enriched(8), disease_enriched(16) | 0.8095833 |
| STR_Leptotrichia | 20 | 65 | health_enriched(6), disease_enriched(14) | 0.8761538 |
| STR_Campylobacter | 17 | 10 | disease_enriched(11), health_enriched(6) | 0.9440000 |
| STR_Porphyromonas | 15 | 2 | disease_enriched(11), health_enriched(4) | 0.5900000 |
| STR_Neisseria | 13 | 14 | health_enriched(9), disease_enriched(4) | 0.6921429 |
| STR_Veillonella | 12 | 81 | health_enriched(8), disease_enriched(4) | 0.9018519 |
| STR_Haemophilus | 12 | 4 | health_enriched(12) | 0.5925000 |
| STR_Tannerella | 10 | 11 | disease_enriched(10) | 0.8990909 |
| STR_Oribacterium | 10 | 10 | health_enriched(8), disease_enriched(2) | 0.7840000 |
| STR_Corynebacterium | 9 | 2 | health_enriched(7), disease_enriched(2) | 0.7400000 |
| STR_Aggregatibacter | 8 | 1 | disease_enriched(6), health_enriched(2) | 0.5100000 |
| STR_Alloprevotella | 8 | 1 | health_enriched(3), disease_enriched(5) | 0.6400000 |
| STR_Rothia | 8 | 1 | health_enriched(8) | 0.9200000 |
| STR_Fusobacterium | 6 | 34 | disease_enriched(4), health_enriched(2) | 0.8488235 |
| STR_Gemella | 6 | 6 | health_enriched(4), disease_enriched(2) | 0.5750000 |
| STR_Candidatus_Gracilibacteria_bacterium | 6 | 1 | health_enriched(4), disease_enriched(2) | 0.6400000 |
| STR_Mogibacterium | 5 | 14 | disease_enriched(3), health_enriched(2) | 0.8785714 |
| STR_Peptostreptococcus | 5 | 11 | disease_enriched(3), health_enriched(2) | 0.8409091 |
| STR_Granulicatella | 4 | 7 | health_enriched(3), disease_enriched(1) | 0.6085714 |
| STR_Streptococcus_sp | 4 | 6 | health_enriched(4) | 0.6516667 |
| STR_Lachnospiraceae_bacterium | 4 | 4 | health_enriched(2), disease_enriched(2) | 0.5800000 |
| STR_Abiotrophia | 4 | 1 | health_enriched(4) | 0.5100000 |
| STR_Bifidobacterium | 4 | 1 | disease_enriched(4) | 0.8100000 |
| STR_Cardiobacterium | 4 | 1 | health_enriched(3), disease_enriched(1) | 0.9200000 |
| STR_Alloprevotella_sp | 2 | 3 | health_enriched(1), disease_enriched(1) | 0.8766667 |
| STR_Megasphaera | 2 | 3 | health_enriched(2) | 0.8366667 |
| STR_Alloscardovia | 2 | 2 | disease_enriched(2) | 0.7450000 |
| STR_Slackia | 2 | 2 | disease_enriched(2) | 0.8600000 |
| STR_Porphyromonas_sp | 2 | 1 | disease_enriched(2) | 0.5400000 |

Phage–host genus mapping summary

## Persist mapping artefacts

<details class="code-fold">
<summary>Code</summary>

``` r
write_csv(df_strict, path_target("mapping_strict.csv"))
write_csv(df_prefix, path_target("mapping_prefix.csv"))
write_csv(df_genus_summary, path_target("mapping_genus_summary.csv"))
```

</details>

## Quick visualization

<details class="code-fold">
<summary>Code</summary>

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

Figure 1

</div>

## Files written

<details class="code-fold">
<summary>Code</summary>

``` r
knitr::kable(qproj::proj_dir_info(path_target(), tz = "CET"))
```

</details>

| path                       | type |   size | modification_time   |
|:---------------------------|:-----|-------:|:--------------------|
| fallback_genus_summary.tsv | file |    594 | 2026-05-05 22:45:42 |
| fig_genus_mapping.png      | file | 47.65K | 2026-05-05 22:45:42 |
| mapping_genus_summary.csv  | file |  2.42K | 2026-05-05 22:45:42 |
| mapping_prefix.csv         | file |  3.54M | 2026-05-05 22:45:42 |
| mapping_strict.csv         | file | 19.94K | 2026-05-05 22:45:42 |
| mapping_with_fallback.tsv  | file |  4.54M | 2026-05-05 22:45:42 |
