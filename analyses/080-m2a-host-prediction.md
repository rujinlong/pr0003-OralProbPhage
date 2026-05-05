# 080-m2a-host-prediction
Your Name
2026-05-05

- [<span class="toc-section-number">1</span> Pull the stashed Yahara
  contig metadata](#pull-the-stashed-yahara-contig-metadata)
- [<span class="toc-section-number">2</span> Parse the CAT taxonomy
  string](#parse-the-cat-taxonomy-string)
- [<span class="toc-section-number">3</span> Restrict to actionable
  ranks (genus / species /
  etc.)](#restrict-to-actionable-ranks-genus--species--etc)
- [<span class="toc-section-number">4</span> Build `phage_host_link`
  rows](#build-phage_host_link-rows)
- [<span class="toc-section-number">5</span> Write to the contract
  layer](#write-to-the-contract-layer)
- [<span class="toc-section-number">6</span> Persist](#persist)
- [<span class="toc-section-number">7</span> Files
  written](#files-written)

**Updated: 2026-05-05 21:13:41 CET.**

Parses the CAT taxonomy column from the Yahara 2021 viral contig table
(stashed by 070) into structured `phage_host_link` rows. CAT assigns
each contig a single best-rank taxonomy string like
`Granulicatella (genus): 0.53`. The numeric tail is the CAT confidence
score (0..1) which we propagate as `confidence`.
`evidence_kind = "published_prediction"` since these calls come from the
paper as-is.

We keep only **genus-level or finer** assignments — superkingdom /
phylum / “no rank” rows are not actionable for host-microbe matching.

**Ripple ladder**:

- `080-r1`: ensemble — re-run host prediction with iPHoP, RaFAH,
  spacepharer CRISPR matches against eHOMD, BLAST against IMG/VR, and
  combine via precision-recall–calibrated voting.
- `080-r2`: deep-learning predictors (vHULK, PHIST) for the long
  contigs.

<details class="code-fold">
<summary>Code</summary>

``` r
suppressPackageStartupMessages({
  library(here)
  library(conflicted)
  library(tidyverse)
  library(jsonlite)
  library(digest)
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

## Pull the stashed Yahara contig metadata

<details class="code-fold">
<summary>Code</summary>

``` r
path_in <- here::here("data", "070-m2a-virome-processing",
                      "phage_meta_for_080.csv")
stopifnot(file.exists(path_in))
df_phage_meta <- readr::read_csv(path_in, show_col_types = FALSE)
cat("Input rows:", nrow(df_phage_meta), "\n")
```

</details>

    Input rows: 1713 

## Parse the CAT taxonomy string

<details class="code-fold">
<summary>Code</summary>

``` r
slugify <- function(x) {
  x |>
    iconv(to = "ASCII//TRANSLIT") |>
    gsub("[^A-Za-z0-9]+", "_", x = _) |>
    gsub("_+$", "", x = _) |>
    gsub("^_+", "", x = _)
}

# CAT format: "<name> (<rank>): <conf>"
re_cat <- "^(.*?)\\s*\\(([^)]+)\\):\\s*([0-9.]+)$"

df_parsed <- df_phage_meta |>
  dplyr::filter(!is.na(cat_taxonomy)) |>
  dplyr::mutate(
    host_name        = sub(re_cat, "\\1", cat_taxonomy),
    host_rank        = sub(re_cat, "\\2", cat_taxonomy),
    confidence       = suppressWarnings(as.numeric(sub(re_cat, "\\3", cat_taxonomy)))
  ) |>
  dplyr::filter(!is.na(confidence))

cat("After parse:", nrow(df_parsed), "\n")
```

</details>

    After parse: 1709 

<details class="code-fold">
<summary>Code</summary>

``` r
cat("Rank distribution:\n")
```

</details>

    Rank distribution:

<details class="code-fold">
<summary>Code</summary>

``` r
print(dplyr::count(df_parsed, host_rank, sort = TRUE))
```

</details>

    # A tibble: 9 × 2
      host_rank        n
      <chr>        <int>
    1 genus          639
    2 superkingdom   379
    3 no rank        167
    4 family         166
    5 species        129
    6 class           84
    7 phylum          80
    8 order           63
    9 subspecies       2

## Restrict to actionable ranks (genus / species / etc.)

<details class="code-fold">
<summary>Code</summary>

``` r
actionable_ranks <- c("genus", "species", "subspecies", "strain")
df_actionable <- df_parsed |>
  dplyr::filter(host_rank %in% actionable_ranks) |>
  dplyr::mutate(
    host_taxon_id = paste0("STR_", slugify(host_name))
  )
cat("Actionable host-link rows:", nrow(df_actionable), "\n")
```

</details>

    Actionable host-link rows: 770 

<details class="code-fold">
<summary>Code</summary>

``` r
print(dplyr::count(df_actionable, host_rank))
```

</details>

    # A tibble: 3 × 2
      host_rank      n
      <chr>      <int>
    1 genus        639
    2 species      129
    3 subspecies     2

## Build `phage_host_link` rows

<details class="code-fold">
<summary>Code</summary>

``` r
ver <- schema_version()

df_phage_host_link <- df_actionable |>
  dplyr::transmute(
    vOTU_id        = vOTU_id,
    host_taxon_id,
    evidence_kind  = "published_prediction",
    confidence,
    support_json   = purrr::pmap_chr(
      list(host_name, host_rank, integrase, virsorter_category, img_vr_cluster),
      function(host_name, host_rank, integrase, virsorter_category, img_vr_cluster) {
        jsonlite::toJSON(
          list(host_name        = host_name,
               host_rank        = host_rank,
               integrase        = integrase,
               virsorter_cat    = virsorter_category,
               img_vr_cluster   = img_vr_cluster,
               method           = "CAT (Yahara 2021)"),
          auto_unbox = TRUE, na = "null"
        )
      }
    ),
    schema_version = ver
  ) |>
  dplyr::mutate(
    link_id = paste0(
      "LINK_",
      vapply(paste(vOTU_id, host_taxon_id, evidence_kind, sep = "|"),
             function(s) substr(digest::digest(s, algo = "xxhash64"), 1, 12),
             character(1))
    )
  ) |>
  dplyr::distinct(link_id, .keep_all = TRUE) |>
  dplyr::select(link_id, vOTU_id, host_taxon_id, evidence_kind,
                confidence, support_json, schema_version)

cat("phage_host_link rows:", nrow(df_phage_host_link),
    "| unique vOTUs:", dplyr::n_distinct(df_phage_host_link$vOTU_id),
    "| unique hosts:", dplyr::n_distinct(df_phage_host_link$host_taxon_id), "\n")
```

</details>

    phage_host_link rows: 770 | unique vOTUs: 770 | unique hosts: 81 

<details class="code-fold">
<summary>Code</summary>

``` r
df_phage_host_link |>
  dplyr::count(host_taxon_id, sort = TRUE) |>
  head(15) |>
  knitr::kable(caption = "Top 15 predicted hosts by phage-link count")
```

</details>

| host_taxon_id                   |   n |
|:--------------------------------|----:|
| STR_Streptococcus               | 251 |
| STR_Veillonella                 |  81 |
| STR_Leptotrichia                |  65 |
| STR_Prevotella                  |  41 |
| STR_Fusobacterium               |  34 |
| STR_Selenomonas                 |  27 |
| STR_Capnocytophaga              |  24 |
| STR_Candidatus_Saccharimonas_sp |  23 |
| STR_Solobacterium_moorei        |  22 |
| STR_Mogibacterium               |  14 |
| STR_Neisseria                   |  14 |
| STR_Peptostreptococcus          |  11 |
| STR_Tannerella                  |  11 |
| STR_Campylobacter               |  10 |
| STR_Oribacterium                |  10 |

Top 15 predicted hosts by phage-link count

## Write to the contract layer

<details class="code-fold">
<summary>Code</summary>

``` r
con <- load_db()
assert_schema(con)
write_table_db(con, "phage_host_link", df_phage_host_link, append = FALSE)
n <- DBI::dbGetQuery(con, "SELECT COUNT(*) AS n FROM phage_host_link")$n
close_db(con)
cat("phage_host_link rows in DB:", n, "\n")
```

</details>

    phage_host_link rows in DB: 770 

## Persist

<details class="code-fold">
<summary>Code</summary>

``` r
write_csv(df_phage_host_link, path_target("phage_host_link.csv"))
```

</details>

## Files written

<details class="code-fold">
<summary>Code</summary>

``` r
knitr::kable(qproj::proj_dir_info(path_target(), tz = "CET"))
```

</details>

| path                | type | size | modification_time   |
|:--------------------|:-----|-----:|:--------------------|
| phage_host_link.csv | file | 335K | 2026-05-05 21:13:42 |
