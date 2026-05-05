Your Name

- [<span class="toc-section-number">1</span>
  040-m1a-data-processing](#040-m1a-data-processing)
  - [<span class="toc-section-number">1.1</span> Load the cohort
    matrices and metadata](#load-the-cohort-matrices-and-metadata)
  - [<span class="toc-section-number">1.2</span> Build the `sample` rows
    (all six cohorts)](#build-the-sample-rows-all-six-cohorts)
  - [<span class="toc-section-number">1.3</span> Build `taxon_profile`
    (long, non-zero only)](#build-taxon_profile-long-non-zero-only)
  - [<span class="toc-section-number">1.4</span> Write to the contract
    layer
    (partition-based)](#write-to-the-contract-layer-partition-based)
  - [<span class="toc-section-number">1.5</span> Persist
    artefacts](#persist-artefacts)
  - [<span class="toc-section-number">1.6</span> Files
    written](#files-written)

# 040-m1a-data-processing

Your Name 2026-05-05

- [<span class="toc-section-number">1</span> Load the cohort matrices
  and metadata](#load-the-cohort-matrices-and-metadata)
- [<span class="toc-section-number">2</span> Build the `sample` rows
  (all six cohorts)](#build-the-sample-rows-all-six-cohorts)
- [<span class="toc-section-number">3</span> Build `taxon_profile`
  (long, non-zero only)](#build-taxon_profile-long-non-zero-only)
- [<span class="toc-section-number">4</span> Write to the contract layer
  (partition-based)](#write-to-the-contract-layer-partition-based)
- [<span class="toc-section-number">5</span> Persist
  artefacts](#persist-artefacts)
- [<span class="toc-section-number">6</span> Files
  written](#files-written)

**Updated: 2026-05-05 21:12:54 CET.**

Loads the six Geng et al. 2024 modeling bacteriome cohorts (PRJDB11203,
PRJNA230363, PRJNA396840, PRJNA678453, PRJNA717815, PRJNA932553 — 223
saliva / oral samples in total, MetaPhlAn-derived species-level relative
abundances) into the contract layer. Populates `sample` (one row per
sequencing run) and `taxon_profile` (non-zero abundances only, long
format).

The DB writes are **partition-based**: rows whose `source_id` starts
with `geng2024_` (resp. profiles whose `sample_id` starts with
`SAMP_geng2024_`) are deleted before the new batch is appended, so the
Yahara virome rows written by module 070 survive across re-renders.

**Ripple-2** (`040-r2`): add HMT/NCBI taxon ID resolution against eHOMD
instead of the slug fallback.

<details class="code-fold">

<summary>

Code
</summary>

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

## Load the cohort matrices and metadata

<details class="code-fold">

<summary>

Code
</summary>

``` r
path_geng_rdata <- here::here("data", "00-raw", "d020-data-source-registry",
                              "wei2024-imeta-repo", "Rdata", "Figure 01.Rdata")
stopifnot(file.exists(path_geng_rdata))
load(path_geng_rdata)  # creates `my_data` (list with feat_list, meta_list)

project_ids <- intersect(names(my_data$feat_list), names(my_data$meta_list))
stopifnot(length(project_ids) >= 1)
cat("Cohorts loaded from Geng Rdata:", paste(project_ids, collapse = ", "), "\n")
```

</details>

    Cohorts loaded from Geng Rdata: PRJDB11203, PRJNA230363, PRJNA396840, PRJNA678453, PRJNA717815, PRJNA932553 

<details class="code-fold">

<summary>

Code
</summary>

``` r
df_cohort_dims <- purrr::map_dfr(project_ids, function(pid) {
  m <- my_data$feat_list[[pid]]
  d <- my_data$meta_list[[pid]]
  tibble::tibble(project_id = pid,
                 n_samples_feat = nrow(m), n_species = ncol(m),
                 n_samples_meta = nrow(d))
})
knitr::kable(df_cohort_dims, caption = "Per-cohort matrix / metadata dims")
```

</details>

| project_id  | n_samples_feat | n_species | n_samples_meta |
|:------------|---------------:|----------:|---------------:|
| PRJDB11203  |             42 |       230 |             42 |
| PRJNA230363 |             28 |       298 |             28 |
| PRJNA396840 |             20 |       219 |             20 |
| PRJNA678453 |             59 |       349 |             59 |
| PRJNA717815 |             26 |       310 |             26 |
| PRJNA932553 |             48 |       345 |             48 |

Per-cohort matrix / metadata dims

## Build the `sample` rows (all six cohorts)

<details class="code-fold">

<summary>

Code
</summary>

``` r
slugify <- function(x) {
  x |>
    iconv(to = "ASCII//TRANSLIT") |>
    gsub("[^A-Za-z0-9]+", "_", x = _) |>
    gsub("_+$", "", x = _) |>
    gsub("^_+", "", x = _)
}

build_sample_one <- function(pid) {
  df_meta_raw <- my_data$meta_list[[pid]]
  df_meta_raw |>
    tibble::rownames_to_column("run_id") |>
    dplyr::mutate(
      sample_id = paste0("SAMP_geng2024_", run_id),
      source_id = paste0("geng2024_", pid),
      subject_id = run_id,
      disease_status = dplyr::case_when(
        Group == "Case"    ~ "periodontitis",
        Group == "Control" ~ "healthy",
        TRUE               ~ "other"
      ),
      body_site = dplyr::if_else(
        is.na(Bodysite) | Bodysite == "", "saliva", as.character(Bodysite)
      ),
      seq_type  = "shotgun",
      meta_json = purrr::pmap_chr(
        list(country, sex, host_age, BMI, disease_stage, Bodysite),
        function(country, sex, host_age, BMI, disease_stage, Bodysite) {
          jsonlite::toJSON(
            list(country = country, sex = sex, host_age = host_age,
                 BMI = BMI, disease_stage = disease_stage,
                 bodysite = Bodysite),
            auto_unbox = TRUE, na = "null"
          )
        }
      )
    ) |>
    dplyr::select(sample_id, source_id, subject_id, disease_status,
                  body_site, seq_type, meta_json)
}

df_sample <- purrr::map_dfr(project_ids, build_sample_one)
stopifnot(!any(duplicated(df_sample$sample_id)))

dplyr::count(df_sample, source_id, disease_status) |>
  tidyr::pivot_wider(names_from = disease_status, values_from = n,
                     values_fill = 0L) |>
  knitr::kable(caption = "Sample disease status distribution per cohort")
```

</details>

| source_id            | healthy | periodontitis |
|:---------------------|--------:|--------------:|
| geng2024_PRJDB11203  |      19 |            23 |
| geng2024_PRJNA230363 |      18 |            10 |
| geng2024_PRJNA396840 |      10 |            10 |
| geng2024_PRJNA678453 |      29 |            30 |
| geng2024_PRJNA717815 |      12 |            14 |
| geng2024_PRJNA932553 |      18 |            30 |

Sample disease status distribution per cohort

## Build `taxon_profile` (long, non-zero only)

<details class="code-fold">

<summary>

Code
</summary>

``` r
build_profile_one <- function(pid) {
  mat_abund <- my_data$feat_list[[pid]]
  mat_abund |>
    as.data.frame() |>
    tibble::rownames_to_column("run_id") |>
    tidyr::pivot_longer(-run_id, names_to = "taxon_name",
                        values_to = "abundance") |>
    dplyr::filter(abundance > 0) |>
    dplyr::mutate(
      sample_id      = paste0("SAMP_geng2024_", run_id),
      taxon_id       = paste0("STR_", slugify(taxon_name)),
      taxon_kind     = "bacterium",
      abundance_kind = "rel_abund"
    )
}

df_long <- purrr::map_dfr(project_ids, build_profile_one)

df_taxon_profile <- df_long |>
  dplyr::mutate(
    profile_id = paste0(
      "PROF_",
      vapply(paste0(sample_id, "|", taxon_id),
             function(s) substr(digest::digest(s, algo = "xxhash64"), 1, 12),
             character(1))
    )
  ) |>
  dplyr::select(profile_id, sample_id, taxon_id, taxon_kind, taxon_name,
                abundance, abundance_kind)

cat("taxon_profile rows:", nrow(df_taxon_profile),
    " | unique taxa:", dplyr::n_distinct(df_taxon_profile$taxon_id),
    " | unique samples:", dplyr::n_distinct(df_taxon_profile$sample_id), "\n")
```

</details>

    taxon_profile rows: 28350  | unique taxa: 408  | unique samples: 223 

<details class="code-fold">

<summary>

Code
</summary>

``` r
stopifnot(!any(duplicated(df_taxon_profile$profile_id)))

df_taxon_profile |>
  dplyr::count(sample_id) |>
  dplyr::summarise(
    median_per_sample = stats::median(n),
    min_per_sample = min(n),
    max_per_sample = max(n),
    n_samples = dplyr::n()
  ) |>
  knitr::kable(caption = "Per-sample non-zero species counts")
```

</details>

| median_per_sample | min_per_sample | max_per_sample | n_samples |
|------------------:|---------------:|---------------:|----------:|
|               127 |              9 |            249 |       223 |

Per-sample non-zero species counts

## Write to the contract layer (partition-based)

We delete only the geng2024 partition before appending, so module 070’s
Yahara virome rows in `sample` and `taxon_profile` are preserved.

<details class="code-fold">

<summary>

Code
</summary>

``` r
con <- load_db()
assert_schema(con)

# Sanity-check: each cohort must already be registered in data_source
n_src <- DBI::dbGetQuery(
  con,
  "SELECT COUNT(*) AS n FROM data_source WHERE source_id LIKE 'geng2024_%'"
)$n
stopifnot(n_src == length(project_ids))

DBI::dbExecute(con,
  "DELETE FROM sample WHERE source_id LIKE 'geng2024_%'")
```

</details>

    [1] 42

<details class="code-fold">

<summary>

Code
</summary>

``` r
DBI::dbWriteTable(con, "sample", df_sample, append = TRUE,
                  row.names = FALSE)

DBI::dbExecute(con,
  "DELETE FROM taxon_profile WHERE sample_id LIKE 'SAMP_geng2024_%'")
```

</details>

    [1] 1912

<details class="code-fold">

<summary>

Code
</summary>

``` r
DBI::dbWriteTable(con, "taxon_profile", df_taxon_profile, append = TRUE,
                  row.names = FALSE)

df_check <- tibble::tibble(
  metric = c(
    "sample (total)", "sample (geng2024_*)", "sample (yahara2021)",
    "taxon_profile (total)", "taxon_profile (bacterium)",
    "taxon_profile (virus)"
  ),
  rows = c(
    DBI::dbGetQuery(con, "SELECT COUNT(*) AS n FROM sample")$n,
    DBI::dbGetQuery(con, "SELECT COUNT(*) AS n FROM sample WHERE source_id LIKE 'geng2024_%'")$n,
    DBI::dbGetQuery(con, "SELECT COUNT(*) AS n FROM sample WHERE source_id = 'yahara2021'")$n,
    DBI::dbGetQuery(con, "SELECT COUNT(*) AS n FROM taxon_profile")$n,
    DBI::dbGetQuery(con, "SELECT COUNT(*) AS n FROM taxon_profile WHERE taxon_kind = 'bacterium'")$n,
    DBI::dbGetQuery(con, "SELECT COUNT(*) AS n FROM taxon_profile WHERE taxon_kind = 'virus'")$n
  )
)
close_db(con)
knitr::kable(df_check, caption = "Contract tables after 040 write")
```

</details>

| metric                    |  rows |
|:--------------------------|------:|
| sample (total)            |   227 |
| sample (geng2024\_\*)     |   223 |
| sample (yahara2021)       |     4 |
| taxon_profile (total)     | 30063 |
| taxon_profile (bacterium) | 28350 |
| taxon_profile (virus)     |  1713 |

Contract tables after 040 write

## Persist artefacts

<details class="code-fold">

<summary>

Code
</summary>

``` r
write_csv(df_sample, path_target("sample.csv"))
write_csv(head(df_taxon_profile, 200), path_target("taxon_profile_head200.csv"))
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

| path                      | type |  size | modification_time   |
|:--------------------------|:-----|------:|:--------------------|
| sample.csv                | file | 43.5K | 2026-05-05 21:12:56 |
| taxon_profile_head200.csv | file | 26.3K | 2026-05-05 21:12:56 |
