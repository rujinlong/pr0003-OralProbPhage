Your Name

- [<span class="toc-section-number">1</span>
  070-m2a-virome-processing](#070-m2a-virome-processing)
  - [<span class="toc-section-number">1.1</span> Load Yahara MOESM7
    (Supp Data 4: VIRSorter
    contigs)](#load-yahara-moesm7-supp-data-4-virsorter-contigs)
  - [<span class="toc-section-number">1.2</span> Filter and
    shape](#filter-and-shape)
  - [<span class="toc-section-number">1.3</span> Build `sample` rows for
    the four Yahara saliva
    samples](#build-sample-rows-for-the-four-yahara-saliva-samples)
  - [<span class="toc-section-number">1.4</span> Build `taxon_profile`
    rows](#build-taxon_profile-rows)
  - [<span class="toc-section-number">1.5</span> Stash phage metadata
    for 080 (host
    prediction)](#stash-phage-metadata-for-080-host-prediction)
  - [<span class="toc-section-number">1.6</span> Write to the contract
    layer (append, since 040 already wrote bacterial samples +
    profiles)](#write-to-the-contract-layer-append-since-040-already-wrote-bacterial-samples--profiles)
  - [<span class="toc-section-number">1.7</span> Files
    written](#files-written)

# 070-m2a-virome-processing

Your Name 2026-05-05

- [<span class="toc-section-number">1</span> Load Yahara MOESM7 (Supp
  Data 4: VIRSorter
  contigs)](#load-yahara-moesm7-supp-data-4-virsorter-contigs)
- [<span class="toc-section-number">2</span> Filter and
  shape](#filter-and-shape)
- [<span class="toc-section-number">3</span> Build `sample` rows for the
  four Yahara saliva
  samples](#build-sample-rows-for-the-four-yahara-saliva-samples)
- [<span class="toc-section-number">4</span> Build `taxon_profile`
  rows](#build-taxon_profile-rows)
- [<span class="toc-section-number">5</span> Stash phage metadata for
  080 (host prediction)](#stash-phage-metadata-for-080-host-prediction)
- [<span class="toc-section-number">6</span> Write to the contract layer
  (append, since 040 already wrote bacterial samples +
  profiles)](#write-to-the-contract-layer-append-since-040-already-wrote-bacterial-samples--profiles)
- [<span class="toc-section-number">7</span> Files
  written](#files-written)

**Updated: 2026-05-05 15:17:56 CET.**

Loads the Yahara et al. 2021 long-read oral virome (PromethION saliva
sequencing, 4 samples, VIRSorter-derived contigs). Populates four
`sample` rows (saliva, healthy donor) and `taxon_profile` rows
(`taxon_kind = "virus"`, `abundance_kind = "presence"` since per-sample
read counts are sparse and inconsistent across the four sheets).

VIRSorter category 3 and category 6 contigs are dropped per the paper’s
convention (low-confidence phage/prophage calls).

**Ripple ladder**:

- `070-r1`: re-process the raw PromethION fastq via ViWrap on SLURM, run
  CheckV for completeness/quality flags, dereplicate at 95% ANI,
  regenerate vOTU IDs.
- `070-r2`: AMG annotation (DRAM-v) and lifestyle (lytic vs temperate)
  calls.

<details class="code-fold">

<summary>

Code
</summary>

``` r
suppressPackageStartupMessages({
  library(here)
  library(conflicted)
  library(tidyverse)
  library(readxl)
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

## Load Yahara MOESM7 (Supp Data 4: VIRSorter contigs)

<details class="code-fold">

<summary>

Code
</summary>

``` r
path_xlsx <- here::here("data", "00-raw", "d020-data-source-registry",
                        "yahara2021", "MOESM7.xlsx")
stopifnot(file.exists(path_xlsx))

sheets <- readxl::excel_sheets(path_xlsx)
df_phage_raw <- purrr::map_dfr(sheets, function(s) {
  d <- readxl::read_excel(path_xlsx, sheet = s, skip = 2,
                          .name_repair = "universal_quiet")
  d$sample_label <- s  # Sample1 .. Sample4
  d
})
cat("Total contigs across 4 sheets:", nrow(df_phage_raw), "\n")
```

</details>

    Total contigs across 4 sheets: 5117 

<details class="code-fold">

<summary>

Code
</summary>

``` r
print(dplyr::count(df_phage_raw, sample_label, name = "n_contigs"))
```

</details>

    # A tibble: 4 × 2
      sample_label n_contigs
      <chr>            <int>
    1 Sample1           1335
    2 Sample2           1200
    3 Sample3           1718
    4 Sample4            864

## Filter and shape

<details class="code-fold">

<summary>

Code
</summary>

``` r
df_phage <- df_phage_raw |>
  dplyr::filter(!is.na(phageID),
                category..VirSorter. %in% c(1, 2, 4, 5)) |>
  dplyr::mutate(
    sample_id      = paste0("SAMP_yahara2021_", sample_label),
    vOTU_index     = sprintf("%05d", dplyr::row_number()),
    taxon_id       = paste0("vOTU_yahara2021_", vOTU_index),
    taxon_kind     = "virus",
    taxon_name     = phageID,
    abundance      = 1.0,
    abundance_kind = "presence"
  )
cat("After category filter:", nrow(df_phage),
    "| samples:", dplyr::n_distinct(df_phage$sample_id),
    "| unique vOTU IDs:", dplyr::n_distinct(df_phage$taxon_id), "\n")
```

</details>

    After category filter: 1713 | samples: 4 | unique vOTU IDs: 1713 

## Build `sample` rows for the four Yahara saliva samples

<details class="code-fold">

<summary>

Code
</summary>

``` r
df_sample_yahara <- tibble::tibble(
  sample_id      = paste0("SAMP_yahara2021_", sheets),
  source_id      = "yahara2021",
  subject_id     = "yahara2021_donor",  # one donor across the four conditions
  disease_status = "healthy",
  body_site      = "saliva",
  seq_type       = "virome",
  meta_json      = jsonlite::toJSON(
    list(platform = "PromethION + HiSeq",
         note     = "single donor x 4 conditions"),
    auto_unbox = TRUE
  ) |> rep(length(sheets))
)
df_sample_yahara
```

</details>

| sample_id | source_id | subject_id | disease_status | body_site | seq_type | meta_json |
|:---|:---|:---|:---|:---|:---|:---|
| SAMP_yahara2021_Sample1 | yahara2021 | yahara2021_donor | healthy | saliva | virome | {“platform”:“PromethION + HiSeq”,“note”:“single donor x 4 conditions”} |
| SAMP_yahara2021_Sample2 | yahara2021 | yahara2021_donor | healthy | saliva | virome | {“platform”:“PromethION + HiSeq”,“note”:“single donor x 4 conditions”} |
| SAMP_yahara2021_Sample3 | yahara2021 | yahara2021_donor | healthy | saliva | virome | {“platform”:“PromethION + HiSeq”,“note”:“single donor x 4 conditions”} |
| SAMP_yahara2021_Sample4 | yahara2021 | yahara2021_donor | healthy | saliva | virome | {“platform”:“PromethION + HiSeq”,“note”:“single donor x 4 conditions”} |

## Build `taxon_profile` rows

<details class="code-fold">

<summary>

Code
</summary>

``` r
df_taxon_profile_v <- df_phage |>
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
stopifnot(!any(duplicated(df_taxon_profile_v$profile_id)))
cat("virus taxon_profile rows:", nrow(df_taxon_profile_v), "\n")
```

</details>

    virus taxon_profile rows: 1713 

## Stash phage metadata for 080 (host prediction)

The full phage table — including the CAT-assigned host taxonomy — is the
input to module 080. Persist it as a CSV that 080 will read.

<details class="code-fold">

<summary>

Code
</summary>

``` r
df_phage_for_080 <- df_phage |>
  dplyr::transmute(
    sample_id,
    vOTU_id              = taxon_id,
    phageID,
    phageLen,
    linear_or_circular   = linear.or.circular,
    integrase            = Integrase,
    virsorter_category   = category..VirSorter.,
    cat_taxonomy         = taxonomy.assigned.to.the.contig.by.CAT,
    img_vr_cluster       = clustered.with.IMG.VR.v2.0,
    assigned_family      = assigned.family
  )
write_csv(df_phage_for_080, path_target("phage_meta_for_080.csv"))
```

</details>

## Write to the contract layer (append, since 040 already wrote bacterial samples + profiles)

<details class="code-fold">

<summary>

Code
</summary>

``` r
con <- load_db()
assert_schema(con)

# Append the 4 Yahara samples (don't truncate the 42 Geng samples 040 wrote)
DBI::dbExecute(con,
  "DELETE FROM sample WHERE source_id = 'yahara2021'")
```

</details>

    [1] 0

<details class="code-fold">

<summary>

Code
</summary>

``` r
DBI::dbWriteTable(con, "sample", df_sample_yahara, append = TRUE,
                  row.names = FALSE)

DBI::dbExecute(con,
  "DELETE FROM taxon_profile WHERE taxon_kind = 'virus'")
```

</details>

    [1] 0

<details class="code-fold">

<summary>

Code
</summary>

``` r
DBI::dbWriteTable(con, "taxon_profile", df_taxon_profile_v, append = TRUE,
                  row.names = FALSE)

df_check <- tibble::tibble(
  table = c("sample", "taxon_profile_bacteria", "taxon_profile_virus"),
  rows  = c(
    DBI::dbGetQuery(con, "SELECT COUNT(*) AS n FROM sample")$n,
    DBI::dbGetQuery(con, "SELECT COUNT(*) AS n FROM taxon_profile WHERE taxon_kind='bacterium'")$n,
    DBI::dbGetQuery(con, "SELECT COUNT(*) AS n FROM taxon_profile WHERE taxon_kind='virus'")$n
  )
)
close_db(con)
knitr::kable(df_check, caption = "Contract tables after 070 write")
```

</details>

| table                  | rows |
|:-----------------------|-----:|
| sample                 |   46 |
| taxon_profile_bacteria | 1912 |
| taxon_profile_virus    | 1713 |

Contract tables after 070 write

## Files written

<details class="code-fold">

<summary>

Code
</summary>

``` r
knitr::kable(qproj::proj_dir_info(path_target(), tz = "CET"))
```

</details>

| path                   | type | size | modification_time   |
|:-----------------------|:-----|-----:|:--------------------|
| phage_meta_for_080.csv | file | 532K | 2026-05-05 15:17:57 |
