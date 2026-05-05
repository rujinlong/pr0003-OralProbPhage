Your Name

- [<span class="toc-section-number">1</span>
  020-data-source-registry](#020-data-source-registry)
  - [<span class="toc-section-number">1.1</span> Verify raw inputs are
    present](#verify-raw-inputs-are-present)
  - [<span class="toc-section-number">1.2</span> Build the registry
    rows](#build-the-registry-rows)
  - [<span class="toc-section-number">1.3</span> Write to the contract
    layer](#write-to-the-contract-layer)
  - [<span class="toc-section-number">1.4</span> Persist the registry
    table for human
    reference](#persist-the-registry-table-for-human-reference)
  - [<span class="toc-section-number">1.5</span> Files
    written](#files-written)

# 020-data-source-registry

Your Name 2026-05-05

- [<span class="toc-section-number">1</span> Verify raw inputs are
  present](#verify-raw-inputs-are-present)
- [<span class="toc-section-number">2</span> Build the registry
  rows](#build-the-registry-rows)
- [<span class="toc-section-number">3</span> Write to the contract
  layer](#write-to-the-contract-layer)
- [<span class="toc-section-number">4</span> Persist the registry table
  for human reference](#persist-the-registry-table-for-human-reference)
- [<span class="toc-section-number">5</span> Files
  written](#files-written)

**Updated: 2026-05-05 21:13:17 CET.**

Registers the public data sources backing the pipeline. Each source is
one row in `data_source`. Raw supplementary tables live under
`data/00-raw/d020-data-source-registry/<source-key>/`.

**Ripple-1 scope** (`020-r1`) — seven sources: all six Geng et
al. modeling bacteriome cohorts plus the Yahara virome.

- `geng2024_<bioproject>` (×6) — periodontitis bacteriome cohorts from
  the Geng et al. 2024 meta-analysis (**2024-Universal_Geng?**);
  cohort-level facts (n.case / n.control, country, instrument) are
  pulled from sheet `Table S1` of the supplementary table.
- `yahara2021` — long-read oral phageome from PromethION saliva
  sequencing (Yahara et al. 2021). 4 saliva samples (one healthy donor,
  four conditions). VIRSorter-derived phage contigs.

Ripple-2 (`020-r2`) will pair each Geng cohort with a body-site- and
disease-matched virome dataset.

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

## Verify raw inputs are present

The supplementary tables were fetched into
`data/00-raw/d020-data-source-registry/`. The render fails fast if
anything is missing.

<details class="code-fold">

<summary>

Code
</summary>

``` r
path_geng_repo <- here::here(path_data, "wei2024-imeta-repo")
path_geng_supp <- here::here(path_geng_repo, "Supplementary table",
                             "Supplementary Table.xlsx")
path_geng_rdata <- here::here(path_geng_repo, "Rdata", "Figure 01.Rdata")
path_yahara_dir <- here::here(path_data, "yahara2021")
path_yahara_supp4 <- here::here(path_yahara_dir, "MOESM7.xlsx")  # Supp Data 4 = viral seqs
path_yahara_supp5 <- here::here(path_yahara_dir, "MOESM8.xlsx")  # Supp Data 5 = Strep phages

stopifnot(
  file.exists(path_geng_supp),
  file.exists(path_geng_rdata),
  file.exists(path_yahara_supp4),
  file.exists(path_yahara_supp5)
)

df_inputs <- tibble::tibble(
  source = c("geng2024", "geng2024", "yahara2021", "yahara2021"),
  file   = c(basename(path_geng_supp), basename(path_geng_rdata),
             basename(path_yahara_supp4), basename(path_yahara_supp5)),
  bytes  = c(file.size(path_geng_supp), file.size(path_geng_rdata),
             file.size(path_yahara_supp4), file.size(path_yahara_supp5))
)
knitr::kable(df_inputs, caption = "Raw supplementary inputs (verified present)")
```

</details>

| source     | file                     |  bytes |
|:-----------|:-------------------------|-------:|
| geng2024   | Supplementary Table.xlsx | 109597 |
| geng2024   | Figure 01.Rdata          | 617042 |
| yahara2021 | MOESM7.xlsx              | 696905 |
| yahara2021 | MOESM8.xlsx              |  18526 |

Raw supplementary inputs (verified present)

## Build the registry rows

The Geng cohorts are read from `Table S1` of the supplementary
spreadsheet. We keep only the six modeling cohorts (the two validation
cohorts are not shipped in the `Figure 01.Rdata` feat_list and are out
of ripple-1 scope).

<details class="code-fold">

<summary>

Code
</summary>

``` r
df_table_s1 <- readxl::read_excel(path_geng_supp, sheet = "Table S1") |>
  dplyr::filter(`modeling or validation` == "modeling")

geng_table_url <- "https://github.com/whchenlab/Oral-Microbiome-Based-Signature-for-Periodontitis"

df_geng <- df_table_s1 |>
  dplyr::transmute(
    source_id        = paste0("geng2024_", project_id),
    citation_key     = "2024-Universal_Geng",
    doi              = "10.1002/imt2.212",
    source_kind      = "bacteriome",
    n_samples        = as.integer(num.case + num.control),
    table_url        = geng_table_url,
    notes            = sprintf(
      "Ripple-1 M1a; %s saliva/oral cohort %s. %d case / %d control. %s.",
      country_fullname, project_id, num.case, num.control, instrument_model
    )
  )

df_yahara <- tibble::tibble(
  source_id    = "yahara2021",
  citation_key = "2021-Longread_Yahara",
  doi          = "10.1038/s41467-020-20199-9",
  source_kind  = "virome",
  n_samples    = 4L,
  table_url    = "https://www.nature.com/articles/s41467-020-20199-9#Sec19",
  notes        = "MVP M2a; long-read PromethION saliva phageome. 1 healthy donor x 4 conditions. VIRSorter-derived contigs."
)

df_sources <- dplyr::bind_rows(df_geng, df_yahara)
df_sources
```

</details>

| source_id | citation_key | doi | source_kind | n_samples | table_url | notes |
|:---|:---|:---|:---|---:|:---|:---|
| geng2024_PRJDB11203 | 2024-Universal_Geng | 10.1002/imt2.212 | bacteriome | 42 | https://github.com/whchenlab/Oral-Microbiome-Based-Signature-for-Periodontitis | Ripple-1 M1a; Japan saliva/oral cohort PRJDB11203. 23 case / 19 control. Illumina MiSeq. |
| geng2024_PRJNA230363 | 2024-Universal_Geng | 10.1002/imt2.212 | bacteriome | 28 | https://github.com/whchenlab/Oral-Microbiome-Based-Signature-for-Periodontitis | Ripple-1 M1a; People’s Republic of China saliva/oral cohort PRJNA230363. 10 case / 18 control. Illumina HiSeq 2000. |
| geng2024_PRJNA396840 | 2024-Universal_Geng | 10.1002/imt2.212 | bacteriome | 20 | https://github.com/whchenlab/Oral-Microbiome-Based-Signature-for-Periodontitis | Ripple-1 M1a; Denmark saliva/oral cohort PRJNA396840. 10 case / 10 control. Illumina HiSeq 2500. |
| geng2024_PRJNA678453 | 2024-Universal_Geng | 10.1002/imt2.212 | bacteriome | 59 | https://github.com/whchenlab/Oral-Microbiome-Based-Signature-for-Periodontitis | Ripple-1 M1a; Denmark saliva/oral cohort PRJNA678453. 30 case / 29 control. Illumina HiSeq 2500. |
| geng2024_PRJNA717815 | 2024-Universal_Geng | 10.1002/imt2.212 | bacteriome | 26 | https://github.com/whchenlab/Oral-Microbiome-Based-Signature-for-Periodontitis | Ripple-1 M1a; Brazil saliva/oral cohort PRJNA717815. 14 case / 12 control. Illumina HiSeq 4000. |
| geng2024_PRJNA932553 | 2024-Universal_Geng | 10.1002/imt2.212 | bacteriome | 48 | https://github.com/whchenlab/Oral-Microbiome-Based-Signature-for-Periodontitis | Ripple-1 M1a; People’s Republic of China saliva/oral cohort PRJNA932553. 30 case / 18 control. Illumina NovaSeq 6000. |
| yahara2021 | 2021-Longread_Yahara | 10.1038/s41467-020-20199-9 | virome | 4 | https://www.nature.com/articles/s41467-020-20199-9#Sec19 | MVP M2a; long-read PromethION saliva phageome. 1 healthy donor x 4 conditions. VIRSorter-derived contigs. |

## Write to the contract layer

<details class="code-fold">

<summary>

Code
</summary>

``` r
con <- load_db()
assert_schema(con)
write_table_db(con, "data_source", df_sources, append = FALSE)
df_after <- read_table_db(con, "data_source")
close_db(con)
knitr::kable(df_after, caption = "data_source after registry write")
```

</details>

| source_id | citation_key | doi | source_kind | n_samples | table_url | notes |
|:---|:---|:---|:---|---:|:---|:---|
| geng2024_PRJDB11203 | 2024-Universal_Geng | 10.1002/imt2.212 | bacteriome | 42 | https://github.com/whchenlab/Oral-Microbiome-Based-Signature-for-Periodontitis | Ripple-1 M1a; Japan saliva/oral cohort PRJDB11203. 23 case / 19 control. Illumina MiSeq. |
| geng2024_PRJNA230363 | 2024-Universal_Geng | 10.1002/imt2.212 | bacteriome | 28 | https://github.com/whchenlab/Oral-Microbiome-Based-Signature-for-Periodontitis | Ripple-1 M1a; People’s Republic of China saliva/oral cohort PRJNA230363. 10 case / 18 control. Illumina HiSeq 2000. |
| geng2024_PRJNA396840 | 2024-Universal_Geng | 10.1002/imt2.212 | bacteriome | 20 | https://github.com/whchenlab/Oral-Microbiome-Based-Signature-for-Periodontitis | Ripple-1 M1a; Denmark saliva/oral cohort PRJNA396840. 10 case / 10 control. Illumina HiSeq 2500. |
| geng2024_PRJNA678453 | 2024-Universal_Geng | 10.1002/imt2.212 | bacteriome | 59 | https://github.com/whchenlab/Oral-Microbiome-Based-Signature-for-Periodontitis | Ripple-1 M1a; Denmark saliva/oral cohort PRJNA678453. 30 case / 29 control. Illumina HiSeq 2500. |
| geng2024_PRJNA717815 | 2024-Universal_Geng | 10.1002/imt2.212 | bacteriome | 26 | https://github.com/whchenlab/Oral-Microbiome-Based-Signature-for-Periodontitis | Ripple-1 M1a; Brazil saliva/oral cohort PRJNA717815. 14 case / 12 control. Illumina HiSeq 4000. |
| geng2024_PRJNA932553 | 2024-Universal_Geng | 10.1002/imt2.212 | bacteriome | 48 | https://github.com/whchenlab/Oral-Microbiome-Based-Signature-for-Periodontitis | Ripple-1 M1a; People’s Republic of China saliva/oral cohort PRJNA932553. 30 case / 18 control. Illumina NovaSeq 6000. |
| yahara2021 | 2021-Longread_Yahara | 10.1038/s41467-020-20199-9 | virome | 4 | https://www.nature.com/articles/s41467-020-20199-9#Sec19 | MVP M2a; long-read PromethION saliva phageome. 1 healthy donor x 4 conditions. VIRSorter-derived contigs. |

data_source after registry write

## Persist the registry table for human reference

<details class="code-fold">

<summary>

Code
</summary>

``` r
write_csv(df_sources, path_target("data_source.csv"))
```

</details>

## Files written

These files have been written to the target directory,
data/020-data-source-registry:

<details class="code-fold">

<summary>

Code
</summary>

``` r
knitr::kable(qproj::proj_dir_info(path_target(), tz = "CET"))
```

</details>

| path            | type |  size | modification_time   |
|:----------------|:-----|------:|:--------------------|
| data_source.csv | file | 1.77K | 2026-05-05 21:13:18 |

<div id="refs" class="references csl-bib-body hanging-indent">

<div id="ref-2021-Longread_Yahara" class="csl-entry">

Yahara, Koji, Masato Suzuki, Aki Hirabayashi, et al. 2021. “Long-Read
Metagenomics Using PromethION Uncovers Oral Bacteriophages and Their
Interaction with Host Bacteria.” *Nature Communications* 12 (1): 27.
<https://doi.org/10.1038/s41467-020-20199-9>.

</div>

</div>
