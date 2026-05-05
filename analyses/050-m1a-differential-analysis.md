# 050-m1a-differential-analysis
Your Name
2026-05-05

- [<span class="toc-section-number">1</span> Pull samples + abundances
  from the contract
  layer](#pull-samples--abundances-from-the-contract-layer)
- [<span class="toc-section-number">2</span> Wilcoxon rank-sum (peri vs
  healthy)](#wilcoxon-rank-sum-peri-vs-healthy)
- [<span class="toc-section-number">3</span> Build Wilcoxon
  `candidate_microbe` rows (MVP baseline,
  preserved)](#build-wilcoxon-candidate_microbe-rows-mvp-baseline-preserved)
- [<span class="toc-section-number">4</span> Helper: rows in the long
  candidate format](#helper-rows-in-the-long-candidate-format)
- [<span class="toc-section-number">5</span> MaAsLin2](#maaslin2)
- [<span class="toc-section-number">6</span> ANCOMBC](#ancombc)
- [<span class="toc-section-number">7</span> DESeq2](#deseq2)
- [<span class="toc-section-number">8</span> Wilcoxon → method-tagged
  rows (full results)](#wilcoxon--method-tagged-rows-full-results)
- [<span class="toc-section-number">9</span> Method agreement
  diagnostics](#method-agreement-diagnostics)
- [<span class="toc-section-number">10</span> Consensus candidate
  set](#consensus-candidate-set)
- [<span class="toc-section-number">11</span> Combine all methods +
  write to the contract
  layer](#combine-all-methods--write-to-the-contract-layer)
- [<span class="toc-section-number">12</span> Persist](#persist)
- [<span class="toc-section-number">13</span> Cross-cohort meta-analysis
  (050-r2)](#cross-cohort-meta-analysis-050-r2)
- [<span class="toc-section-number">14</span> Files
  written](#files-written)

**Updated: 2026-05-05 22:37:53 CET.**

Multi-method, single-cohort differential abundance analysis on the Geng
PRJDB11203 saliva sub-cohort. The MVP baseline (Wilcoxon rank-sum +
Benjamini–Hochberg FDR on TSS-normalized relative abundances) is
preserved verbatim. Ripple `050-r1` adds three additional methods
(MaAsLin2, ANCOMBC, DESeq2) and emits a consensus candidate set
requiring agreement (q \< 0.1 + consistent direction) in ≥ 2 of the 4
methods. Output goes into `candidate_microbe` (phenotype =
“periodontitis”).

**Ripple ladder**:

- `050-r1` (this revision): MaAsLin2 + ANCOMBC + DESeq2 alongside
  Wilcoxon, with consensus requiring q \< 0.1 in ≥ 2 of 4 methods +
  matching direction.
- `050-r2`: cross-cohort meta-analysis (random-effects) over the seven
  Geng cohorts; emit candidates with `n_cohorts ≥ 3`.
- `050-r3`: pathway-level candidates from HUMAnN3 in addition to
  species.

<details class="code-fold">
<summary>Code</summary>

``` r
suppressPackageStartupMessages({
  library(here)
  library(conflicted)
  library(tidyverse)
  library(digest)
  library(DBI)
  library(Maaslin2)
  library(ANCOMBC)
  library(DESeq2)
  library(phyloseq)
  devtools::load_all()
})
conflicted::conflicts_prefer(
  dplyr::filter,
  dplyr::lag,
  dplyr::first,
  dplyr::last,
  dplyr::between,
  dplyr::desc,
  dplyr::collapse,
  dplyr::count,
  dplyr::rename,
  dplyr::select,
  dplyr::slice,
  dplyr::combine,
  base::setdiff,
  base::intersect,
  base::union,
  base::Reduce,
  base::Position,
  base::Find,
  base::Filter,
  base::Map,
  base::strsplit,
  stats::IQR,
  base::which,
  .quiet = TRUE
)
```

</details>

## Pull samples + abundances from the contract layer

<details class="code-fold">
<summary>Code</summary>

``` r
con <- load_db()
df_sample <- read_table_db(con, "sample") |>
  dplyr::filter(source_id == "geng2024_PRJDB11203")
df_profile <- read_table_db(con, "taxon_profile") |>
  dplyr::filter(sample_id %in% df_sample$sample_id, taxon_kind == "bacterium")
close_db(con)

df_wide <- df_profile |>
  dplyr::select(sample_id, taxon_id, taxon_name, abundance) |>
  tidyr::pivot_wider(id_cols = c(taxon_id, taxon_name),
                     names_from = sample_id, values_from = abundance,
                     values_fill = 0)

cat("Wide matrix:", nrow(df_wide), "species x",
    ncol(df_wide) - 2, "samples\n")
```

</details>

    Wide matrix: 230 species x 42 samples

## Wilcoxon rank-sum (peri vs healthy)

<details class="code-fold">
<summary>Code</summary>

``` r
ids_case <- df_sample$sample_id[df_sample$disease_status == "periodontitis"]
ids_ctrl <- df_sample$sample_id[df_sample$disease_status == "healthy"]
stopifnot(length(ids_case) >= 5, length(ids_ctrl) >= 5)

# species must be present (>0) in >= 10% of all samples to be tested
mat_abund <- as.matrix(df_wide[, c(ids_case, ids_ctrl)])
rownames(mat_abund) <- df_wide$taxon_id
prevalence <- rowMeans(mat_abund > 0)
keep <- prevalence >= 0.10

df_test <- tibble::tibble(
  taxon_id   = df_wide$taxon_id[keep],
  taxon_name = df_wide$taxon_name[keep]
) |>
  dplyr::mutate(
    mean_case = rowMeans(mat_abund[keep, ids_case, drop = FALSE]),
    mean_ctrl = rowMeans(mat_abund[keep, ids_ctrl, drop = FALSE]),
    p_value   = purrr::map2_dbl(
      asplit(mat_abund[keep, ids_case, drop = FALSE], 1),
      asplit(mat_abund[keep, ids_ctrl, drop = FALSE], 1),
      function(a, b) {
        suppressWarnings(stats::wilcox.test(a, b, exact = FALSE)$p.value)
      }
    ),
    effect_size = log2((mean_case + 1e-6) / (mean_ctrl + 1e-6)),
    direction   = ifelse(effect_size > 0, "disease_enriched", "health_enriched"),
    q_value     = stats::p.adjust(p_value, method = "BH")
  )
cat("Tested species:", nrow(df_test), "| q<0.1:", sum(df_test$q_value < 0.1, na.rm = TRUE), "\n")
```

</details>

    Tested species: 136 | q<0.1: 10 

## Build Wilcoxon `candidate_microbe` rows (MVP baseline, preserved)

We keep all tested species with `q_value < 0.1` plus, if fewer than 30
pass, top-up with the smallest-p species so the MVP downstream has
signal. (Ripple will tighten this back to a strict q\<0.1 once
meta-analysis is in place.)

<details class="code-fold">
<summary>Code</summary>

``` r
df_pass <- df_test |>
  dplyr::arrange(p_value) |>
  dplyr::filter(q_value < 0.1)

if (nrow(df_pass) < 30) {
  df_topup <- df_test |>
    dplyr::arrange(p_value) |>
    dplyr::filter(!taxon_id %in% df_pass$taxon_id) |>
    head(30 - nrow(df_pass))
  df_candidates <- dplyr::bind_rows(df_pass, df_topup)
} else {
  df_candidates <- df_pass
}

ver <- schema_version()

df_candidate_microbe <- df_candidates |>
  dplyr::mutate(
    phenotype      = "periodontitis",
    n_cohorts      = 1L,
    method         = "wilcoxon",
    schema_version = ver,
    candidate_id   = paste0(
      "CAND_",
      vapply(paste(taxon_id, phenotype, method, sep = "|"),
             function(s) substr(digest::digest(s, algo = "xxhash64"), 1, 12),
             character(1))
    )
  ) |>
  dplyr::select(candidate_id, taxon_id, taxon_name, phenotype, direction,
                effect_size, p_value, q_value, n_cohorts, method, schema_version)

stopifnot(!any(duplicated(df_candidate_microbe$candidate_id)))
df_candidate_microbe |>
  head(20) |>
  knitr::kable(caption = "Top 20 candidate microbes (smallest p)")
```

</details>

| candidate_id | taxon_id | taxon_name | phenotype | direction | effect_size | p_value | q_value | n_cohorts | method | schema_version |
|:---|:---|:---|:---|:---|---:|---:|---:|---:|:---|:---|
| CAND_6c6398dc311c | STR_Eubacterium_nodatum | Eubacterium_nodatum | periodontitis | disease_enriched | 4.9856648 | 0.0006074 | 0.0535455 | 1 | wilcoxon | 0.1 |
| CAND_87cf61d58fe5 | STR_GGB10852_SGB17523 | GGB10852_SGB17523 | periodontitis | disease_enriched | 4.2609913 | 0.0009587 | 0.0535455 | 1 | wilcoxon | 0.1 |
| CAND_006977768ec7 | STR_Corynebacterium_durum | Corynebacterium_durum | periodontitis | health_enriched | -3.0387145 | 0.0011811 | 0.0535455 | 1 | wilcoxon | 0.1 |
| CAND_68cc8c342e43 | STR_Lautropia_mirabilis | Lautropia_mirabilis | periodontitis | health_enriched | -2.2156737 | 0.0029322 | 0.0848636 | 1 | wilcoxon | 0.1 |
| CAND_5d908f9df765 | STR_Treponema_maltophilum | Treponema_maltophilum | periodontitis | disease_enriched | 3.9140699 | 0.0031200 | 0.0848636 | 1 | wilcoxon | 0.1 |
| CAND_852ca01d05a6 | STR_Porphyromonas_endodontalis | Porphyromonas_endodontalis | periodontitis | disease_enriched | 3.2351803 | 0.0037988 | 0.0861058 | 1 | wilcoxon | 0.1 |
| CAND_4b93285f265a | STR_Filifactor_alocis | Filifactor_alocis | periodontitis | disease_enriched | 4.5283311 | 0.0060398 | 0.0972879 | 1 | wilcoxon | 0.1 |
| CAND_0b5027932aa5 | STR_Fretibacterium_fastidiosum | Fretibacterium_fastidiosum | periodontitis | disease_enriched | 5.0319361 | 0.0060546 | 0.0972879 | 1 | wilcoxon | 0.1 |
| CAND_985e72991e25 | STR_Tannerella_forsythia | Tannerella_forsythia | periodontitis | disease_enriched | 1.7712688 | 0.0067574 | 0.0972879 | 1 | wilcoxon | 0.1 |
| CAND_bb075dbb2347 | STR_Actinomyces_naeslundii | Actinomyces_naeslundii | periodontitis | health_enriched | -1.5256359 | 0.0071535 | 0.0972879 | 1 | wilcoxon | 0.1 |
| CAND_ae12cb2a6107 | STR_Treponema_SGB69443 | Treponema_SGB69443 | periodontitis | disease_enriched | 4.4555455 | 0.0091839 | 0.1083419 | 1 | wilcoxon | 0.1 |
| CAND_c08c8562883e | STR_GGB1025_SGB1319 | GGB1025_SGB1319 | periodontitis | disease_enriched | 3.5913726 | 0.0095596 | 0.1083419 | 1 | wilcoxon | 0.1 |
| CAND_7f28f3c7edb6 | STR_Rothia_aeria | Rothia_aeria | periodontitis | health_enriched | -2.5013694 | 0.0110358 | 0.1154511 | 1 | wilcoxon | 0.1 |
| CAND_91283a05d3f5 | STR_Treponema_denticola | Treponema_denticola | periodontitis | disease_enriched | 3.9003382 | 0.0135171 | 0.1313094 | 1 | wilcoxon | 0.1 |
| CAND_64bf0ec09c7e | STR_Peptostreptococcaceae_bacterium_oral_taxon_113 | Peptostreptococcaceae_bacterium_oral_taxon_113 | periodontitis | disease_enriched | 11.3408312 | 0.0189745 | 0.1624344 | 1 | wilcoxon | 0.1 |
| CAND_0ddb1c305c8b | STR_Actinobaculum_sp_oral_taxon_183 | Actinobaculum_sp_oral_taxon_183 | periodontitis | disease_enriched | 0.3788957 | 0.0191099 | 0.1624344 | 1 | wilcoxon | 0.1 |
| CAND_643a3fc93962 | STR_Mogibacterium_timidum | Mogibacterium_timidum | periodontitis | disease_enriched | 3.0975171 | 0.0226061 | 0.1751742 | 1 | wilcoxon | 0.1 |
| CAND_1bb8facd8dfa | STR_Streptococcus_sanguinis | Streptococcus_sanguinis | periodontitis | health_enriched | -4.0603883 | 0.0231848 | 0.1751742 | 1 | wilcoxon | 0.1 |
| CAND_8aa44ec72dcd | STR_Dialister_invisus | Dialister_invisus | periodontitis | disease_enriched | 1.9967555 | 0.0258216 | 0.1848280 | 1 | wilcoxon | 0.1 |
| CAND_257badcacb3f | STR_GGB4333_SGB5935 | GGB4333_SGB5935 | periodontitis | disease_enriched | 9.1379478 | 0.0324581 | 0.2055235 | 1 | wilcoxon | 0.1 |

Top 20 candidate microbes (smallest p)

<details class="code-fold">
<summary>Code</summary>

``` r
dplyr::count(df_candidate_microbe, direction) |>
  knitr::kable(caption = "Direction balance among Wilcoxon candidates")
```

</details>

| direction        |   n |
|:-----------------|----:|
| disease_enriched |  19 |
| health_enriched  |  11 |

Direction balance among Wilcoxon candidates

## Helper: rows in the long candidate format

<details class="code-fold">
<summary>Code</summary>

``` r
make_candidate_rows <- function(df, method_name) {
  df |>
    dplyr::mutate(
      phenotype      = "periodontitis",
      n_cohorts      = 1L,
      method         = method_name,
      schema_version = ver,
      candidate_id   = paste0(
        "CAND_",
        vapply(paste(taxon_id, phenotype, method, sep = "|"),
               function(s) substr(digest::digest(s, algo = "xxhash64"), 1, 12),
               character(1))
      )
    ) |>
    dplyr::select(candidate_id, taxon_id, taxon_name, phenotype, direction,
                  effect_size, p_value, q_value, n_cohorts, method,
                  schema_version)
}

# Reused inputs for every method
mat_filt   <- mat_abund[keep, , drop = FALSE]                   # filtered abundance
df_features <- df_wide |> dplyr::filter(taxon_id %in% rownames(mat_filt)) |>
  dplyr::select(taxon_id, taxon_name)
df_meta    <- tibble::tibble(
  sample_id      = c(ids_case, ids_ctrl),
  disease_status = c(rep("periodontitis", length(ids_case)),
                     rep("healthy",       length(ids_ctrl)))
) |>
  dplyr::mutate(disease_status = factor(disease_status,
                                        levels = c("healthy", "periodontitis")))
df_meta_df <- as.data.frame(df_meta); rownames(df_meta_df) <- df_meta$sample_id
mat_for_methods <- mat_filt[, df_meta$sample_id, drop = FALSE]
```

</details>

## MaAsLin2

<details class="code-fold">
<summary>Code</summary>

``` r
path_maaslin <- path_target("maaslin2_out")
if (dir.exists(path_maaslin)) unlink(path_maaslin, recursive = TRUE)

# Maaslin2 wants samples as rows, features as cols
mat_maaslin <- t(mat_for_methods)

fit_maaslin <- Maaslin2::Maaslin2(
  input_data        = as.data.frame(mat_maaslin),
  input_metadata    = df_meta_df,
  output            = path_maaslin,
  fixed_effects     = "disease_status",
  reference         = "disease_status,healthy",
  normalization     = "TSS",
  transform         = "LOG",
  analysis_method   = "LM",
  min_prevalence    = 0.10,
  min_abundance     = 0,
  max_significance  = 1,
  plot_heatmap      = FALSE,
  plot_scatter      = FALSE,
  standardize       = FALSE
)
```

</details>

    [1] "Creating output folder"
    [1] "Creating output feature tables folder"
    [1] "Creating output fits folder"
    2026-05-05 21:13:25.649998 INFO::Writing function arguments to log file
    2026-05-05 21:13:25.655989 INFO::Verifying options selected are valid
    2026-05-05 21:13:25.667205 INFO::Determining format of input files
    2026-05-05 21:13:25.667456 INFO::Input format is data samples as rows and metadata samples as rows
    2026-05-05 21:13:25.668539 INFO::Formula for fixed effects: expr ~  disease_status
    2026-05-05 21:13:25.668892 INFO::Filter data based on min abundance and min prevalence
    2026-05-05 21:13:25.669068 INFO::Total samples in data: 42
    2026-05-05 21:13:25.669241 INFO::Min samples required with min abundance for a feature not to be filtered: 4.200000
    2026-05-05 21:13:25.669846 INFO::Total filtered features: 0
    2026-05-05 21:13:25.670059 INFO::Filtered feature names from abundance and prevalence filtering:
    2026-05-05 21:13:25.671058 INFO::Total filtered features with variance filtering: 0
    2026-05-05 21:13:25.671264 INFO::Filtered feature names from variance filtering:
    2026-05-05 21:13:25.671434 INFO::Running selected normalization method: TSS
    2026-05-05 21:13:25.672371 INFO::Bypass z-score application to metadata
    2026-05-05 21:13:25.67258 INFO::Running selected transform method: LOG
    2026-05-05 21:13:25.673424 INFO::Running selected analysis method: LM
    2026-05-05 21:13:25.674913 INFO::Fitting model to feature number 1, STR_Corynebacterium_matruchotii
    2026-05-05 21:13:25.678436 INFO::Fitting model to feature number 2, STR_Capnocytophaga_gingivalis
    2026-05-05 21:13:25.679102 INFO::Fitting model to feature number 3, STR_Capnocytophaga_sputigena
    2026-05-05 21:13:25.679724 INFO::Fitting model to feature number 4, STR_Capnocytophaga_granulosa
    2026-05-05 21:13:25.680338 INFO::Fitting model to feature number 5, STR_Candidatus_Nanosynsacchari_sp_TM7_ANC_38_39_G1_1
    2026-05-05 21:13:25.680932 INFO::Fitting model to feature number 6, STR_Fusobacterium_nucleatum
    2026-05-05 21:13:25.681531 INFO::Fitting model to feature number 7, STR_Lautropia_mirabilis
    2026-05-05 21:13:25.682115 INFO::Fitting model to feature number 8, STR_Aggregatibacter_sp_oral_taxon_458
    2026-05-05 21:13:25.68269 INFO::Fitting model to feature number 9, STR_Tannerella_sp_oral_taxon_HOT_286
    2026-05-05 21:13:25.683296 INFO::Fitting model to feature number 10, STR_Cardiobacterium_hominis
    2026-05-05 21:13:25.683871 INFO::Fitting model to feature number 11, STR_Pseudopropionibacterium_propionicum
    2026-05-05 21:13:25.684467 INFO::Fitting model to feature number 12, STR_Campylobacter_gracilis
    2026-05-05 21:13:25.685043 INFO::Fitting model to feature number 13, STR_Prevotella_conceptionensis
    2026-05-05 21:13:25.68561 INFO::Fitting model to feature number 14, STR_Capnocytophaga_leadbetteri
    2026-05-05 21:13:25.686189 INFO::Fitting model to feature number 15, STR_Actinomyces_naeslundii
    2026-05-05 21:13:25.686742 INFO::Fitting model to feature number 16, STR_Ottowia_sp_Marseille_P4747
    2026-05-05 21:13:25.687325 INFO::Fitting model to feature number 17, STR_Selenomonas_noxia
    2026-05-05 21:13:25.687889 INFO::Fitting model to feature number 18, STR_Selenomonas_artemidis
    2026-05-05 21:13:25.688468 INFO::Fitting model to feature number 19, STR_Neisseria_sp_oral_taxon_014
    2026-05-05 21:13:25.689091 INFO::Fitting model to feature number 20, STR_Actinomyces_massiliensis
    2026-05-05 21:13:25.689674 INFO::Fitting model to feature number 21, STR_Treponema_sp_OMZ_804
    2026-05-05 21:13:25.690269 INFO::Fitting model to feature number 22, STR_Porphyromonas_catoniae
    2026-05-05 21:13:25.690887 INFO::Fitting model to feature number 23, STR_Actinomyces_dentalis
    2026-05-05 21:13:25.691544 INFO::Fitting model to feature number 24, STR_Neisseria_subflava
    2026-05-05 21:13:25.692128 INFO::Fitting model to feature number 25, STR_Actinobaculum_sp_oral_taxon_183
    2026-05-05 21:13:25.692693 INFO::Fitting model to feature number 26, STR_GGB1022_SGB1315
    2026-05-05 21:13:25.693268 INFO::Fitting model to feature number 27, STR_Leptotrichia_hongkongensis
    2026-05-05 21:13:25.693825 INFO::Fitting model to feature number 28, STR_Porphyromonas_pasteri
    2026-05-05 21:13:25.694431 INFO::Fitting model to feature number 29, STR_Prevotella_loescheii
    2026-05-05 21:13:25.695002 INFO::Fitting model to feature number 30, STR_Treponema_socranskii
    2026-05-05 21:13:25.695576 INFO::Fitting model to feature number 31, STR_Capnocytophaga_SGB2480
    2026-05-05 21:13:25.696159 INFO::Fitting model to feature number 32, STR_GGB1843_SGB2524
    2026-05-05 21:13:25.696712 INFO::Fitting model to feature number 33, STR_TM7_phylum_sp_oral_taxon_348
    2026-05-05 21:13:25.697308 INFO::Fitting model to feature number 34, STR_Eikenella_corrodens
    2026-05-05 21:13:25.697872 INFO::Fitting model to feature number 35, STR_Actinomyces_gerencseriae
    2026-05-05 21:13:25.698477 INFO::Fitting model to feature number 36, STR_Alloprevotella_tannerae
    2026-05-05 21:13:25.699066 INFO::Fitting model to feature number 37, STR_Porphyromonas_sp_oral_taxon_278
    2026-05-05 21:13:25.699641 INFO::Fitting model to feature number 38, STR_Prevotella_oris
    2026-05-05 21:13:25.700241 INFO::Fitting model to feature number 39, STR_Actinomyces_johnsonii
    2026-05-05 21:13:25.700811 INFO::Fitting model to feature number 40, STR_Arachnia_SGB15899
    2026-05-05 21:13:25.701396 INFO::Fitting model to feature number 41, STR_Corynebacterium_durum
    2026-05-05 21:13:25.701981 INFO::Fitting model to feature number 42, STR_Actinomyces_oris
    2026-05-05 21:13:25.702573 INFO::Fitting model to feature number 43, STR_Campylobacter_SGB19337
    2026-05-05 21:13:25.70318 INFO::Fitting model to feature number 44, STR_Campylobacter_SGB19317
    2026-05-05 21:13:25.703757 INFO::Fitting model to feature number 45, STR_Streptococcus_oralis
    2026-05-05 21:13:25.704356 INFO::Fitting model to feature number 46, STR_Streptococcus_cristatus
    2026-05-05 21:13:25.704927 INFO::Fitting model to feature number 47, STR_Prevotella_oulorum
    2026-05-05 21:13:25.705521 INFO::Fitting model to feature number 48, STR_Rothia_dentocariosa
    2026-05-05 21:13:25.706097 INFO::Fitting model to feature number 49, STR_Leptotrichia_buccalis
    2026-05-05 21:13:25.706692 INFO::Fitting model to feature number 50, STR_Selenomonas_SGB5891
    2026-05-05 21:13:25.707265 INFO::Fitting model to feature number 51, STR_Candidatus_Saccharibacteria_bacterium_oral_taxon_488
    2026-05-05 21:13:25.707857 INFO::Fitting model to feature number 52, STR_Rothia_aeria
    2026-05-05 21:13:25.708433 INFO::Fitting model to feature number 53, STR_Prevotella_maculosa
    2026-05-05 21:13:25.709045 INFO::Fitting model to feature number 54, STR_Veillonella_parvula
    2026-05-05 21:13:25.709635 INFO::Fitting model to feature number 55, STR_Prevotella_nigrescens
    2026-05-05 21:13:25.710254 INFO::Fitting model to feature number 56, STR_Treponema_maltophilum
    2026-05-05 21:13:25.710837 INFO::Fitting model to feature number 57, STR_Prevotella_saccharolytica
    2026-05-05 21:13:25.711429 INFO::Fitting model to feature number 58, STR_GGB4400_SGB6074
    2026-05-05 21:13:25.712052 INFO::Fitting model to feature number 59, STR_Actinomyces_israelii
    2026-05-05 21:13:25.712642 INFO::Fitting model to feature number 60, STR_Prevotella_oralis
    2026-05-05 21:13:25.713227 INFO::Fitting model to feature number 61, STR_Cardiobacterium_valvarum
    2026-05-05 21:13:25.713834 INFO::Fitting model to feature number 62, STR_Lachnospiraceae_bacterium_oral_taxon_500
    2026-05-05 21:13:25.714488 INFO::Fitting model to feature number 63, STR_Leptotrichia_sp_oral_taxon_212
    2026-05-05 21:13:25.715094 INFO::Fitting model to feature number 64, STR_Streptococcus_sanguinis
    2026-05-05 21:13:25.7157 INFO::Fitting model to feature number 65, STR_Actinomyces_sp_oral_taxon_448
    2026-05-05 21:13:25.716278 INFO::Fitting model to feature number 66, STR_Bacteroidetes_oral_taxon_274
    2026-05-05 21:13:25.716866 INFO::Fitting model to feature number 67, STR_Porphyromonas_gingivalis
    2026-05-05 21:13:25.717438 INFO::Fitting model to feature number 68, STR_Campylobacter_curvus
    2026-05-05 21:13:25.718029 INFO::Fitting model to feature number 69, STR_Neisseria_sicca
    2026-05-05 21:13:25.718614 INFO::Fitting model to feature number 70, STR_GGB10022_SGB15896
    2026-05-05 21:13:25.719179 INFO::Fitting model to feature number 71, STR_Haemophilus_parainfluenzae
    2026-05-05 21:13:25.719783 INFO::Fitting model to feature number 72, STR_Neisseria_elongata
    2026-05-05 21:13:25.720347 INFO::Fitting model to feature number 73, STR_Streptococcus_gordonii
    2026-05-05 21:13:25.720942 INFO::Fitting model to feature number 74, STR_Capnocytophaga_ochracea
    2026-05-05 21:13:25.721517 INFO::Fitting model to feature number 75, STR_Leptotrichia_hofstadii
    2026-05-05 21:13:25.722107 INFO::Fitting model to feature number 76, STR_Kytococcus_sedentarius
    2026-05-05 21:13:25.72271 INFO::Fitting model to feature number 77, STR_Prevotella_veroralis
    2026-05-05 21:13:25.72329 INFO::Fitting model to feature number 78, STR_Pauljensenia_hongkongensis
    2026-05-05 21:13:25.7239 INFO::Fitting model to feature number 79, STR_Isoptericola_variabilis
    2026-05-05 21:13:25.7245 INFO::Fitting model to feature number 80, STR_Granulicatella_adiacens
    2026-05-05 21:13:25.725112 INFO::Fitting model to feature number 81, STR_candidate_division_SR1_bacterium_MGEHA
    2026-05-05 21:13:25.725757 INFO::Fitting model to feature number 82, STR_Selenomonas_sputigena
    2026-05-05 21:13:25.726375 INFO::Fitting model to feature number 83, STR_GGB1201_SGB1566
    2026-05-05 21:13:25.727013 INFO::Fitting model to feature number 84, STR_Anaeroglobus_geminatus
    2026-05-05 21:13:25.727666 INFO::Fitting model to feature number 85, STR_GGB6675_SGB9425
    2026-05-05 21:13:25.72828 INFO::Fitting model to feature number 86, STR_Arachnia_SGB15898
    2026-05-05 21:13:25.728917 INFO::Fitting model to feature number 87, STR_Peptidiphaga_gingivicola
    2026-05-05 21:13:25.729531 INFO::Fitting model to feature number 88, STR_Lautropia_dentalis
    2026-05-05 21:13:25.730149 INFO::Fitting model to feature number 89, STR_Porphyromonas_endodontalis
    2026-05-05 21:13:25.730794 INFO::Fitting model to feature number 90, STR_GGB12796_SGB19893
    2026-05-05 21:13:25.731396 INFO::Fitting model to feature number 91, STR_Dialister_invisus
    2026-05-05 21:13:25.732029 INFO::Fitting model to feature number 92, STR_GGB4533_SGB6246
    2026-05-05 21:13:25.732655 INFO::Fitting model to feature number 93, STR_GGB12790_SGB19844
    2026-05-05 21:13:25.733282 INFO::Fitting model to feature number 94, STR_Olsenella_sp_oral_taxon_807
    2026-05-05 21:13:25.733907 INFO::Fitting model to feature number 95, STR_Tannerella_forsythia
    2026-05-05 21:13:25.734524 INFO::Fitting model to feature number 96, STR_Campylobacter_rectus
    2026-05-05 21:13:25.735136 INFO::Fitting model to feature number 97, STR_Eubacterium_brachy
    2026-05-05 21:13:25.735774 INFO::Fitting model to feature number 98, STR_Treponema_denticola
    2026-05-05 21:13:25.736377 INFO::Fitting model to feature number 99, STR_GGB12790_SGB19845
    2026-05-05 21:13:25.737013 INFO::Fitting model to feature number 100, STR_Lachnoanaerobaculum_saburreum
    2026-05-05 21:13:25.737639 INFO::Fitting model to feature number 101, STR_Actinomyces_sp_oral_taxon_897
    2026-05-05 21:13:25.738239 INFO::Fitting model to feature number 102, STR_Prevotella_denticola
    2026-05-05 21:13:25.738878 INFO::Fitting model to feature number 103, STR_Candidatus_Nanoperiomorbus_periodonticus
    2026-05-05 21:13:25.739507 INFO::Fitting model to feature number 104, STR_Prevotella_intermedia
    2026-05-05 21:13:25.740194 INFO::Fitting model to feature number 105, STR_GGB3385_SGB4472
    2026-05-05 21:13:25.740804 INFO::Fitting model to feature number 106, STR_Selenomonas_massiliensis
    2026-05-05 21:13:25.741443 INFO::Fitting model to feature number 107, STR_Centipeda_periodontii
    2026-05-05 21:13:25.742069 INFO::Fitting model to feature number 108, STR_Mogibacterium_timidum
    2026-05-05 21:13:25.742712 INFO::Fitting model to feature number 109, STR_Gemella_morbillorum
    2026-05-05 21:13:25.743351 INFO::Fitting model to feature number 110, STR_Kingella_denitrificans
    2026-05-05 21:13:25.744012 INFO::Fitting model to feature number 111, STR_Candidatus_Gracilibacteria_bacterium_GN02_873
    2026-05-05 21:13:25.744657 INFO::Fitting model to feature number 112, STR_Selenomonas_sp_oral_taxon_892
    2026-05-05 21:13:25.7453 INFO::Fitting model to feature number 113, STR_Selenomonas_sp_oral_taxon_126
    2026-05-05 21:13:25.745958 INFO::Fitting model to feature number 114, STR_Catonella_morbi
    2026-05-05 21:13:25.746595 INFO::Fitting model to feature number 115, STR_Parvimonas_micra
    2026-05-05 21:13:25.747233 INFO::Fitting model to feature number 116, STR_Bacteroidetes_unclassified_SGB1343
    2026-05-05 21:13:25.747882 INFO::Fitting model to feature number 117, STR_Campylobacter_showae
    2026-05-05 21:13:25.748536 INFO::Fitting model to feature number 118, STR_Dialister_pneumosintes
    2026-05-05 21:13:25.749166 INFO::Fitting model to feature number 119, STR_Prevotella_baroniae
    2026-05-05 21:13:25.749821 INFO::Fitting model to feature number 120, STR_Streptococcus_anginosus
    2026-05-05 21:13:25.750444 INFO::Fitting model to feature number 121, STR_Peptostreptococcus_stomatis
    2026-05-05 21:13:25.75109 INFO::Fitting model to feature number 122, STR_GGB1025_SGB1319
    2026-05-05 21:13:25.751751 INFO::Fitting model to feature number 123, STR_GGB10852_SGB17523
    2026-05-05 21:13:25.752373 INFO::Fitting model to feature number 124, STR_Filifactor_alocis
    2026-05-05 21:13:25.753028 INFO::Fitting model to feature number 125, STR_Eubacterium_yurii
    2026-05-05 21:13:25.75366 INFO::Fitting model to feature number 126, STR_Prevotella_sp_oral_taxon_473
    2026-05-05 21:13:25.75432 INFO::Fitting model to feature number 127, STR_Treponema_SGB69443
    2026-05-05 21:13:25.754954 INFO::Fitting model to feature number 128, STR_GGB49229_SGB69060
    2026-05-05 21:13:25.755601 INFO::Fitting model to feature number 129, STR_Treponema_sp_Marseille_Q3903
    2026-05-05 21:13:25.756255 INFO::Fitting model to feature number 130, STR_GGB1611_SGB2208
    2026-05-05 21:13:25.756886 INFO::Fitting model to feature number 131, STR_Eubacterium_nodatum
    2026-05-05 21:13:25.757537 INFO::Fitting model to feature number 132, STR_Fretibacterium_fastidiosum
    2026-05-05 21:13:25.758186 INFO::Fitting model to feature number 133, STR_Prevotella_dentalis
    2026-05-05 21:13:25.75881 INFO::Fitting model to feature number 134, STR_Anaerolineaceae_bacterium_oral_taxon_439
    2026-05-05 21:13:25.759436 INFO::Fitting model to feature number 135, STR_GGB4333_SGB5935
    2026-05-05 21:13:25.760056 INFO::Fitting model to feature number 136, STR_Peptostreptococcaceae_bacterium_oral_taxon_113
    2026-05-05 21:13:25.765822 INFO::Counting total values for each feature
    2026-05-05 21:13:25.767678 INFO::Writing filtered data to file /Volumes/ssd2t/github/rujinlong/pr0003-OralProbPhage/analyses/data/050-m1a-differential-analysis/maaslin2_out/features/filtered_data.tsv
    2026-05-05 21:13:25.770061 INFO::Writing filtered, normalized data to file /Volumes/ssd2t/github/rujinlong/pr0003-OralProbPhage/analyses/data/050-m1a-differential-analysis/maaslin2_out/features/filtered_data_norm.tsv
    2026-05-05 21:13:25.77254 INFO::Writing filtered, normalized, transformed data to file /Volumes/ssd2t/github/rujinlong/pr0003-OralProbPhage/analyses/data/050-m1a-differential-analysis/maaslin2_out/features/filtered_data_norm_transformed.tsv
    2026-05-05 21:13:25.775571 INFO::Writing residuals to file /Volumes/ssd2t/github/rujinlong/pr0003-OralProbPhage/analyses/data/050-m1a-differential-analysis/maaslin2_out/fits/residuals.rds
    2026-05-05 21:13:25.776387 INFO::Writing fitted values to file /Volumes/ssd2t/github/rujinlong/pr0003-OralProbPhage/analyses/data/050-m1a-differential-analysis/maaslin2_out/fits/fitted.rds
    2026-05-05 21:13:25.77692 INFO::Writing all results to file (ordered by increasing q-values): /Volumes/ssd2t/github/rujinlong/pr0003-OralProbPhage/analyses/data/050-m1a-differential-analysis/maaslin2_out/all_results.tsv
    2026-05-05 21:13:25.777687 INFO::Writing the significant results (those which are less than or equal to the threshold of 1.000000 ) to file (ordered by increasing q-values): /Volumes/ssd2t/github/rujinlong/pr0003-OralProbPhage/analyses/data/050-m1a-differential-analysis/maaslin2_out/significant_results.tsv

<details class="code-fold">
<summary>Code</summary>

``` r
df_maaslin <- tibble::as_tibble(fit_maaslin$results) |>
  dplyr::filter(metadata == "disease_status",
                value    == "periodontitis") |>
  dplyr::transmute(
    taxon_id    = feature,
    effect_size = coef,
    p_value     = pval,
    q_value     = qval,
    direction   = ifelse(coef > 0, "disease_enriched", "health_enriched")
  ) |>
  dplyr::left_join(df_features, by = "taxon_id") |>
  dplyr::filter(!is.na(p_value))

cat("MaAsLin2: tested", nrow(df_maaslin),
    "| q<0.1:", sum(df_maaslin$q_value < 0.1, na.rm = TRUE), "\n")
```

</details>

    MaAsLin2: tested 136 | q<0.1: 13 

<details class="code-fold">
<summary>Code</summary>

``` r
df_cand_maaslin <- df_maaslin |>
  dplyr::filter(q_value < 0.1) |>
  make_candidate_rows("maaslin2")
```

</details>

## ANCOMBC

<details class="code-fold">
<summary>Code</summary>

``` r
# ANCOMBC needs counts; round(rel_abund * 1e6) is the standard trick.
mat_counts <- round(mat_for_methods * 1e6)
storage.mode(mat_counts) <- "integer"

ps_otu  <- phyloseq::otu_table(mat_counts, taxa_are_rows = TRUE)
ps_samp <- phyloseq::sample_data(df_meta_df)
ps_obj  <- phyloseq::phyloseq(ps_otu, ps_samp)

fit_ancom <- ANCOMBC::ancombc(
  data         = ps_obj,
  formula      = "disease_status",
  p_adj_method = "BH",
  prv_cut      = 0.10,
  group        = "disease_status",
  struc_zero   = FALSE,
  neg_lb       = FALSE,
  conserve     = TRUE,
  verbose      = FALSE
)

df_ancom_raw <- tibble::tibble(
  taxon_id    = fit_ancom$res$lfc$taxon,
  effect_size = fit_ancom$res$lfc[["disease_statusperiodontitis"]],
  p_value     = fit_ancom$res$p_val[["disease_statusperiodontitis"]],
  q_value     = fit_ancom$res$q_val[["disease_statusperiodontitis"]]
) |>
  dplyr::mutate(
    direction = ifelse(effect_size > 0, "disease_enriched", "health_enriched")
  ) |>
  dplyr::left_join(df_features, by = "taxon_id") |>
  dplyr::filter(!is.na(p_value))

cat("ANCOMBC: tested", nrow(df_ancom_raw),
    "| q<0.1:", sum(df_ancom_raw$q_value < 0.1, na.rm = TRUE), "\n")
```

</details>

    ANCOMBC: tested 136 | q<0.1: 14 

<details class="code-fold">
<summary>Code</summary>

``` r
df_cand_ancom <- df_ancom_raw |>
  dplyr::filter(q_value < 0.1) |>
  make_candidate_rows("ancombc")
```

</details>

## DESeq2

<details class="code-fold">
<summary>Code</summary>

``` r
# Same count trick. DESeq2 requires integer counts.
dds <- DESeq2::DESeqDataSetFromMatrix(
  countData = mat_counts,
  colData   = df_meta_df,
  design    = ~ disease_status
)
# poscounts size factors handle zero-rich microbiome data better than RLE
dds <- DESeq2::estimateSizeFactors(dds, type = "poscounts")
dds <- DESeq2::DESeq(dds, quiet = TRUE)
res <- DESeq2::results(dds,
                      contrast = c("disease_status", "periodontitis", "healthy"))
df_deseq <- tibble::tibble(
  taxon_id    = rownames(res),
  effect_size = res$log2FoldChange,
  p_value     = res$pvalue,
  q_value     = stats::p.adjust(res$pvalue, method = "BH")
) |>
  dplyr::mutate(
    direction = ifelse(effect_size > 0, "disease_enriched", "health_enriched")
  ) |>
  dplyr::left_join(df_features, by = "taxon_id") |>
  dplyr::filter(!is.na(p_value))

cat("DESeq2: tested", nrow(df_deseq),
    "| q<0.1:", sum(df_deseq$q_value < 0.1, na.rm = TRUE), "\n")
```

</details>

    DESeq2: tested 136 | q<0.1: 27 

<details class="code-fold">
<summary>Code</summary>

``` r
df_cand_deseq <- df_deseq |>
  dplyr::filter(q_value < 0.1) |>
  make_candidate_rows("deseq2")
```

</details>

## Wilcoxon → method-tagged rows (full results)

<details class="code-fold">
<summary>Code</summary>

``` r
df_wilcox_full <- df_test |>
  dplyr::transmute(taxon_id, taxon_name, effect_size, p_value, q_value,
                   direction)
```

</details>

## Method agreement diagnostics

<details class="code-fold">
<summary>Code</summary>

``` r
df_methods <- list(
  wilcoxon = df_wilcox_full,
  maaslin2 = df_maaslin,
  ancombc  = df_ancom_raw,
  deseq2   = df_deseq
)

# significant set per method (q<0.1, with direction)
sig_sets <- purrr::imap(df_methods, function(df, m) {
  df |>
    dplyr::filter(q_value < 0.1) |>
    dplyr::transmute(taxon_id, direction)
})

df_agree <- tidyr::expand_grid(a = names(sig_sets), b = names(sig_sets)) |>
  dplyr::filter(a < b) |>
  dplyr::mutate(
    n_a       = vapply(a, function(x) nrow(sig_sets[[x]]), integer(1)),
    n_b       = vapply(b, function(x) nrow(sig_sets[[x]]), integer(1)),
    n_overlap = purrr::map2_int(a, b, function(x, y) {
      length(intersect(sig_sets[[x]]$taxon_id, sig_sets[[y]]$taxon_id))
    }),
    jaccard = purrr::map2_dbl(a, b, function(x, y) {
      ix <- intersect(sig_sets[[x]]$taxon_id, sig_sets[[y]]$taxon_id)
      un <- union(sig_sets[[x]]$taxon_id, sig_sets[[y]]$taxon_id)
      if (length(un) == 0) NA_real_ else length(ix) / length(un)
    })
  )
knitr::kable(df_agree, digits = 3,
             caption = "Pairwise overlap among q<0.1 sets across methods")
```

</details>

| a        | b        | n_a | n_b | n_overlap | jaccard |
|:---------|:---------|----:|----:|----------:|--------:|
| maaslin2 | wilcoxon |  13 |  10 |        10 |   0.769 |
| ancombc  | wilcoxon |  14 |  10 |        10 |   0.714 |
| ancombc  | maaslin2 |  14 |  13 |        10 |   0.588 |
| ancombc  | deseq2   |  14 |  27 |         2 |   0.051 |
| deseq2   | wilcoxon |  27 |  10 |         1 |   0.028 |
| deseq2   | maaslin2 |  27 |  13 |         2 |   0.053 |

Pairwise overlap among q\<0.1 sets across methods

## Consensus candidate set

A species enters the consensus set if it is significant (q\<0.1) in ≥ 2
of the 4 methods AND all flagging methods agree on direction. The
consensus row’s effect-size is the mean of the flagging methods’
effect-sizes, and the q-value is the **maximum** across flagging methods
(conservative).

<details class="code-fold">
<summary>Code</summary>

``` r
df_long_sig <- purrr::imap_dfr(df_methods, function(df, m) {
  df |>
    dplyr::filter(q_value < 0.1) |>
    dplyr::mutate(method = m)
})

df_consensus <- df_long_sig |>
  dplyr::group_by(taxon_id, taxon_name) |>
  dplyr::summarise(
    n_methods         = dplyr::n(),
    n_dir_disease     = sum(direction == "disease_enriched"),
    n_dir_health      = sum(direction == "health_enriched"),
    effect_size       = mean(effect_size, na.rm = TRUE),
    q_value           = max(q_value,      na.rm = TRUE),
    methods_flagging  = paste(sort(method), collapse = ","),
    .groups = "drop"
  ) |>
  dplyr::filter(n_methods >= 2,
                n_dir_disease == 0 | n_dir_health == 0) |>
  dplyr::mutate(
    direction = ifelse(n_dir_disease > 0, "disease_enriched", "health_enriched"),
    p_value   = NA_real_
  )

cat("Consensus species (q<0.1 in >=2 methods + concordant direction):",
    nrow(df_consensus), "\n")
```

</details>

    Consensus species (q<0.1 in >=2 methods + concordant direction): 12 

<details class="code-fold">
<summary>Code</summary>

``` r
df_cand_consensus <- df_consensus |>
  dplyr::select(taxon_id, taxon_name, effect_size, p_value, q_value,
                direction) |>
  make_candidate_rows("consensus")

if (nrow(df_consensus) > 0) {
  df_consensus |>
    dplyr::arrange(q_value) |>
    head(15) |>
    dplyr::select(taxon_id, taxon_name, n_methods, direction,
                  effect_size, q_value, methods_flagging) |>
    knitr::kable(digits = 3,
                 caption = "Top 15 consensus species (sorted by q_value)")
}
```

</details>

| taxon_id | taxon_name | n_methods | direction | effect_size | q_value | methods_flagging |
|:---|:---|---:|:---|---:|---:|:---|
| STR_Eubacterium_nodatum | Eubacterium_nodatum | 3 | disease_enriched | 3.894 | 0.054 | ancombc,maaslin2,wilcoxon |
| STR_GGB10852_SGB17523 | GGB10852_SGB17523 | 3 | disease_enriched | 3.369 | 0.054 | ancombc,maaslin2,wilcoxon |
| STR_Corynebacterium_durum | Corynebacterium_durum | 3 | health_enriched | -3.211 | 0.055 | ancombc,maaslin2,wilcoxon |
| STR_GGB4333_SGB5935 | GGB4333_SGB5935 | 2 | disease_enriched | 16.477 | 0.066 | ancombc,deseq2 |
| STR_Lautropia_mirabilis | Lautropia_mirabilis | 3 | health_enriched | -3.016 | 0.085 | ancombc,maaslin2,wilcoxon |
| STR_Treponema_maltophilum | Treponema_maltophilum | 3 | disease_enriched | 3.396 | 0.085 | ancombc,maaslin2,wilcoxon |
| STR_Porphyromonas_endodontalis | Porphyromonas_endodontalis | 3 | disease_enriched | 3.691 | 0.086 | ancombc,maaslin2,wilcoxon |
| STR_Actinomyces_massiliensis | Actinomyces_massiliensis | 2 | health_enriched | -3.547 | 0.095 | deseq2,maaslin2 |
| STR_Actinomyces_naeslundii | Actinomyces_naeslundii | 3 | health_enriched | -2.685 | 0.097 | ancombc,maaslin2,wilcoxon |
| STR_Filifactor_alocis | Filifactor_alocis | 3 | disease_enriched | 3.620 | 0.097 | ancombc,maaslin2,wilcoxon |
| STR_Fretibacterium_fastidiosum | Fretibacterium_fastidiosum | 4 | disease_enriched | 4.763 | 0.097 | ancombc,deseq2,maaslin2,wilcoxon |
| STR_Tannerella_forsythia | Tannerella_forsythia | 3 | disease_enriched | 3.095 | 0.097 | ancombc,maaslin2,wilcoxon |

Top 15 consensus species (sorted by q_value)

## Combine all methods + write to the contract layer

<details class="code-fold">
<summary>Code</summary>

``` r
df_all_methods <- dplyr::bind_rows(
  df_candidate_microbe,    # wilcoxon (preserved MVP block)
  df_cand_maaslin,
  df_cand_ancom,
  df_cand_deseq,
  df_cand_consensus
)

stopifnot(!any(duplicated(df_all_methods$candidate_id)))

con <- load_db()
assert_schema(con)
DBI::dbExecute(con, "DELETE FROM candidate_microbe WHERE phenotype = 'periodontitis'")
```

</details>

    [1] 96

<details class="code-fold">
<summary>Code</summary>

``` r
DBI::dbWriteTable(con, "candidate_microbe", df_all_methods,
                  append = TRUE, row.names = FALSE)
df_summary <- DBI::dbGetQuery(
  con,
  "SELECT method, direction, COUNT(*) AS n
     FROM candidate_microbe
    WHERE phenotype = 'periodontitis'
 GROUP BY method, direction
 ORDER BY method, direction"
)
n_total <- DBI::dbGetQuery(con,
  "SELECT COUNT(*) AS n FROM candidate_microbe")$n
close_db(con)

knitr::kable(df_summary,
             caption = "candidate_microbe rows by method × direction (in DB)")
```

</details>

| method    | direction        |   n |
|:----------|:-----------------|----:|
| ancombc   | disease_enriched |  10 |
| ancombc   | health_enriched  |   4 |
| consensus | disease_enriched |   8 |
| consensus | health_enriched  |   4 |
| deseq2    | disease_enriched |  21 |
| deseq2    | health_enriched  |   6 |
| maaslin2  | disease_enriched |   8 |
| maaslin2  | health_enriched  |   5 |
| wilcoxon  | disease_enriched |  19 |
| wilcoxon  | health_enriched  |  11 |

candidate_microbe rows by method × direction (in DB)

<details class="code-fold">
<summary>Code</summary>

``` r
cat("Total candidate_microbe rows in DB:", n_total, "\n")
```

</details>

    Total candidate_microbe rows in DB: 96 

## Persist

<details class="code-fold">
<summary>Code</summary>

``` r
write_csv(df_test,           path_target("wilcoxon_full_results.csv"))
write_csv(df_maaslin,        path_target("maaslin2_full_results.csv"))
write_csv(df_ancom_raw,      path_target("ancombc_full_results.csv"))
write_csv(df_deseq,          path_target("deseq2_full_results.csv"))
write_csv(df_all_methods,    path_target("candidate_microbe.csv"))
write_csv(df_agree,          path_target("method_agreement.csv"))
write_csv(df_consensus,      path_target("consensus_species.csv"))
```

</details>

## Cross-cohort meta-analysis (050-r2)

Combine per-cohort Wilcoxon statistics across the six Geng 2024
sub-cohorts via Stouffer’s Z (sqrt-n weights, sign tracking on log2FC)
and a DerSimonian–Laird random-effects pool of rank-biserial
correlations. Outputs two new `candidate_microbe` method tags:
`meta_stouffer` and `meta_random`, with `n_cohorts >= 2`.

<details class="code-fold">
<summary>Code</summary>

``` r
con <- load_db()
df_sample_all <- read_table_db(con, "sample") |>
  dplyr::filter(grepl("^geng2024_", source_id),
                disease_status %in% c("periodontitis", "healthy"))
df_profile_all <- read_table_db(con, "taxon_profile") |>
  dplyr::filter(sample_id %in% df_sample_all$sample_id,
                taxon_kind == "bacterium")
close_db(con)

cohorts_meta <- sort(unique(df_sample_all$source_id))
cat("Meta-analysis cohorts:", length(cohorts_meta), "\n")
```

</details>

    Meta-analysis cohorts: 6 

<details class="code-fold">
<summary>Code</summary>

``` r
per_cohort_wilcox <- function(source_id_x) {
  df_s <- df_sample_all |> dplyr::filter(source_id == source_id_x)
  df_p <- df_profile_all |> dplyr::filter(sample_id %in% df_s$sample_id)
  df_w <- df_p |>
    dplyr::select(sample_id, taxon_id, taxon_name, abundance) |>
    tidyr::pivot_wider(id_cols = c(taxon_id, taxon_name),
                       names_from = sample_id, values_from = abundance,
                       values_fill = 0)
  ids_case_x <- df_s$sample_id[df_s$disease_status == "periodontitis"]
  ids_ctrl_x <- df_s$sample_id[df_s$disease_status == "healthy"]
  if (length(ids_case_x) < 5 || length(ids_ctrl_x) < 5) return(NULL)
  mat_x <- as.matrix(df_w[, c(ids_case_x, ids_ctrl_x)])
  rownames(mat_x) <- df_w$taxon_id
  prev_x <- rowMeans(mat_x > 0)
  keep_x <- prev_x >= 0.10
  if (sum(keep_x) == 0) return(NULL)
  mat_case <- mat_x[keep_x, ids_case_x, drop = FALSE]
  mat_ctrl <- mat_x[keep_x, ids_ctrl_x, drop = FALSE]
  n1 <- length(ids_case_x); n2 <- length(ids_ctrl_x); n_total <- n1 + n2
  taxa_x <- rownames(mat_case)
  out <- purrr::map_dfr(seq_along(taxa_x), function(i) {
    a <- mat_case[i, ]; b <- mat_ctrl[i, ]
    if (length(unique(c(a, b))) < 2) {
      return(tibble::tibble(taxon_id = taxa_x[i],
                            p_value = NA_real_, U = NA_real_,
                            r = NA_real_, var_r = NA_real_,
                            log2fc = 0, mean_case = mean(a),
                            mean_ctrl = mean(b)))
    }
    wt <- suppressWarnings(stats::wilcox.test(a, b, exact = FALSE))
    U <- unname(wt$statistic)
    r <- 2 * U / (n1 * n2) - 1        # positive r = case-enriched
    var_r <- (1 - r^2)^2 / max(n_total - 1, 1)
    if (var_r < 1e-8) var_r <- 1e-8
    tibble::tibble(taxon_id = taxa_x[i],
                   p_value = wt$p.value, U = U,
                   r = r, var_r = var_r,
                   log2fc = log2((mean(a) + 1e-6) / (mean(b) + 1e-6)),
                   mean_case = mean(a), mean_ctrl = mean(b))
  })
  out |>
    dplyr::left_join(df_w |> dplyr::select(taxon_id, taxon_name),
                     by = "taxon_id") |>
    dplyr::mutate(source_id = source_id_x, n_case = n1, n_ctrl = n2,
                  n_total = n_total) |>
    dplyr::filter(!is.na(p_value))
}

df_per_cohort <- purrr::map_dfr(cohorts_meta, per_cohort_wilcox)

df_per_cohort |>
  dplyr::count(source_id, name = "n_species_tested") |>
  knitr::kable(caption = "Per-cohort Wilcoxon: species tested (prevalence >= 10%)")
```

</details>

| source_id            | n_species_tested |
|:---------------------|-----------------:|
| geng2024_PRJDB11203  |              136 |
| geng2024_PRJNA230363 |              280 |
| geng2024_PRJNA396840 |              210 |
| geng2024_PRJNA678453 |              293 |
| geng2024_PRJNA717815 |              283 |
| geng2024_PRJNA932553 |              333 |

Per-cohort Wilcoxon: species tested (prevalence \>= 10%)

<details class="code-fold">
<summary>Code</summary>

``` r
df_per_cohort_clean <- df_per_cohort |>
  dplyr::filter(is.finite(p_value), p_value > 0, p_value < 1,
                is.finite(r), is.finite(var_r))

df_grouped <- df_per_cohort_clean |>
  dplyr::group_by(taxon_id) |>
  dplyr::summarise(
    n_cohorts  = dplyr::n_distinct(source_id),
    taxon_name = dplyr::first(stats::na.omit(taxon_name)),
    .groups = "drop"
  ) |>
  dplyr::filter(n_cohorts >= 2)

cat("Species testable in >= 2 cohorts:", nrow(df_grouped), "\n")
```

</details>

    Species testable in >= 2 cohorts: 334 

<details class="code-fold">
<summary>Code</summary>

``` r
stouffer_one <- function(df_one) {
  pv  <- df_one$p_value
  fc  <- df_one$log2fc
  n   <- df_one$n_total
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

df_meta_stouffer <- df_per_cohort_clean |>
  dplyr::semi_join(df_grouped, by = "taxon_id") |>
  dplyr::group_by(taxon_id) |>
  dplyr::group_modify(~ stouffer_one(.x)) |>
  dplyr::ungroup() |>
  dplyr::left_join(df_grouped, by = "taxon_id") |>
  dplyr::mutate(q_value = stats::p.adjust(p_value, method = "BH"))

cat("Stouffer meta: ", nrow(df_meta_stouffer),
    " species | q<0.1: ", sum(df_meta_stouffer$q_value < 0.1, na.rm = TRUE),
    "\n", sep = "")
```

</details>

    Stouffer meta: 334 species | q<0.1: 141

<details class="code-fold">
<summary>Code</summary>

``` r
has_metafor <- requireNamespace("metafor", quietly = TRUE)
random_one <- function(df_one) {
  yi <- df_one$r; vi <- df_one$var_r
  if (has_metafor) {
    fit <- try(metafor::rma(yi = yi, vi = vi, method = "DL"), silent = TRUE)
    if (!inherits(fit, "try-error")) {
      return(tibble::tibble(effect_size = unname(fit$b[1, 1]),
                            se          = fit$se,
                            p_value     = fit$pval,
                            tau2        = fit$tau2,
                            I2          = fit$I2))
    }
  }
  # DL fallback: IVW
  w <- 1 / vi
  est <- sum(w * yi) / sum(w)
  se  <- sqrt(1 / sum(w))
  z   <- est / se
  tibble::tibble(effect_size = est, se = se,
                 p_value     = 2 * stats::pnorm(-abs(z)),
                 tau2 = NA_real_, I2 = NA_real_)
}

df_meta_random <- df_per_cohort_clean |>
  dplyr::semi_join(df_grouped, by = "taxon_id") |>
  dplyr::group_by(taxon_id) |>
  dplyr::group_modify(~ random_one(.x)) |>
  dplyr::ungroup() |>
  dplyr::left_join(df_grouped, by = "taxon_id") |>
  dplyr::mutate(q_value   = stats::p.adjust(p_value, method = "BH"),
                direction = ifelse(effect_size >= 0, "disease_enriched",
                                                     "health_enriched"))

cat("Random-effects meta (",
    ifelse(has_metafor, "DL via metafor", "DL fallback: IVW"),
    "): ", nrow(df_meta_random),
    " species | q<0.1: ", sum(df_meta_random$q_value < 0.1, na.rm = TRUE),
    "\n", sep = "")
```

</details>

    Random-effects meta (DL via metafor): 334 species | q<0.1: 119

<details class="code-fold">
<summary>Code</summary>

``` r
df_cand_meta_stouffer <- df_meta_stouffer |>
  dplyr::transmute(taxon_id, taxon_name, direction,
                   effect_size, p_value, q_value, n_cohorts) |>
  dplyr::mutate(
    phenotype      = "periodontitis",
    method         = "meta_stouffer",
    schema_version = ver,
    candidate_id   = paste0(
      "CAND_",
      vapply(paste(taxon_id, phenotype, method, sep = "|"),
             function(s) substr(digest::digest(s, algo = "xxhash64"), 1, 12),
             character(1))
    )
  ) |>
  dplyr::select(candidate_id, taxon_id, taxon_name, phenotype, direction,
                effect_size, p_value, q_value, n_cohorts, method,
                schema_version)

df_cand_meta_random <- df_meta_random |>
  dplyr::transmute(taxon_id, taxon_name, direction,
                   effect_size, p_value, q_value, n_cohorts) |>
  dplyr::mutate(
    phenotype      = "periodontitis",
    method         = "meta_random",
    schema_version = ver,
    candidate_id   = paste0(
      "CAND_",
      vapply(paste(taxon_id, phenotype, method, sep = "|"),
             function(s) substr(digest::digest(s, algo = "xxhash64"), 1, 12),
             character(1))
    )
  ) |>
  dplyr::select(candidate_id, taxon_id, taxon_name, phenotype, direction,
                effect_size, p_value, q_value, n_cohorts, method,
                schema_version)

df_cand_meta <- dplyr::bind_rows(df_cand_meta_stouffer, df_cand_meta_random)
stopifnot(!any(duplicated(df_cand_meta$candidate_id)))

con <- load_db()
DBI::dbExecute(con, "DELETE FROM candidate_microbe WHERE method LIKE 'meta_%'")
```

</details>

    [1] 668

<details class="code-fold">
<summary>Code</summary>

``` r
DBI::dbWriteTable(con, "candidate_microbe", df_cand_meta,
                  append = TRUE, row.names = FALSE)
df_meta_summary <- DBI::dbGetQuery(con,
  "SELECT method, COUNT(*) n, MIN(n_cohorts) min_c, MAX(n_cohorts) max_c
     FROM candidate_microbe GROUP BY method ORDER BY method")
close_db(con)

knitr::kable(df_meta_summary,
             caption = "candidate_microbe by method after meta merge")
```

</details>

| method        |   n | min_c | max_c |
|:--------------|----:|------:|------:|
| ancombc       |  14 |     1 |     1 |
| consensus     |  12 |     1 |     1 |
| deseq2        |  27 |     1 |     1 |
| maaslin2      |  13 |     1 |     1 |
| meta_random   | 334 |     2 |     6 |
| meta_stouffer | 334 |     2 |     6 |
| wilcoxon      |  30 |     1 |     1 |

candidate_microbe by method after meta merge

<details class="code-fold">
<summary>Code</summary>

``` r
df_meta_stouffer |>
  dplyr::arrange(q_value) |>
  head(10) |>
  dplyr::select(taxon_id, taxon_name, n_cohorts, direction,
                z_pool, p_value, q_value) |>
  knitr::kable(digits = 4,
               caption = "Top 10 species by Stouffer meta q_value")
```

</details>

| taxon_id | taxon_name | n_cohorts | direction | z_pool | p_value | q_value |
|:---|:---|---:|:---|---:|---:|---:|
| STR_Filifactor_alocis | Filifactor_alocis | 6 | disease_enriched | 7.9548 | 0 | 0 |
| STR_Porphyromonas_gingivalis | Porphyromonas_gingivalis | 6 | disease_enriched | 7.5043 | 0 | 0 |
| STR_Eubacterium_nodatum | Eubacterium_nodatum | 6 | disease_enriched | 7.4103 | 0 | 0 |
| STR_Treponema_denticola | Treponema_denticola | 6 | disease_enriched | 7.4217 | 0 | 0 |
| STR_Prevotella_intermedia | Prevotella_intermedia | 6 | disease_enriched | 7.3400 | 0 | 0 |
| STR_GGB4333_SGB5935 | GGB4333_SGB5935 | 6 | disease_enriched | 7.2086 | 0 | 0 |
| STR_Tannerella_forsythia | Tannerella_forsythia | 6 | disease_enriched | 7.1770 | 0 | 0 |
| STR_GGB10852_SGB17523 | GGB10852_SGB17523 | 6 | disease_enriched | 7.0037 | 0 | 0 |
| STR_Eubacterium_saphenum | Eubacterium_saphenum | 5 | disease_enriched | 6.7626 | 0 | 0 |
| STR_Porphyromonas_endodontalis | Porphyromonas_endodontalis | 6 | disease_enriched | 6.5941 | 0 | 0 |

Top 10 species by Stouffer meta q_value

<details class="code-fold">
<summary>Code</summary>

``` r
df_meta_random |>
  dplyr::arrange(q_value) |>
  head(10) |>
  dplyr::select(taxon_id, taxon_name, n_cohorts, direction,
                effect_size, se, p_value, q_value) |>
  knitr::kable(digits = 4,
               caption = "Top 10 species by random-effects meta q_value")
```

</details>

| taxon_id | taxon_name | n_cohorts | direction | effect_size | se | p_value | q_value |
|:---|:---|---:|:---|---:|---:|---:|---:|
| STR_Eubacterium_nodatum | Eubacterium_nodatum | 6 | disease_enriched | 0.5526 | 0.0470 | 0 | 0 |
| STR_Filifactor_alocis | Filifactor_alocis | 6 | disease_enriched | 0.6232 | 0.0548 | 0 | 0 |
| STR_Treponema_denticola | Treponema_denticola | 6 | disease_enriched | 0.5434 | 0.0477 | 0 | 0 |
| STR_Porphyromonas_gingivalis | Porphyromonas_gingivalis | 6 | disease_enriched | 0.5562 | 0.0540 | 0 | 0 |
| STR_Porphyromonas_endodontalis | Porphyromonas_endodontalis | 6 | disease_enriched | 0.5623 | 0.0581 | 0 | 0 |
| STR_Prevotella_koreensis | Prevotella_koreensis | 4 | disease_enriched | 0.5638 | 0.0625 | 0 | 0 |
| STR_Fretibacterium_fastidiosum | Fretibacterium_fastidiosum | 6 | disease_enriched | 0.4618 | 0.0531 | 0 | 0 |
| STR_Treponema_SGB69443 | Treponema_SGB69443 | 5 | disease_enriched | 0.4576 | 0.0560 | 0 | 0 |
| STR_GGB4333_SGB5935 | GGB4333_SGB5935 | 6 | disease_enriched | 0.5268 | 0.0681 | 0 | 0 |
| STR_Treponema_SGB3604 | Treponema_SGB3604 | 5 | disease_enriched | 0.4498 | 0.0598 | 0 | 0 |

Top 10 species by random-effects meta q_value

<details class="code-fold">
<summary>Code</summary>

``` r
write_csv(df_per_cohort,    path_target("meta_per_cohort_wilcoxon.csv"))
write_csv(df_meta_stouffer, path_target("meta_stouffer.csv"))
write_csv(df_meta_random,   path_target("meta_random.csv"))
```

</details>

## Files written

<details class="code-fold">
<summary>Code</summary>

``` r
knitr::kable(qproj::proj_dir_info(path_target(), tz = "CET"))
```

</details>

| path                      | type      |   size | modification_time   |
|:--------------------------|:----------|-------:|:--------------------|
| ancombc_full_results.csv  | file      | 16.86K | 2026-05-05 21:13:30 |
| candidate_microbe.csv     | file      | 16.05K | 2026-05-05 21:13:30 |
| consensus_species.csv     | file      |  1.74K | 2026-05-05 21:13:30 |
| deseq2_full_results.csv   | file      | 16.84K | 2026-05-05 21:13:30 |
| maaslin2_full_results.csv | file      | 16.87K | 2026-05-05 21:13:30 |
| maaslin2_out              | directory |    224 | 2026-05-05 21:13:25 |
| method_agreement.csv      | file      |    298 | 2026-05-05 21:13:30 |
| wilcoxon_full_results.csv | file      | 22.37K | 2026-05-05 21:13:30 |
