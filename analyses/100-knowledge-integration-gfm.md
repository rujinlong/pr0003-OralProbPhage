Your Name

- [<span class="toc-section-number">1</span>
  100-knowledge-integration](#100-knowledge-integration)
  - [<span class="toc-section-number">1.1</span> Pull the prefix mapping
    from 090 + full contract
    tables](#pull-the-prefix-mapping-from-090--full-contract-tables)
  - [<span class="toc-section-number">1.2</span> Parse `support_json` to
    get lifestyle (`integrase`) and
    `virsorter_cat`](#parse-support_json-to-get-lifestyle-integrase-and-virsorter_cat)
  - [<span class="toc-section-number">1.3</span> Compute compatibility
    score](#compute-compatibility-score)
  - [<span class="toc-section-number">1.4</span> Build
    packets](#build-packets)
  - [<span class="toc-section-number">1.5</span> Score distribution +
    top/bottom packets](#score-distribution--topbottom-packets)
  - [<span class="toc-section-number">1.6</span> Write](#write)
  - [<span class="toc-section-number">1.7</span> Persist](#persist)
  - [<span class="toc-section-number">1.8</span> Files
    written](#files-written)

# 100-knowledge-integration

Your Name 2026-05-05

- [<span class="toc-section-number">1</span> Pull the prefix mapping
  from 090 + full contract
  tables](#pull-the-prefix-mapping-from-090--full-contract-tables)
- [<span class="toc-section-number">2</span> Parse `support_json` to get
  lifestyle (`integrase`) and
  `virsorter_cat`](#parse-support_json-to-get-lifestyle-integrase-and-virsorter_cat)
- [<span class="toc-section-number">3</span> Compute compatibility
  score](#compute-compatibility-score)
- [<span class="toc-section-number">4</span> Build
  packets](#build-packets)
- [<span class="toc-section-number">5</span> Score distribution +
  top/bottom packets](#score-distribution--topbottom-packets)
- [<span class="toc-section-number">6</span> Write](#write)
- [<span class="toc-section-number">7</span> Persist](#persist)
- [<span class="toc-section-number">8</span> Files
  written](#files-written)

**Updated: 2026-05-05 21:13:37 CET.**

Builds `evidence_packet` rows by joining `candidate_microbe` ×
`phage_host_link` on the genus-prefix mapping from 090. Each row of the
join becomes one packet that bundles a candidate species and one of the
phages predicted to infect its genus.

**Ripple `100-r1` (this version)**: a real `compatibility_score` is now
computed per packet that captures therapeutic priority. Four orthogonal
terms are summed in raw space (each capped) and squashed through
`plogis()` to land in (0, 1):

- **Direction** (±1.0): `+1.0` if the host is `disease_enriched` (we
  want to lyse it); `-1.0` if `health_enriched` (lysing would harm the
  patient).
- **Lifestyle** (0 / +0.5): `+0.5` for lytic-leaning phages
  (`integrase = "no"` or NA); `0.0` for temperate/integrase-bearing
  phages to deprioritize HGT-prone candidates.
- **Confidence** (0 / +1.0): scaled directly from CAT host-prediction
  confidence (already in 0-1).
- **Effect size** (-0.5..+0.5):
  `sign(effect_size) * min(0.5, 0.1 * abs(effect_size))`.

`summary_text` now embeds the score and key drivers so downstream
readers (or the manuscript stub in 990) can read the packet without
joining tables.

**Ripple ladder**:

- `100-r1`: real compatibility scoring — combine direction, lifestyle,
  host-prediction confidence, and effect-size magnitude (this version).
- `100-r2`: therapeutic ranking — boost packets where the host is a
  known K12/M18-suppressed taxon or a halitosis-VSC producer; add CRISPR
  spacer and receptor evidence as additional score terms.

<details class="code-fold">

<summary>

Code
</summary>

``` r
suppressPackageStartupMessages({
  library(here)
  library(conflicted)
  library(tidyverse)
  library(digest)
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

## Pull the prefix mapping from 090 + full contract tables

The 090 mapping CSV gives us the candidate × link join keys; we re-read
`candidate_microbe` and `phage_host_link` from the contract DB so that
we have access to every column (in particular `support_json`, which
carries the lifestyle/integrase signal that the score depends on).

<details class="code-fold">

<summary>

Code
</summary>

``` r
path_in <- here::here("data", "090-m2a-mapping", "mapping_prefix.csv")
stopifnot(file.exists(path_in))
df_pairs <- readr::read_csv(path_in, show_col_types = FALSE)

con <- load_db()
df_candidate <- read_table_db(con, "candidate_microbe")
df_link      <- read_table_db(con, "phage_host_link")
close_db(con)

cat("Pairs to package:", nrow(df_pairs),
    "| candidate_microbe rows:", nrow(df_candidate),
    "| phage_host_link rows:", nrow(df_link), "\n")
```

</details>

    Pairs to package: 344 | candidate_microbe rows: 96 | phage_host_link rows: 770 

## Parse `support_json` to get lifestyle (`integrase`) and `virsorter_cat`

`support_json` is a per-link TEXT blob written by 080. We parse it once,
row-wise, and keep the two fields that feed scoring (`integrase`,
`virsorter_cat`). Missing or malformed JSON is tolerated — those rows
fall back to `NA`, which the scoring step treats as the neutral /
lytic-leaning default per the spec.

<details class="code-fold">

<summary>

Code
</summary>

``` r
parse_support <- function(s) {
  if (is.null(s) || is.na(s) || !nzchar(s)) {
    return(list(integrase = NA_character_, virsorter_cat = NA_integer_))
  }
  out <- tryCatch(jsonlite::fromJSON(s, simplifyVector = TRUE),
                  error = function(e) NULL)
  if (is.null(out)) {
    return(list(integrase = NA_character_, virsorter_cat = NA_integer_))
  }
  list(
    integrase     = if (!is.null(out$integrase))     as.character(out$integrase) else NA_character_,
    virsorter_cat = if (!is.null(out$virsorter_cat)) suppressWarnings(as.integer(out$virsorter_cat)) else NA_integer_
  )
}

df_support <- df_link |>
  dplyr::mutate(
    parsed = purrr::map(support_json, parse_support),
    integrase     = vapply(parsed, function(x) x$integrase,     character(1)),
    virsorter_cat = vapply(parsed, function(x) x$virsorter_cat, integer(1))
  ) |>
  dplyr::select(link_id, integrase, virsorter_cat)

cat("Support-json parsed for", nrow(df_support),
    "links | integrase distribution:\n")
```

</details>

    Support-json parsed for 770 links | integrase distribution:

<details class="code-fold">

<summary>

Code
</summary>

``` r
print(dplyr::count(df_support, integrase))
```

</details>

    # A tibble: 3 × 2
      integrase     n
      <chr>     <int>
    1 no           29
    2 yes          32
    3 <NA>        709

## Compute compatibility score

For each packet (candidate × link) we compute four independent terms in
raw space, sum them, and squash through `plogis()` so the score lands in
(0, 1). Higher = more therapeutically promising.

| Term | Range | Logic |
|----|----|----|
| direction | ±1.0 | +1.0 disease_enriched, -1.0 health_enriched |
| lifestyle | 0.0 / +0.5 | +0.5 if integrase no / NA (lytic-leaning), 0 if yes |
| confidence | 0.0 .. +1.0 | CAT confidence (already in 0-1) |
| effect-size | -0.5 .. +0.5 | sign(es) \* min(0.5, 0.1 \* |

If any input is missing for a packet, the missing term defaults to 0 and
a warning is emitted. The final `compatibility_score` is finite numeric
for every packet (no NA leaks).

<details class="code-fold">

<summary>

Code
</summary>

``` r
direction_term <- function(direction) {
  out <- dplyr::case_when(
    direction == "disease_enriched" ~  1.0,
    direction == "health_enriched"  ~ -1.0,
    TRUE                            ~  0.0
  )
  if (any(is.na(direction) | !direction %in% c("disease_enriched", "health_enriched"))) {
    n_bad <- sum(is.na(direction) | !direction %in% c("disease_enriched", "health_enriched"))
    warning(sprintf("direction_term: %d packet(s) had unrecognized direction; defaulted to 0", n_bad))
  }
  out
}

lifestyle_term <- function(integrase) {
  # +0.5 for lytic-leaning ("no" or NA). 0.0 for temperate ("yes").
  # Anything else falls to 0 with a warning (treated as temperate-equivalent
  # to be conservative — we only credit lifestyle when it's clearly lytic).
  out <- dplyr::case_when(
    is.na(integrase)        ~  0.5,
    integrase == "no"       ~  0.5,
    integrase == "yes"      ~  0.0,
    TRUE                    ~  0.0
  )
  bad <- !is.na(integrase) & !integrase %in% c("yes", "no")
  if (any(bad)) {
    warning(sprintf("lifestyle_term: %d packet(s) had unrecognized integrase value; defaulted to 0", sum(bad)))
  }
  out
}

confidence_term <- function(confidence) {
  out <- dplyr::if_else(is.finite(confidence), as.numeric(confidence), 0.0)
  if (any(!is.finite(confidence))) {
    warning(sprintf("confidence_term: %d packet(s) had non-finite confidence; defaulted to 0",
                    sum(!is.finite(confidence))))
  }
  # Clamp into [0, 1] just in case.
  pmin(pmax(out, 0.0), 1.0)
}

effect_size_term <- function(effect_size) {
  out <- ifelse(is.finite(effect_size),
                sign(effect_size) * pmin(0.5, 0.1 * abs(effect_size)),
                0.0)
  if (any(!is.finite(effect_size))) {
    warning(sprintf("effect_size_term: %d packet(s) had non-finite effect_size; defaulted to 0",
                    sum(!is.finite(effect_size))))
  }
  out
}

df_pairs_enriched <- df_pairs |>
  dplyr::left_join(df_support, by = "link_id") |>
  dplyr::mutate(
    term_direction  = direction_term(direction),
    term_lifestyle  = lifestyle_term(integrase),
    term_confidence = confidence_term(confidence),
    term_effect     = effect_size_term(effect_size),
    raw_score       = term_direction + term_lifestyle + term_confidence + term_effect,
    compatibility_score = plogis(raw_score)
  )

stopifnot(all(is.finite(df_pairs_enriched$compatibility_score)))
stopifnot(all(df_pairs_enriched$compatibility_score > 0 &
              df_pairs_enriched$compatibility_score < 1))

cat("Compatibility-score summary:\n")
```

</details>

    Compatibility-score summary:

<details class="code-fold">

<summary>

Code
</summary>

``` r
print(summary(df_pairs_enriched$compatibility_score))
```

</details>

       Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
     0.3056  0.4314  0.4598  0.5012  0.4991  0.9399 

## Build packets

`packet_id` is a deterministic hash of `candidate_id|link_id` (xxhash64,
12-char prefix), so re-renders are stable: the same packet always gets
the same ID across runs. The summary text now embeds the score and the
four key drivers (direction, host genus, CAT confidence, lifestyle).

<details class="code-fold">

<summary>

Code
</summary>

``` r
generated_at <- format(Sys.time(), "%Y-%m-%dT%H:%M:%S%z")

df_evidence_packet <- df_pairs_enriched |>
  dplyr::mutate(
    lifestyle_label = dplyr::case_when(
      is.na(integrase)   ~ "lytic",
      integrase == "no"  ~ "lytic",
      integrase == "yes" ~ "temperate",
      TRUE               ~ "unknown"
    ),
    summary_text = sprintf(
      "score=%.2f - %s (%s in %s, q=%.3g, log2FC=%.2f) <-> phage %s (CAT confidence %.2f, %s) on host %s",
      compatibility_score,
      taxon_name,
      direction,
      phenotype,
      q_value,
      effect_size,
      vOTU_id,
      confidence,
      lifestyle_label,
      sub("^STR_", "", host_taxon_id)
    ),
    generated_at = generated_at,
    packet_id = paste0(
      "PKT_",
      vapply(paste(candidate_id, link_id, sep = "|"),
             function(s) substr(digest::digest(s, algo = "xxhash64"), 1, 12),
             character(1))
    )
  ) |>
  dplyr::distinct(packet_id, .keep_all = TRUE) |>
  dplyr::select(packet_id, candidate_id, link_id, compatibility_score,
                summary_text, generated_at)

cat("evidence_packet rows:", nrow(df_evidence_packet), "\n")
```

</details>

    evidence_packet rows: 344 

<details class="code-fold">

<summary>

Code
</summary>

``` r
df_evidence_packet |> head(8) |>
  dplyr::transmute(packet_id, compatibility_score, summary_text) |>
  knitr::kable(caption = "First 8 evidence packets",
               digits = c(NA, 3, NA))
```

</details>

| packet_id | compatibility_score | summary_text |
|:---|---:|:---|
| PKT_ff2433f223b2 | 0.437 | score=0.44 - Corynebacterium_durum (health_enriched in periodontitis, q=0.0535, log2FC=-3.04) \<-\> phage vOTU_yahara2021_00394 (CAT confidence 0.55, lytic) on host Corynebacterium |
| PKT_5835f3a7459c | 0.531 | score=0.53 - Corynebacterium_durum (health_enriched in periodontitis, q=0.0535, log2FC=-3.04) \<-\> phage vOTU_yahara2021_01663 (CAT confidence 0.93, lytic) on host Corynebacterium |
| PKT_8980697bea1d | 0.911 | score=0.91 - Treponema_maltophilum (disease_enriched in periodontitis, q=0.0849, log2FC=3.91) \<-\> phage vOTU_yahara2021_00333 (CAT confidence 0.93, temperate) on host Treponema |
| PKT_627aa2c0500e | 0.917 | score=0.92 - Porphyromonas_endodontalis (disease_enriched in periodontitis, q=0.0861, log2FC=3.24) \<-\> phage vOTU_yahara2021_00146 (CAT confidence 0.58, lytic) on host Porphyromonas |
| PKT_2c92a357fa93 | 0.919 | score=0.92 - Porphyromonas_endodontalis (disease_enriched in periodontitis, q=0.0861, log2FC=3.24) \<-\> phage vOTU_yahara2021_01698 (CAT confidence 0.60, lytic) on host Porphyromonas |
| PKT_bbfb6415fea1 | 0.931 | score=0.93 - Tannerella_forsythia (disease_enriched in periodontitis, q=0.0973, log2FC=1.77) \<-\> phage vOTU_yahara2021_00295 (CAT confidence 0.92, lytic) on host Tannerella |
| PKT_3dc0b5a5ec8b | 0.928 | score=0.93 - Tannerella_forsythia (disease_enriched in periodontitis, q=0.0973, log2FC=1.77) \<-\> phage vOTU_yahara2021_00422 (CAT confidence 0.88, lytic) on host Tannerella |
| PKT_23759afda343 | 0.932 | score=0.93 - Tannerella_forsythia (disease_enriched in periodontitis, q=0.0973, log2FC=1.77) \<-\> phage vOTU_yahara2021_00443 (CAT confidence 0.94, lytic) on host Tannerella |

First 8 evidence packets

## Score distribution + top/bottom packets

<details class="code-fold">

<summary>

Code
</summary>

``` r
df_score_long <- df_pairs_enriched |>
  dplyr::transmute(
    packet_hash = paste(candidate_id, link_id, sep = "|"),
    direction,
    integrase,
    confidence,
    effect_size,
    term_direction, term_lifestyle, term_confidence, term_effect,
    raw_score,
    compatibility_score
  )

p_score <- ggplot(df_score_long,
                  aes(x = compatibility_score, fill = direction)) +
  geom_histogram(bins = 40, alpha = 0.85, position = "stack") +
  scale_fill_manual(values = c(disease_enriched = "#c84a4a",
                               health_enriched  = "#3f86b8")) +
  labs(x = "compatibility_score",
       y = "packets",
       title = "Distribution of compatibility scores",
       subtitle = sprintf("n = %d packets across %d candidates and %d phages",
                          nrow(df_score_long),
                          dplyr::n_distinct(df_pairs_enriched$candidate_id),
                          dplyr::n_distinct(df_pairs_enriched$vOTU_id))) +
  theme_minimal(base_size = 11)
ggsave(path_target("fig_score_distribution.png"),
       p_score, width = 6.5, height = 3.5, dpi = 150)
print(p_score)
```

</details>

![](100-knowledge-integration_files/figure-commonmark/score-dist-1.png)

<details class="code-fold">

<summary>

Code
</summary>

``` r
df_stats <- tibble::tibble(
  metric = c("min", "median", "mean", "max", "sd"),
  value  = c(min(df_pairs_enriched$compatibility_score),
             median(df_pairs_enriched$compatibility_score),
             mean(df_pairs_enriched$compatibility_score),
             max(df_pairs_enriched$compatibility_score),
             sd(df_pairs_enriched$compatibility_score))
)
knitr::kable(df_stats, digits = 4,
             caption = "compatibility_score summary statistics")
```

</details>

| metric |  value |
|:-------|-------:|
| min    | 0.3056 |
| median | 0.4598 |
| mean   | 0.5012 |
| max    | 0.9399 |
| sd     | 0.1470 |

compatibility_score summary statistics

<details class="code-fold">

<summary>

Code
</summary>

``` r
df_top10 <- df_pairs_enriched |>
  dplyr::arrange(dplyr::desc(compatibility_score)) |>
  dplyr::slice_head(n = 10) |>
  dplyr::transmute(taxon_name,
                   vOTU_id,
                   integrase = ifelse(is.na(integrase), "NA", integrase),
                   confidence,
                   effect_size,
                   direction,
                   compatibility_score)
knitr::kable(df_top10, digits = 3,
             caption = "Top 10 packets by compatibility_score")
```

</details>

| taxon_name | vOTU_id | integrase | confidence | effect_size | direction | compatibility_score |
|:---|:---|:---|---:|---:|:---|---:|
| Mogibacterium_timidum | vOTU_yahara2021_01545 | NA | 0.94 | 3.098 | disease_enriched | 0.940 |
| Mogibacterium_timidum | vOTU_yahara2021_01710 | NA | 0.94 | 3.098 | disease_enriched | 0.940 |
| Mogibacterium_timidum | vOTU_yahara2021_00189 | NA | 0.92 | 3.098 | disease_enriched | 0.939 |
| Mogibacterium_timidum | vOTU_yahara2021_00361 | NA | 0.92 | 3.098 | disease_enriched | 0.939 |
| Mogibacterium_timidum | vOTU_yahara2021_00550 | NA | 0.92 | 3.098 | disease_enriched | 0.939 |
| Mogibacterium_timidum | vOTU_yahara2021_00611 | NA | 0.92 | 3.098 | disease_enriched | 0.939 |
| Mogibacterium_timidum | vOTU_yahara2021_00635 | NA | 0.92 | 3.098 | disease_enriched | 0.939 |
| Mogibacterium_timidum | vOTU_yahara2021_00636 | NA | 0.92 | 3.098 | disease_enriched | 0.939 |
| Mogibacterium_timidum | vOTU_yahara2021_00644 | NA | 0.92 | 3.098 | disease_enriched | 0.939 |
| Mogibacterium_timidum | vOTU_yahara2021_00568 | NA | 0.89 | 3.098 | disease_enriched | 0.937 |

Top 10 packets by compatibility_score

<details class="code-fold">

<summary>

Code
</summary>

``` r
df_bot5 <- df_pairs_enriched |>
  dplyr::arrange(compatibility_score) |>
  dplyr::slice_head(n = 5) |>
  dplyr::transmute(taxon_name,
                   vOTU_id,
                   integrase = ifelse(is.na(integrase), "NA", integrase),
                   confidence,
                   effect_size,
                   direction,
                   compatibility_score)
knitr::kable(df_bot5, digits = 3,
             caption = "Bottom 5 packets by compatibility_score")
```

</details>

| taxon_name | vOTU_id | integrase | confidence | effect_size | direction | compatibility_score |
|:---|:---|:---|---:|---:|:---|---:|
| Actinomyces_massiliensis | vOTU_yahara2021_00915 | yes | 0.59 | -4.107 | health_enriched | 0.306 |
| Streptococcus_sanguinis | vOTU_yahara2021_00917 | yes | 0.60 | -4.060 | health_enriched | 0.309 |
| Streptococcus_sanguinis | vOTU_yahara2021_00926 | yes | 0.60 | -4.060 | health_enriched | 0.309 |
| Streptococcus_sanguinis | vOTU_yahara2021_01512 | yes | 0.62 | -4.060 | health_enriched | 0.313 |
| Streptococcus_sanguinis | vOTU_yahara2021_00938 | yes | 0.65 | -4.060 | health_enriched | 0.320 |

Bottom 5 packets by compatibility_score

<details class="code-fold">

<summary>

Code
</summary>

``` r
readr::write_csv(df_top10, path_target("score_top10.csv"))
readr::write_csv(df_bot5,  path_target("score_bottom5.csv"))
readr::write_csv(df_stats, path_target("score_stats.csv"))
```

</details>

## Write

Partition-based write: explicitly `DELETE FROM evidence_packet`, then
append the freshly built rows. This keeps the table contents in lockstep
with the current scoring run while preserving the table definition (no
DROP, no schema churn).

<details class="code-fold">

<summary>

Code
</summary>

``` r
con <- load_db()
assert_schema(con)
DBI::dbExecute(con, "DELETE FROM evidence_packet")
```

</details>

    [1] 344

<details class="code-fold">

<summary>

Code
</summary>

``` r
DBI::dbWriteTable(con, "evidence_packet", df_evidence_packet,
                  append = TRUE, row.names = FALSE)

# FK spot-check
n_orphan_cand <- DBI::dbGetQuery(con,
  "SELECT COUNT(*) AS n FROM evidence_packet ep
   LEFT JOIN candidate_microbe cm ON ep.candidate_id = cm.candidate_id
   WHERE cm.candidate_id IS NULL")$n
n_orphan_link <- DBI::dbGetQuery(con,
  "SELECT COUNT(*) AS n FROM evidence_packet ep
   LEFT JOIN phage_host_link phl ON ep.link_id = phl.link_id
   WHERE phl.link_id IS NULL")$n
stopifnot(n_orphan_cand == 0L, n_orphan_link == 0L)

# Score-shape sanity: every packet must have a finite score in (0, 1).
db_score_check <- DBI::dbGetQuery(con,
  "SELECT MIN(compatibility_score) AS min_s,
          AVG(compatibility_score) AS avg_s,
          MAX(compatibility_score) AS max_s,
          SUM(CASE WHEN compatibility_score IS NULL THEN 1 ELSE 0 END) AS n_null,
          COUNT(*) AS n_total
     FROM evidence_packet")
stopifnot(db_score_check$n_null == 0L,
          db_score_check$min_s > 0, db_score_check$max_s < 1)

n <- db_score_check$n_total
close_db(con)
cat("evidence_packet rows in DB:", n,
    "(0 orphan candidates, 0 orphan links, 0 NA scores)\n")
```

</details>

    evidence_packet rows in DB: 344 (0 orphan candidates, 0 orphan links, 0 NA scores)

<details class="code-fold">

<summary>

Code
</summary>

``` r
cat(sprintf("score range: [%.4f, %.4f] | mean: %.4f\n",
            db_score_check$min_s, db_score_check$max_s, db_score_check$avg_s))
```

</details>

    score range: [0.3056, 0.9399] | mean: 0.5012

## Persist

<details class="code-fold">

<summary>

Code
</summary>

``` r
write_csv(df_evidence_packet, path_target("evidence_packet.csv"))
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

| path                       | type |   size | modification_time   |
|:---------------------------|:-----|-------:|:--------------------|
| evidence_packet.csv        | file | 93.08K | 2026-05-05 21:13:15 |
| fig_score_distribution.png | file | 27.03K | 2026-05-05 21:13:15 |
| score_bottom5.csv          | file |    628 | 2026-05-05 21:13:15 |
| score_stats.csv            | file |    132 | 2026-05-05 21:13:15 |
| score_top10.csv            | file |  1.12K | 2026-05-05 21:13:15 |
