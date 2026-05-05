# 010-schema-design
Your Name
2026-05-05

- [<span class="toc-section-number">1</span> Schema
  definition](#schema-definition)
- [<span class="toc-section-number">2</span> Controlled
  vocabularies](#controlled-vocabularies)
- [<span class="toc-section-number">3</span> ID
  conventions](#id-conventions)
- [<span class="toc-section-number">4</span> Initialize the
  database](#initialize-the-database)
- [<span class="toc-section-number">5</span> Files
  written](#files-written)

**Updated: 2026-05-05 21:13:20 CET.**

Defines the contract layer for the project — six SQLite tables in
`pr0003.sqlite` that materialize the five core objects (Sample,
TaxonProfile, CandidateMicrobe, PhageHostLink, EvidencePacket) plus a
`data_source` registry. Schema version v0.1 (MVP). Downstream modules
read/write these tables only.

<details class="code-fold">
<summary>Code</summary>

``` r
suppressPackageStartupMessages({
  library(here)
  library(conflicted)
  library(tidyverse)
  library(data.table)
  library(yaml)
  library(DBI)
  library(RSQLite)
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

## Schema definition

The canonical schema is built here as an R list and serialized to YAML.
The YAML file under `data/010-schema-design/schema-v0.1.yml` is the
human-readable copy that `R/utils.R::init_schema()` and
`assert_schema()` read.

<details class="code-fold">
<summary>Code</summary>

``` r
schema <- list(
  schema_version = "0.1",
  tables = list(
    data_source = list(
      pk = "source_id",
      columns = list(
        source_id     = "TEXT",
        citation_key  = "TEXT",
        doi           = "TEXT",
        source_kind   = "TEXT",
        n_samples     = "INTEGER",
        table_url     = "TEXT",
        notes         = "TEXT"
      )
    ),
    sample = list(
      pk = "sample_id",
      columns = list(
        sample_id      = "TEXT",
        source_id      = "TEXT",
        subject_id     = "TEXT",
        disease_status = "TEXT",
        body_site      = "TEXT",
        seq_type       = "TEXT",
        meta_json      = "TEXT"
      )
    ),
    taxon_profile = list(
      pk = "profile_id",
      columns = list(
        profile_id     = "TEXT",
        sample_id      = "TEXT",
        taxon_id       = "TEXT",
        taxon_kind     = "TEXT",
        taxon_name     = "TEXT",
        abundance      = "REAL",
        abundance_kind = "TEXT"
      )
    ),
    candidate_microbe = list(
      pk = "candidate_id",
      columns = list(
        candidate_id   = "TEXT",
        taxon_id       = "TEXT",
        taxon_name     = "TEXT",
        phenotype      = "TEXT",
        direction      = "TEXT",
        effect_size    = "REAL",
        p_value        = "REAL",
        q_value        = "REAL",
        n_cohorts      = "INTEGER",
        method         = "TEXT",
        schema_version = "TEXT"
      )
    ),
    phage_host_link = list(
      pk = "link_id",
      columns = list(
        link_id        = "TEXT",
        vOTU_id        = "TEXT",
        host_taxon_id  = "TEXT",
        evidence_kind  = "TEXT",
        confidence     = "REAL",
        support_json   = "TEXT",
        schema_version = "TEXT"
      )
    ),
    evidence_packet = list(
      pk = "packet_id",
      columns = list(
        packet_id            = "TEXT",
        candidate_id         = "TEXT",
        link_id              = "TEXT",
        compatibility_score  = "REAL",
        summary_text         = "TEXT",
        generated_at         = "TEXT"
      )
    )
  )
)

path_schema <- path_target(paste0("schema-v", schema$schema_version, ".yml"))
yaml::write_yaml(schema, path_schema)
```

</details>

## Controlled vocabularies

Columns whose values are restricted to a small set. SQLite does not
enforce ENUMs natively; downstream modules should validate before
insert.

<details class="code-fold">
<summary>Code</summary>

``` r
df_vocab <- tibble::tribble(
  ~table,              ~column,          ~allowed_values,
  "data_source",       "source_kind",    "bacteriome | virome",
  "sample",            "disease_status", "healthy | periodontitis | gingivitis | caries | halitosis | other",
  "sample",            "body_site",      "saliva | supragingival_plaque | subgingival_plaque | mucosa | tongue | other",
  "sample",            "seq_type",       "16S | shotgun | virome",
  "taxon_profile",     "taxon_kind",     "bacterium | virus",
  "taxon_profile",     "abundance_kind", "rel_abund | count | rpkm | tpm | presence",
  "candidate_microbe", "phenotype",      "periodontitis | gingivitis | caries | halitosis",
  "candidate_microbe", "direction",      "health_enriched | disease_enriched",
  "candidate_microbe", "method",         "wilcoxon | deseq2 | maaslin2 | ancombc | meta_re",
  "phage_host_link",   "evidence_kind",  "crispr | similarity | tetra | ml | curated | published_prediction"
)
knitr::kable(df_vocab, caption = "Controlled vocabularies (MVP)")
```

</details>

| table | column | allowed_values |
|:---|:---|:---|
| data_source | source_kind | bacteriome \| virome |
| sample | disease_status | healthy \| periodontitis \| gingivitis \| caries \| halitosis \| other |
| sample | body_site | saliva \| supragingival_plaque \| subgingival_plaque \| mucosa \| tongue \| other |
| sample | seq_type | 16S \| shotgun \| virome |
| taxon_profile | taxon_kind | bacterium \| virus |
| taxon_profile | abundance_kind | rel_abund \| count \| rpkm \| tpm \| presence |
| candidate_microbe | phenotype | periodontitis \| gingivitis \| caries \| halitosis |
| candidate_microbe | direction | health_enriched \| disease_enriched |
| candidate_microbe | method | wilcoxon \| deseq2 \| maaslin2 \| ancombc \| meta_re |
| phage_host_link | evidence_kind | crispr \| similarity \| tetra \| ml \| curated \| published_prediction |

Controlled vocabularies (MVP)

## ID conventions

- **Bacterial taxa**: `HMT_<id>` (eHOMD) preferred → `NCBI_<taxid>`
  fallback → `STR_<slug>` last resort.
- **Viral taxa (vOTU)**: `vOTU_<source>_<index>` — preserves source
  provenance.
- **Sample**: `SAMP_<source>_<sra_or_index>`.
- **Candidate / link / packet**: synthetic short hash from natural keys.

## Initialize the database

<details class="code-fold">
<summary>Code</summary>

``` r
con <- load_db()
init_schema(con, schema_yaml = path_schema)
assert_schema(con, schema_yaml = path_schema)
df_tables <- tibble::tibble(
  table = DBI::dbListTables(con)
) |>
  dplyr::mutate(
    n_columns = purrr::map_int(table, ~ length(DBI::dbListFields(con, .x))),
    n_rows    = purrr::map_int(table, ~ as.integer(
      DBI::dbGetQuery(con, sprintf('SELECT COUNT(*) AS n FROM "%s"', .x))$n
    ))
  )
close_db(con)
knitr::kable(df_tables, caption = "Contract tables in pr0003.sqlite after init")
```

</details>

| table             | n_columns | n_rows |
|:------------------|----------:|-------:|
| candidate_microbe |        11 |     96 |
| data_source       |         7 |      7 |
| evidence_packet   |         6 |    344 |
| phage_host_link   |         7 |    770 |
| sample            |         7 |    227 |
| taxon_profile     |         7 |  30063 |

Contract tables in pr0003.sqlite after init

All six contract tables exist and are empty — downstream modules will
populate them. `assert_schema()` is the gate; subsequent qmds should
call it before writing.

## Files written

These files have been written to the target directory,
data/010-schema-design:

<details class="code-fold">
<summary>Code</summary>

``` r
knitr::kable(qproj::proj_dir_info(path_target(), tz = "CET"))
```

</details>

| path            | type |  size | modification_time   |
|:----------------|:-----|------:|:--------------------|
| schema-v0.1.yml | file | 1.32K | 2026-05-05 21:13:21 |
