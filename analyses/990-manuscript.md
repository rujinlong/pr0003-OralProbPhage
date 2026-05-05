# Oral probiotic phage knowledge layer — MVP data report
Your Name
2026-05-05

- [<span class="toc-section-number">1</span> Data layer](#data-layer)
- [<span class="toc-section-number">2</span> Sources used in the
  MVP](#sources-used-in-the-mvp)
- [<span class="toc-section-number">3</span> MVP figures](#mvp-figures)
- [<span class="toc-section-number">4</span> Methods
  (skeleton)](#methods-skeleton)
- [<span class="toc-section-number">5</span> What’s deferred to ripple
  (manuscript-blocking
  items)](#whats-deferred-to-ripple-manuscript-blocking-items)
- [References](#references)
  - [<span class="toc-section-number">5.1</span> Files
    written](#files-written)

> [!NOTE]
>
> ### MVP scope
>
> This is the manuscript **stub**. It references the schema, summarises
> the contract-table populations, and links to the figures from 110. No
> introduction / discussion is written yet — that is `990-r1`. When
> ripples land enough cohorts (`020-r1`, `040-r1`) and methods
> (`050-r1`, `080-r1`), the prose will fold around this stub.

# Data layer

The contract layer defines five core objects plus a data-source registry
(\[Sample, TaxonProfile, CandidateMicrobe, PhageHostLink,
EvidencePacket\]). Schema version: 0.1; canonical YAML lives at
`analyses/data/010-schema-design/schema-v0.1.yml`.

| table             | pk           | n_columns |
|:------------------|:-------------|----------:|
| data_source       | source_id    |         7 |
| sample            | sample_id    |         7 |
| taxon_profile     | profile_id   |         7 |
| candidate_microbe | candidate_id |        11 |
| phage_host_link   | link_id      |         7 |
| evidence_packet   | packet_id    |         6 |

Schema v0.1 — six contract tables

# Sources used in the MVP

| source_id | citation_key | doi | source_kind | n_samples |
|:---|:---|:---|:---|---:|
| geng2024_PRJDB11203 | 2024-Universal_Geng | 10.1002/imt2.212 | bacteriome | 42 |
| geng2024_PRJNA230363 | 2024-Universal_Geng | 10.1002/imt2.212 | bacteriome | 28 |
| geng2024_PRJNA396840 | 2024-Universal_Geng | 10.1002/imt2.212 | bacteriome | 20 |
| geng2024_PRJNA678453 | 2024-Universal_Geng | 10.1002/imt2.212 | bacteriome | 59 |
| geng2024_PRJNA717815 | 2024-Universal_Geng | 10.1002/imt2.212 | bacteriome | 26 |
| geng2024_PRJNA932553 | 2024-Universal_Geng | 10.1002/imt2.212 | bacteriome | 48 |
| yahara2021 | 2021-Longread_Yahara | 10.1038/s41467-020-20199-9 | virome | 4 |

Registered MVP data sources

| table             |  rows |
|:------------------|------:|
| data_source       |     7 |
| sample            |   227 |
| taxon_profile     | 30063 |
| candidate_microbe |    96 |
| phage_host_link   |   770 |
| evidence_packet   |   344 |

Contract-layer populations after MVP

The bacteriome layer comes from the Geng et al. 2024 universal
periodontitis signature (**2024-Universal_Geng?**), specifically the
PRJDB11203 saliva sub-cohort. The virome layer comes from the Yahara et
al. 2021 long-read PromethION oral phageome (Yahara et al. 2021).

# MVP figures

- `analyses/data/090-m2a-mapping/fig_genus_mapping.png` — phages per
  M1a-implicated host genus.
- `analyses/data/110-evaluation/fig1_samples.png` — sample × source
  breakdown.
- `analyses/data/110-evaluation/fig2_volcano.png` — M1a Wilcoxon
  volcano.
- `analyses/data/110-evaluation/fig3_phage_coverage.png` —
  phage-coverage of M1a candidates.

# Methods (skeleton)

- M1a differential: Wilcoxon rank-sum on TSS relative abundances, BH
  FDR; filter species with prevalence ≥ 10%. See
  `analyses/050-m1a-differential-analysis.qmd`.
- M2a host prediction: CAT taxonomy of VIRSorter contigs from the
  published long-read assembly, restricted to genus-or-finer ranks. See
  `analyses/080-m2a-host-prediction.qmd`.
- Mapping: candidate species ↔ phage host via genus-prefix join. See
  `analyses/090-m2a-mapping.qmd`.

# What’s deferred to ripple (manuscript-blocking items)

- multi-cohort meta-analysis (`050-r2`) → upgrade `n_cohorts ≥ 1` MVP
  candidates to robust signature.
- Real `compatibility_score` (`100-r1`) → currently `NA`.
- Cross-cohort virome (`020-r2` + `070-r1` on SLURM) → currently a
  single long-read paper.

# References

<div id="refs" class="references csl-bib-body hanging-indent">

<div id="ref-2021-Longread_Yahara" class="csl-entry">

Yahara, Koji, Masato Suzuki, Aki Hirabayashi, et al. 2021. “Long-Read
Metagenomics Using PromethION Uncovers Oral Bacteriophages and Their
Interaction with Host Bacteria.” *Nature Communications* 12 (1): 27.
<https://doi.org/10.1038/s41467-020-20199-9>.

</div>

</div>

## Files written

<details class="code-fold">
<summary>Code</summary>

``` r
knitr::kable(qproj::proj_dir_info(path_target(), tz = "CET"))
```

</details>

| path | type | size | modification_time |
|:-----|:-----|-----:|:------------------|
