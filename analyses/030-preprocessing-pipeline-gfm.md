Your Name

- [<span class="toc-section-number">1</span>
  030-preprocessing-pipeline](#030-preprocessing-pipeline)
  - [<span class="toc-section-number">1.1</span> Status](#status)
  - [<span class="toc-section-number">1.2</span> Files
    written](#files-written)

# 030-preprocessing-pipeline

Your Name 2026-05-05

- [<span class="toc-section-number">1</span> Status](#status)
- [<span class="toc-section-number">2</span> Files
  written](#files-written)

**Updated: 2026-05-05 15:17:49 CET.**

> \[!NOTE\]
>
> ### Deferred to ripple
>
> This module is intentionally a stub in the MVP. The MVP loads
> **published abundance and vOTU tables** directly (modules 040 / 070),
> so no fastq → species pipeline is required to validate the contract
> layer.
>
> **Ripple plan** (`030-r1` / `030-r2`):
>
> - `030-r1`: Nextflow pipeline on SLURM that re-processes one published
>   cohort (`geng2024_PRJDB11203` SRA fastq) with vpipe / Kraken2 /
>   MetaPhlAn4 → species table. Validate by comparing species relative
>   abundances against the published `Figure 01.Rdata` matrix; flag any
>   species with \> 25% disagreement.
> - `030-r2`: extend to virome on SLURM (ViWrap / ViOTUcluster / CheckV)
>   and reproduce the Yahara 2021 vOTU catalog from raw long reads.
>
> See `IMPLEMENTATION_PLAN.md` and
> `/Users/cmbjx/.claude/plans/mvp-linear-ripple.md` §Ripple roadmap for
> the full ripple ladder.

## Status

<details class="code-fold">

<summary>

Code
</summary>

``` r
tibble::tibble(
  module = "030-preprocessing-pipeline",
  status = "deferred-to-ripple",
  ripple_targets = "030-r1 (bacteriome fastq), 030-r2 (virome fastq)",
  next_action = "After MVP DoD passes, scaffold a Nextflow pipeline under nextflow/"
) |> knitr::kable(caption = "Module status")
```

</details>

| module | status | ripple_targets | next_action |
|:---|:---|:---|:---|
| 030-preprocessing-pipeline | deferred-to-ripple | 030-r1 (bacteriome fastq), 030-r2 (virome fastq) | After MVP DoD passes, scaffold a Nextflow pipeline under nextflow/ |

Module status

## Files written

These files have been written to the target directory,
data/030-preprocessing-pipeline:

<details class="code-fold">

<summary>

Code
</summary>

``` r
knitr::kable(qproj::proj_dir_info(path_target(), tz = "CET"))
```

</details>

| path | type | size | modification_time |
|:-----|:-----|-----:|:------------------|
