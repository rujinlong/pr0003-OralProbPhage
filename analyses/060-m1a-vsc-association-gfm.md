Your Name

- [<span class="toc-section-number">1</span>
  060-m1a-vsc-association](#060-m1a-vsc-association)
  - [<span class="toc-section-number">1.1</span> Status](#status)
  - [<span class="toc-section-number">1.2</span> Files
    written](#files-written)

# 060-m1a-vsc-association

Your Name 2026-05-05

- [<span class="toc-section-number">1</span> Status](#status)
- [<span class="toc-section-number">2</span> Files
  written](#files-written)

**Updated: 2026-05-05 15:17:55 CET.**

> \[!NOTE\]
>
> ### Deferred to ripple
>
> VSC (volatile sulfur compound) and clinical-marker association is
> **not in the MVP scope**. The PRJDB11203 sub-cohort used in 040/050 is
> a periodontitis case-control without VSC measurements, so there is
> nothing to correlate at this point in the pipeline.
>
> **Ripple plan** (`060-r1`):
>
> - Identify cohorts with paired VSC / pocket depth / inflammation
>   measurements (the halitosis–VSC reviews referenced in
>   `analyses/background.md` give a starting list).
> - Register them in `data_source` (`020-r1`) and ingest in `040-r1`.
> - Run quantitative association: Spearman rank-correlation against VSC
>   concentration; partial correlations adjusting for age/sex; mediation
>   analysis (microbe → VSC → halitosis severity).
> - Write candidates with `phenotype = "halitosis"` into
>   `candidate_microbe`.

## Status

<details class="code-fold">

<summary>

Code
</summary>

``` r
tibble::tibble(
  module          = "060-m1a-vsc-association",
  status          = "deferred-to-ripple",
  ripple_targets  = "060-r1 (Spearman + partial correlation against VSC / pocket depth)",
  blocker_reason  = "MVP cohort has no VSC measurements",
  next_action     = "Register a halitosis cohort with VSC in 020-r1, then run this module"
) |> knitr::kable(caption = "Module status")
```

</details>

| module | status | ripple_targets | blocker_reason | next_action |
|:---|:---|:---|:---|:---|
| 060-m1a-vsc-association | deferred-to-ripple | 060-r1 (Spearman + partial correlation against VSC / pocket depth) | MVP cohort has no VSC measurements | Register a halitosis cohort with VSC in 020-r1, then run this module |

Module status

## Files written

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
