# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Project type

This is a **qproj-based R/Quarto research pipeline** that doubles as an R package skeleton (`DESCRIPTION` + `NAMESPACE` + `R/`). It is not a standalone library — the package skeleton exists so `devtools::load_all()` can expose the helpers in `R/utils.R` to every `analyses/*.qmd` step. Treat the codebase as a numbered, sandboxed pipeline first and an R package second.

## Pipeline architecture

The pipeline lives entirely in `analyses/` and is organized as numbered `.qmd` modules. The number encodes both execution order and pipeline phase:

| Range | Phase | Purpose |
|-------|-------|---------|
| `010`–`030` | Infrastructure | Schema design (M0a's 5 core objects), data source registry, unified preprocessing |
| `040`–`060` | M1a — bacteriome | Data ingest, cross-cohort differential analysis, VSC/clinical association |
| `070`–`090` | M2a — virome | vOTU catalog, phage–host prediction, mapping back to M1a hits |
| `100`–`110` | Integration | Knowledge integration into `EvidencePacket`s, internal evaluation |
| `990` | Manuscript | Data paper draft (`gfm` output) |

The five core objects produced by the pipeline are: **Sample**, **TaxonProfile**, **CandidateMicrobe**, **PhageHostLink**, **EvidencePacket**. Anything new should be modeled in terms of these — they form the "data API" of the project.

See `analyses/background.md` for the scientific rationale and `IMPLEMENTATION_PLAN.md` (note: the filename has a leading space) for the phased roadmap.

## Sandbox / data lineage rules

These constraints come from qproj and are non-obvious:

- **One writer per directory.** Each step writes only to its own `analyses/data/<NNN-step-name>/` subdirectory. Reading from another step's data dir is fine; writing into it is a lineage violation.
- **`analyses/data/00-raw/` is read-only input.** Steps must not deposit derived files there. User-numbered steps start at `010-` (the `00-` prefix is reserved for the framework).
- **Persistent state lives in `analyses/pr0003.sqlite`.** This SQLite DB is the canonical store for derived tables that need to flow between steps. It is gitignored. Use the helpers in `R/utils.R`:
  - `load_db()` — open a connection (also registers it for the RStudio Connections pane)
  - `close_db(con)` — disconnect safely
  - `write_table_db(con, "stepNN_snake_case", df)` — write with the `stepNN_*` naming convention
- All `analyses/data/*` and `*.sqlite` are gitignored. Rendered `*.md` outputs from Quarto are also ignored (only `*.qmd` and `README.md` are tracked).

## qmd module conventions

Every `analyses/*.qmd` follows the same three-block header (`params` / `setup` / `packages`). When creating a new step or editing an existing one:

- Keep the `here::i_am(..., uuid = ...)` line — the UUID is per-file and pins the project root.
- Keep the `qproj::proj_create_dir_target(params$name, clean = FALSE)` call — this materializes the step's output dir.
- `path_target()` / `path_source()` / `path_raw` / `path_resource` / `path_data` are the conventional path helpers. Prefer them over hardcoded paths.
- The setup block calls `devtools::load_all()`, which loads `R/utils.R`. Add new shared helpers there rather than redefining inside qmds.
- `conflicted::conflicts_prefer(...)` resolves common dplyr/data.table/stats clashes — extend the existing call instead of suppressing warnings ad-hoc.

## Common commands

```r
# Set up the workflow (one-time, already done)
qproj::proj_use_workflow("analyses")

# Create a new analysis module (numbering encodes phase — see table above)
qproj::use_qmd("NNN-name", "analyses")

# Render a single module (run from project root)
quarto::quarto_render("analyses/050-m1a-differential-analysis.qmd")

# Render the entire workflow
quarto::quarto_render("analyses")

# Open the project DB for ad-hoc inspection
con <- load_db(); DBI::dbListTables(con); close_db(con)
```

Output format is `gfm` (per `analyses/_quarto.yml`); rendered `.md` files are gitignored.

## Bibliography

`analyses/references.bib` is a **symlink** to `/Users/cmbjx/storage/references.bib` (a shared, machine-local Zotero export). Don't replace it with a real file or commit a copy — manuscripts on other machines will resolve the symlink to the same shared library. If the symlink target is missing, ask before creating a placeholder.

## Status

The repo is at "infrastructure ready, no analyses run yet": all 12 qmd modules are scaffolded but empty (just the boilerplate header with a placeholder `data-import` chunk). The next concrete work is filling in `010-schema-design.qmd`, `020-data-source-registry.qmd`, then `030-preprocessing-pipeline.qmd`. Don't assume any step has produced data yet — check the SQLite DB and `analyses/data/` before referencing prior outputs.
