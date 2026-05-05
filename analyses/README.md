# analyses

This directory contains the analysis workflow.

## Structure

- Each analysis file writes only to its own subdirectory of `data/`.
- Files are named with a numeric prefix to indicate execution order,
  e.g. `01-import.qmd`, `02-clean.qmd`, `99-summary.qmd`.
  The `00-` prefix is reserved for the framework's `data/00-raw/`
  input region -- start your own steps at `01-` or higher.

## Usage

Create a new analysis file:

```r
qproj::use_qmd("01-import")
```
