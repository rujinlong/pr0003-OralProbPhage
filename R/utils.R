#' Load project SQLite database
#'
#' Assumes the caller has run the qproj qmd boilerplate (`setwd(here::here("analyses"))` +
#' `here::i_am(<qmd-name>)`), which pins the here root at `analyses/`. The DB
#' therefore lives at `here::here("pr0003.sqlite")`.
#' @return A DBI connection to the project SQLite database
load_db <- function() {
  con <- DBI::dbConnect(RSQLite::SQLite(), here::here("pr0003.sqlite"))
  DBI::dbExecute(con, "PRAGMA foreign_keys = ON;")
  if (requireNamespace("connections", quietly = TRUE)) {
    connections::connection_view(con)
  }
  con
}

#' Close database connection
#' @param con A DBI connection
close_db <- function(con) {
  if (!missing(con) && !is.null(con) && DBI::dbIsValid(con)) {
    DBI::dbDisconnect(con)
  }
}

#' Write table to project database
#' @param con A DBI connection
#' @param name Table name (must already exist in the contract schema)
#' @param df Data frame to write
#' @param append If TRUE, append rows; if FALSE (default), truncate-then-insert
#' @return The input data frame (invisible)
write_table_db <- function(con, name, df, append = FALSE) {
  if (append) {
    DBI::dbWriteTable(con, name, df, append = TRUE, row.names = FALSE)
  } else {
    if (DBI::dbExistsTable(con, name)) {
      DBI::dbExecute(con, sprintf('DELETE FROM "%s"', name))
    }
    DBI::dbWriteTable(con, name, df, append = TRUE, row.names = FALSE)
  }
  invisible(df)
}

#' Read a contract table as tibble
#' @param con A DBI connection
#' @param name Table name
#' @return A tibble
read_table_db <- function(con, name) {
  tibble::as_tibble(DBI::dbReadTable(con, name))
}

#' Path to the canonical project schema YAML
#'
#' MVP version 0.2 lives under the 010 module's data target dir. Other helpers
#' default to this path when no explicit schema is given.
#' @param version Schema version string (default "0.2")
schema_path <- function(version = "0.2") {
  here::here("data", "010-schema-design",
             paste0("schema-v", version, ".yml"))
}

#' Initialize the contract schema from a YAML definition
#'
#' Creates the project tables in pr0003.sqlite via CREATE TABLE IF NOT EXISTS.
#' This is idempotent — running on an already-initialized DB is a no-op.
#' @param con A DBI connection (creates one via load_db() if NULL — caller closes)
#' @param schema_yaml Path to schema YAML
#' @return The DBI connection (invisible)
init_schema <- function(con = NULL, schema_yaml = schema_path()) {
  if (is.null(con)) con <- load_db()
  schema <- yaml::read_yaml(schema_yaml)
  for (tname in names(schema$tables)) {
    tdef <- schema$tables[[tname]]
    pk <- tdef$pk %||% ""
    col_defs <- vapply(names(tdef$columns), function(cname) {
      ctype <- tdef$columns[[cname]]
      suffix <- if (cname == pk) " PRIMARY KEY" else ""
      sprintf('"%s" %s%s', cname, ctype, suffix)
    }, character(1))
    sql <- sprintf('CREATE TABLE IF NOT EXISTS "%s" (%s)',
                   tname, paste(col_defs, collapse = ", "))
    DBI::dbExecute(con, sql)
  }
  invisible(con)
}

#' Assert that the live DB matches the YAML schema
#'
#' Checks: every YAML table exists in the DB, and every YAML column exists in
#' the live table. Extra columns in the DB (e.g. a future migration in
#' progress) do not fail. Throws on the first violation.
#' @param con A DBI connection
#' @param schema_yaml Path to schema YAML
#' @return TRUE invisibly on success
assert_schema <- function(con, schema_yaml = schema_path()) {
  schema <- yaml::read_yaml(schema_yaml)
  db_tables <- DBI::dbListTables(con)
  for (tname in names(schema$tables)) {
    if (!tname %in% db_tables) {
      stop(sprintf("Missing table in DB: %s (declared in %s)",
                   tname, schema_yaml), call. = FALSE)
    }
    expected_cols <- names(schema$tables[[tname]]$columns)
    actual_cols <- DBI::dbListFields(con, tname)
    missing <- setdiff(expected_cols, actual_cols)
    if (length(missing) > 0) {
      stop(sprintf("Table '%s' missing columns: %s",
                   tname, paste(missing, collapse = ", ")), call. = FALSE)
    }
  }
  invisible(TRUE)
}

#' Project schema version (read from canonical YAML)
schema_version <- function(schema_yaml = schema_path()) {
  yaml::read_yaml(schema_yaml)$schema_version
}

`%||%` <- function(a, b) if (is.null(a)) b else a
