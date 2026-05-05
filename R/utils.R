#' Load project SQLite database
#' @return A DBI connection to the project SQLite database
load_db <- function() {
  con <- DBI::dbConnect(RSQLite::SQLite(), here::here("analyses", "pr0003.sqlite"))
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
#' @param name Table name (format: stepNN_snake_case)
#' @param df Data frame to write
#' @param overwrite Whether to overwrite existing table
#' @return The input data frame (invisible)
write_table_db <- function(con, name, df, overwrite = TRUE) {
  DBI::dbWriteTable(con, name, df, overwrite = overwrite)
  invisible(df)
}
