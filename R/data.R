#' Get path to example dataset(s)
#'
#' `playbase.ingest` comes bundled with a number of sample files in its
#' `inst/extdata` directory. This function makes them easy to access. This
#' function was taken from tidyverse/readr.
#'
#' @param file string. Name of file. If `NULL`, the example files will
#'   be listed.
#' @return File path, or the file names when `file` is NULL.
#' @examples
#' example_file()
#' example_file("counts.csv")
#' @export
example_file <- function(file = NULL) {
  if (is.null(file)) {
    dir(system.file("extdata", package = "playbase.ingest"))
  } else {
    system.file("extdata", file, package = "playbase.ingest", mustWork = TRUE)
  }
}

#' Checks performed by pgx.checkINPUT
#'
#' @format ## `data.frame`
#' rows are checks, columns are description of the check performed.
"PGX_CHECKS"
