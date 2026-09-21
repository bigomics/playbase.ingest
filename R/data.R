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
