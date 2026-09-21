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

