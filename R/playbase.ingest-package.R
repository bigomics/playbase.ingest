#' playbase.ingest: data ingestion for OmicsPlayground
#'
#' Readers and input validation split out of playbase.
#'
#' @import data.table
#' @importFrom utils head tail type.convert read.csv untar unzip
#' @importFrom stats as.formula sd
#' @importFrom methods is
#' @keywords internal
"_PACKAGE"

#' Checks performed by pgx.checkINPUT
#'
#' @format ## `data.frame`
#' rows are checks, columns are description of the check performed.
"PGX_CHECKS"
