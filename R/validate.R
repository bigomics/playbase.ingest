getError <- function(e, what = "Description") {
  ERROR_MSG <- playbase.ingest::PGX_CHECKS
  if (!e %in% ERROR_MSG$error) {
    return(paste("unknown error", e))
  }
  ERROR_MSG[match(e, ERROR_MSG$error), what]
}

#' Validate counts data
#'
#' Counts data is valid if:
#'  - no duplicate rows
#'  - no empty rows
#'  - no duplicate cols
#'
#' @param df Matrix or data frame, or path to a file read with \code{read.as_matrix()}.
#'
#' @return boolean. true if data is valid
#' @export
validate_counts <- function(df) {
  if (is.character(df) && is.null(dim(df))) {
    df <- read.as_matrix(df)
  }
  chk <- pgx.checkINPUT(df, "COUNTS")
  ERROR_MSG <- playbase.ingest::PGX_CHECKS
  err <- names(chk$checks)
  if (length(err)) {
    msg <- lapply(err, function(e) getError(e))
    msg <- paste(msg, collapse = "; ")
    message("WARNING: ", msg)
  }
  chk$PASS
}

#' Validate samples data
#'
#' Samples data is valid if:
#'  - no duplicate rows
#'  - no empty rows
#'  - no duplicate cols
#'  - contains less than the max samples allowed
#'
#' @param df Matrix or data frame, or path to a file read with \code{read.as_matrix()}.
#'
#' @return boolean. true if data is valid
#' @export
validate_samples <- function(df) {
  if (is.character(df) && is.null(dim(df))) {
    df <- read.as_matrix(df)
  }
  chk <- pgx.checkINPUT(df, "SAMPLES")
  ERROR_MSG <- playbase.ingest::PGX_CHECKS
  err <- names(chk$checks)
  if (length(err)) {
    msg <- lapply(err, function(e) getError(e))
    msg <- paste(msg, collapse = "; ")
    message("WARNING: ", msg)
  }
  chk$PASS
}

#' Validate contrasts data
#'
#' Contrasts data is valid if:
#'  - no duplicate rows
#'  - no empty rows
#'  - no duplicate cols
#'  - contains only cols with "_vs_" in names
#'
#' @param df Matrix or data frame, or path to a file read with \code{read.as_matrix()}.
#'
#' @return boolean. true if data is valid
#' @export
validate_contrasts <- function(df) {
  if (is.character(df) && is.null(dim(df))) {
    df <- read.as_matrix(df)
  }
  chk <- pgx.checkINPUT(df, "CONTRASTS")
  ERROR_MSG <- playbase.ingest::PGX_CHECKS
  err <- names(chk$checks)
  if (length(err)) {
    msg <- lapply(err, function(e) getError(e))
    msg <- paste(msg, collapse = "; ")
    message("WARNING: ", msg)
  }
  chk$PASS
}
