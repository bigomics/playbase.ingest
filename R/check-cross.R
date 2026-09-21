#' Cross check input files for pgx.computePGX
#'
#' @param SAMPLE data.frame. The data frame corresponding to the input file as in playbase::SAMPLES
#' @param COUNTS data.frame. The data frame corresponding to the input file as in playbase::COUNTS
#' @param CONTRASTS data.frame. The data frame corresponding to the input file as in playbase::CONTRASTS
#'
#' @return a list with FIVE elements: SAMPLES, COUNTS and CONTRASTS that are the cleaned version of the
#'  input data frames, `checks` which contains the status of the checks, and
#'  `PASS` which contains the overall status of the check.
#' @export
pgx.crosscheckINPUT <- function(
  SAMPLES = NULL,
  COUNTS = NULL,
  CONTRASTS = NULL,
  PASS = TRUE
) {
  samples <- SAMPLES
  counts <- COUNTS
  contrasts <- CONTRASTS
  PASS <- PASS
  check_return <- list()

  if (!is.null(samples) && !is.null(counts)) {
    # Check that rownames(samples) match colnames(counts)
    COUNTS_NAMES_NOT_MATCHING_SAMPLES <- colnames(counts)[!colnames(counts) %in% rownames(samples)]

    # if there are not matches between samples and counts, return error e25
    if (length(COUNTS_NAMES_NOT_MATCHING_SAMPLES) == length(colnames(counts)) && PASS) {
      check_return$e25 <- COUNTS_NAMES_NOT_MATCHING_SAMPLES
      PASS <- FALSE
    }

    if (length(COUNTS_NAMES_NOT_MATCHING_SAMPLES) > 0 && PASS) {
      check_return$e21 <- COUNTS_NAMES_NOT_MATCHING_SAMPLES
      PASS <- TRUE
      # align counts columns with samples rownames
      samples_in_counts <- rownames(samples)[rownames(samples) %in% colnames(counts)]
      counts <- counts[, samples_in_counts, drop = FALSE]
    }
    SAMPLE_NAMES_NOT_MATCHING_COUNTS <- rownames(samples)[!rownames(samples) %in% colnames(counts)]

    if (length(SAMPLE_NAMES_NOT_MATCHING_COUNTS) > 0 && PASS) {
      check_return$e16 <- SAMPLE_NAMES_NOT_MATCHING_COUNTS
      PASS <- TRUE
      # align samples rows with counts colnames
      counts_in_samples <- colnames(counts)[colnames(counts) %in% rownames(samples)]
      samples <- samples[counts_in_samples, , drop = FALSE]
    }

    # Check that rownames(samples) match colnames(counts)
    SAMPLE_NAMES_PARTIAL_MATCHING_COUNTS <- intersect(
      rownames(samples),
      colnames(counts)
    )

    nsamples <- max(ncol(counts), nrow(samples))

    if (
      length(SAMPLE_NAMES_PARTIAL_MATCHING_COUNTS) > 0 &&
        length(SAMPLE_NAMES_PARTIAL_MATCHING_COUNTS) < nsamples &&
        PASS
    ) {
      TOTAL_SAMPLES_NAMES <- unique(c(rownames(samples), colnames(counts)))
      check_return$e19 <- TOTAL_SAMPLES_NAMES[!TOTAL_SAMPLES_NAMES %in% SAMPLE_NAMES_PARTIAL_MATCHING_COUNTS]
      samples <- samples[SAMPLE_NAMES_PARTIAL_MATCHING_COUNTS, , drop = FALSE]
      counts <- counts[, SAMPLE_NAMES_PARTIAL_MATCHING_COUNTS, drop = FALSE]
    }

    # Check that counts have the same order as samples.

    MATCH_SAMPLES_COUNTS_ORDER <- all(diff(match(rownames(samples), colnames(counts))) > 0)

    # in case no matches are found, we get an NA, which should be converted to FALSE
    if (is.na(MATCH_SAMPLES_COUNTS_ORDER)) {
      MATCH_SAMPLES_COUNTS_ORDER <- FALSE
    }

    if (!MATCH_SAMPLES_COUNTS_ORDER && PASS) {
      check_return$e18 <- "We will reorder your samples and counts."
      counts <- counts[, match(rownames(samples), colnames(counts))]
    }
  }

  if (!is.null(samples) && !is.null(contrasts)) {
    # Conver contrasts and Check that rows names of contrasts match rownames of samples.
    contrasts_check_results <- contrasts_conversion_check(samples, contrasts, PASS)

    if (contrasts_check_results$PASS == FALSE && PASS) {
      PASS <- FALSE
      check_return$e20 <- rownames(contrasts)[!rownames(contrasts) %in% rownames(samples)]
    }

    contrasts <- contrasts_check_results$CONTRASTS

    # Check that rownames(samples) match long contrast rownames.

    SAMPLE_NAMES_NOT_MATCHING_CONTRASTS <- NULL

    if (dim(contrasts)[1] > dim(samples)[1] && PASS) { # check that contrasts are in long format
      SAMPLE_NAMES_NOT_MATCHING_CONTRASTS <- c(
        setdiff(rownames(samples), rownames(contrasts)),
        setdiff(rownames(contrasts), rownames(samples))
      )
    }
    if (length(SAMPLE_NAMES_NOT_MATCHING_CONTRASTS) > 0 && PASS) {
      check_return$e17 <- SAMPLE_NAMES_NOT_MATCHING_CONTRASTS
    }
  }

  return(
    list(
      SAMPLES = samples,
      COUNTS = counts,
      CONTRASTS = contrasts,
      checks = check_return,
      PASS = PASS
    )
  )
}

#' Convert contrasts for OPG
#'
#' @param SAMPLE data.frame. The data frame corresponding to the input file as in playbase::SAMPLES
#' @param CONTRASTS data.frame. The data frame corresponding to the input file as in playbase::CONTRASTS
#' @param PASS boolean. The status of the checks.
#' @return converted contrast df
#' @export
contrasts_conversion_check <- function(SAMPLES, CONTRASTS, PASS) {
  samples1 <- SAMPLES
  contrasts1 <- contrasts.convertToLabelMatrix(CONTRASTS, SAMPLES)

  if (is.null(contrasts1)) {
    message("[contrasts_conversion_check] WARNING: could not convert contrasts!")
    return(list(CONTRASTS = CONTRASTS, PASS = FALSE))
  }

  ok.contrast <- length(intersect(rownames(samples1), rownames(contrasts1))) > 0
  if (ok.contrast && NCOL(contrasts1) > 0 && PASS) {
    # check that dimentions of contrasts match samples
    if (dim(contrasts1)[1] != dim(samples1)[1] && PASS) {
      message("[contrasts_conversion_check] WARNING: numrows of contrast1 do not match samples!")
      return(list(CONTRASTS = contrasts1, PASS = FALSE))
    }
    rownames(contrasts1) <- rownames(samples1)
    for (i in 1:ncol(contrasts1)) {
      isz <- (contrasts1[, i] %in% c(NA, "NA", "NA ", "", " ", "  ", "   ", " NA"))
      if (length(isz)) contrasts1[isz, i] <- NA
    }
  }
  return(list(CONTRASTS = contrasts1, PASS = PASS))
}
