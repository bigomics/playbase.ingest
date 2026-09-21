#' Check an uploaded input table before a pgx object is built
#'
#' @param df data.frame. The data frame corresponding to the input file as
#'   described in Omics Playground documentation
#' @param type character. One of "SAMPLES", "COUNTS", "EXPRESSION", "CONTRASTS" - The type
#'  of data file to check.
#'
#' @return a list with two elements: `checks` which contains the status of the checks, and
#'  `PASS` which contains the overall status of the check.
#' @export
pgx.checkINPUT <- function(
  df,
  type = c("SAMPLES", "COUNTS", "EXPRESSION", "CONTRASTS")
) {
  datatype <- match.arg(type)
  df_clean <- df
  PASS <- TRUE
  check_return <- list()

  if (datatype == "COUNTS" || datatype == "EXPRESSION") {
    ## Sparse matrices (e.g. dgCMatrix from h5 scRNA-seq reads) cannot have character
    ## columns by definition, so skip the character-class check entirely. This avoids
    ## both the O(n_cols) vapply overhead and R's internal as.matrix() coercion that
    ## apply() triggers on any S4 object.
    if (!inherits(df_clean, "sparseMatrix")) {
      df.class <- apply(df_clean, 2, class)
      if (any(df.class == "character")) {
        # remove special characters before converting to numeric (keep
        # commas, dots) as they define the decimal separator
        jj <- which(df.class == "character")
        df_clean[, jj] <- apply(df_clean[, jj], 2, function(x) gsub("[^0-9,.]", "", x))

        # convert matrix from character to numeric, sometimes we receive
        # character matrix from read.csv function
        df_clean[, jj] <- apply(df_clean[, jj], 2, as.numeric, simplify = TRUE)
      }
    }

    rownames(df_clean) <- rownames(df)
    sample_names <- colnames(df_clean)

    # check if any value in df had a non-NA value that was converted to NA in df_clean.
    # Skip for sparse matrices: the character-stripping conversion that introduces NAs
    # is skipped for sparse input, so the result is always empty and !is.na() on a
    # sparse matrix densifies (structural FALSE -> TRUE covers the full matrix).
    if (!inherits(df_clean, "sparseMatrix")) {
      ANY_NON_NUMERIC <- which(!is.na(df) & is.na(df_clean), arr.ind = TRUE)
      if (length(ANY_NON_NUMERIC) > 0 && PASS) {
        check_return$e27 <- paste("gene:", rownames(ANY_NON_NUMERIC), " and ", "sample:", colnames(df_clean)[ANY_NON_NUMERIC[, 2]])
      }
    }

    # check if there are any infinite values. replace infinite values
    # in counts by NA.
    # For sparse matrices scan only the stored non-zero values (@x slot) to avoid
    # a full sparse->dense comparison that creates a dense logical of the same size.
    if (inherits(df_clean, "sparseMatrix")) {
      inf_in_x <- which(is.infinite(df_clean@x))
      ANY_INFINITE <- if (length(inf_in_x) > 0) {
        Matrix::which(is.infinite(df_clean), arr.ind = TRUE)
      } else {
        integer(0)
      }
    } else {
      ANY_INFINITE <- which(df_clean == Inf | df_clean == -Inf, arr.ind = TRUE)
    }

    if (length(ANY_INFINITE) > 0 && PASS) {
      ## replace Inf with NA
      df_clean[ANY_INFINITE] <- NA
      check_return$e28 <- paste(
        "gene:", rownames(ANY_INFINITE), " and ",
        "sample:", colnames(df_clean)[ANY_INFINITE[, 2]]
      )
    }

    # check for duplicated colnanes (gives error)
    ANY_DUPLICATED <- unique(sample_names[which(duplicated(sample_names))])

    if (length(ANY_DUPLICATED) > 0 && PASS) {
      PASS <- FALSE
      check_return$e6 <- ANY_DUPLICATED
    }

    feature_names <- rownames(df_clean)

    # check for duplicated rownames (but pass)
    ANY_DUPLICATED <- unique(feature_names[which(duplicated(feature_names))])

    if (length(ANY_DUPLICATED) > 0 && PASS) {
      check_return$e7 <- ANY_DUPLICATED
    }

    # check for zero count rows, remove them
    ANY_ROW_ZERO <- which(Matrix::rowSums(df_clean, na.rm = TRUE) == 0, )

    if (length(ANY_ROW_ZERO) > 0 && PASS) {
      # get the row names with all zeros
      zero.rows <- names(ANY_ROW_ZERO)

      # remove the rownames with all zeros
      # df_clean <- df_clean[!(rownames(df_clean) %in% zero.rows), , drop = FALSE]

      nzerorows <- length(ANY_ROW_ZERO)
      err.mesg <- zero.rows
      check_return$e9 <- err.mesg
    }

    # check for zero count columns, remove them
    ANY_COLUMN_ZERO <- which(Matrix::colSums(df_clean) == 0)

    if (length(ANY_COLUMN_ZERO) > 0 && PASS) {
      check_return$e10 <- names(ANY_COLUMN_ZERO)
      # remove the column names with all zeros by using check_return$e9
      df_clean <- df_clean[, !(colnames(df_clean) %in% check_return$e10), drop = FALSE]
    }

    # check if counts are log transformed
    check.log <- is_logged(df_clean)
    if (check.log) {
      check_return$e29 <- "Possible log transformed counts detected."
    }

    # check min amount of features
    n.features <- nrow(df_clean)
    if (n.features < 3) {
      check_return$e31 <- "Too few features detected. Minimum is 3."
      PASS <- FALSE
    }
  }

  if (datatype == "SAMPLES") {
    feature_names <- rownames(df_clean)

    # check for duplicated rownames
    ANY_DUPLICATED <- unique(feature_names[which(duplicated(feature_names))])

    if (length(x = ANY_DUPLICATED) > 0 && PASS) {
      check_return$e1 <- ANY_DUPLICATED
      PASS <- FALSE
    }
  }

  if (datatype == "CONTRASTS") {
    feature_names <- rownames(df_clean)

    # check that contrasts has at least one column

    COMPARISONS_WITHOUT_COLUMNS <- dim(df_clean)[2] == 0

    if (COMPARISONS_WITHOUT_COLUMNS && PASS) {
      check_return$e26 <- "No columns provided in comparisons."
      PASS <- FALSE
    }

    # check for duplicated rownames (but pass)
    ANY_DUPLICATED <- unique(feature_names[which(duplicated(feature_names))])

    if (length(x = ANY_DUPLICATED) > 0 && PASS) {
      check_return$e11 <- ANY_DUPLICATED
      PASS <- FALSE
    }

    ## check that numerator_vs_denominator is in the contrasts
    if (all(grepl(" vs ", colnames(df_clean)))) {
      colnames(df_clean) <- gsub(" vs ", "_vs_", colnames(df_clean))
    }
    has_no_vs <- which(!grepl("_vs_", colnames(df_clean)))
    if (length(has_no_vs) > 0 && PASS) {
      check_return$e24 <- colnames(df_clean)[has_no_vs]
      PASS <- FALSE
    }

    if (PASS) {
      # Split the column names at "_vs_"
      split_names <- strsplit(colnames(df_clean), "_vs_")

      # Get the numerators and denominators
      numerators <- sapply(split_names, "[[", 1)
      split_numerators <- strsplit(numerators, ":")

      # if colon is present in numerators, keep the elements after colon
      numerators <- sapply(split_numerators, function(x) {
        if (length(x) > 1) {
          x[2]
        } else {
          x[1]
        }
      })

      denominators <- sapply(split_names, "[[", 2)
      ## Check if all elements in the matrix are character
      all_numeric <- any(apply(df_clean, c(1, 2), is.numeric))

      if (!all_numeric && PASS) {
        ## only run if we have characters in matrix
        escape_special <- function(x) { # Special characters that grep considers operators
          special_chars <- c(".", "\\", "|", "(", ")", "[", "]", "{", "}", "^", "$", "*", "+", "?")
          gsub("([\\.^$*+?(){}|\\[\\]])", "\\\\\\1", x, perl = TRUE)
        }

        numerators <- sapply(numerators, escape_special)
        denominators <- sapply(denominators, escape_special)

        COLUMN_IN_GROUPS <- sapply(1:length(denominators), function(i) {
          vv <- setdiff(df_clean[, i], c(NA, "", " ", "NA"))
          all(grepl(paste0("^", numerators[i], "|^", denominators[i]), vv))
        })
        CONTRASTS_IN_GROUPS <- COLUMN_IN_GROUPS

        if (all(!CONTRASTS_IN_GROUPS) && PASS) {
          check_return$e23 <- "All comparisons were invalid."
          PASS <- FALSE
        }

        if (any(!CONTRASTS_IN_GROUPS) && PASS) {
          check_return$e22 <- colnames(df_clean)[!CONTRASTS_IN_GROUPS]
          df_clean <- df_clean[, CONTRASTS_IN_GROUPS, drop = FALSE]
        }
      }
    } ## if PASS
  }

  # general checks for all data datatypes

  # check for empty df
  IS_DF_EMPTY <- any(dim(df_clean) == 0)

  if (!IS_DF_EMPTY && !inherits(df_clean, "sparseMatrix")) {
    df_clean <- as.matrix(df_clean)
  }

  if (IS_DF_EMPTY && PASS) {
    check_return$e15 <- "empty dataframe"
    PASS <- FALSE
  }

  return(
    list(
      df = df_clean,
      checks = check_return,
      PASS = PASS
    )
  )
}
