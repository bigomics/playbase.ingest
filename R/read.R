#' Read all CSV file in folder or current directory
#'
#' @export
read_files <- function(dir = ".", pattern = NULL) {
  ff <- dir(dir, pattern = pattern)
  file1 <- head(grep("count|expression|abundance|concentration|intensity",
    ff,
    value = TRUE
  ), 1)
  file2 <- head(grep("sample", ff, value = TRUE), 1)
  file3 <- head(grep("contrast|comparison", ff, value = TRUE), 1)
  counts <- samples <- contrasts <- NULL
  if (length(file1)) counts <- read_counts(file.path(dir, file1))
  if (length(file2)) samples <- read_samples(file.path(dir, file2))
  if (length(file3)) contrasts <- read_contrasts(file.path(dir, file3))

  list(
    counts = counts,
    samples = samples,
    contrasts = contrasts
  )
}

#' Read counts data from file
#'
#' @param file string. path to file
#' @param drop_na_rows boolean. drop rows without rownames
#' @param first boolean. drop multiple feature names (separated by ;)
#' @param unique boolean. make duplicated rows unique by pasting a number
#'
#' @details This function reads a count matrix using \code{read.as_matrix()},
#' validates it with \code{validate_counts()}, and optionally converts row names from IDs to
#' gene symbols.
#'
#' It removes rows with NA, blank or invalid symbols, and collapses any duplicate symbols by
#' summing counts across rows.
#'
#' @examples
#' \dontrun{
#' counts <- read_counts(playbase.ingest::example_file("counts.csv"))
#' }
#' @export
read_counts <- function(file, first = FALSE, unique = FALSE, paste_char = "_") {
  if (is.character(file)) {
    df <- read.as_matrix(file, as.char = FALSE)
  } else if (is.matrix(file) || is.data.frame(file)) {
    df <- file
  } else {
    message("[read_counts] ERROR: Input must be filename or matrix")
    return(NULL)
  }

  ## determine column types (NEED RETHINK!)
  df1 <- type.convert(data.frame(head(df, 20), check.names = FALSE), as.is = TRUE)
  col.type <- sapply(df1, class)
  xannot.names <- "gene|^id$$|metabolite|compound|position|phospo.*site" ## possible numeric annotations
  is.xannot <- grepl(xannot.names, tolower(colnames(df)))
  char.cols <- which(col.type == "character" | is.xannot)
  last.charcol <- tail(char.cols, 1)
  last.charcol

  ## NEED RETHINK. if the rownames are not unique, and some more
  ## character columns exists, then search for best column and paste
  ## after rownames.
  if (length(char.cols) > 0 && sum(duplicated(rownames(df)))) {
    ndup <- sapply(char.cols, function(k) {
      sum(duplicated(paste0(rownames(df), "_", df[, k])))
    })
    ndup
    sel <- names(which.min(ndup))
    rownames(df) <- paste0(rownames(df), "_", df[, sel])
  }

  ## As expression values we take all columns after the last character
  ## column (if any).
  if (length(last.charcol)) {
    message("[read_counts] extra annotation columns = ", paste(1:last.charcol, collapse = " "))
    df <- df[, (last.charcol + 1):ncol(df)]
  }

  ## convert to numeric if needed (probably yes...)
  is.numeric.matrix <- all(apply(head(df, 20), 2, is.numeric))
  if (!is.numeric.matrix) {
    message("[read_counts] force to numeric values")
    rn <- rownames(df)
    suppressWarnings(df <- apply(df, 2, as.numeric))
    rownames(df) <- rn
  }

  ## when multiple feature names, should we take only first feature?
  if (first) rownames(df) <- first_feature(rownames(df))

  ## if rownames are duplicated, we append a number behing
  if (unique) rownames(df) <- make_unique(rownames(df))
  return(df)
}

#' Read samples data from file
#'
#' @param file Path to input sample data file. Should be a matrix with samples as rows and metadata as columns.
#'
#' @return dataframe the file with the data
#' @details This function reads the sample matrix with the meta-data information
#' of the counts and converts it to a dataframe.
#'
#' @examples
#' samples <- read_samples(playbase.ingest::example_file("samples.csv"))
#' @export
read_samples <- function(file) {
  df <- read.as_matrix(file)
  is_valid <- validate_samples(df)
  if (!is_valid) message("[read_samples] WARNING: Samples file has errors")
  return(as.data.frame(df))
}

#' Read contrasts data from file
#'
#' @param file string. path to file
#'
#' @return matrix. the file with the data
#'
#' @examples
#' contrasts <- read_contrasts(playbase.ingest::example_file("contrasts.csv"))
#' @export
read_contrasts <- function(file) {
  df <- read.as_matrix(file)
  is_valid <- validate_contrasts(df)
  if (!is_valid) message("[read_contrasts] WARNING: Contrasts file has errors")
  return(df)
}

#' Read gene/probe annotation file
#' @export
read_annot <- function(file, unique = TRUE) {
  if (is.character(file)) {
    ## we read without rownames because we want to retain full header
    df <- read.as_matrix(file, row.names = NULL)
    rownames(df) <- df[, 1]
  } else if (is.matrix(file) || is.data.frame(file)) {
    df <- file
  } else {
    message("[read_annot] ERROR: input not valid.")
    return(NULL)
  }

  ## add column title
  if (colnames(df)[1] == "") colnames(df)[1] <- "row.names"

  ## determine last character column
  df1 <- type.convert(data.frame(head(df, 20), check.names = FALSE), as.is = TRUE)
  col.type <- sapply(df1, class)
  xannot.names <- "gene|symbol|protein|compound|title|description|name|position"
  is.xannot <- grepl(xannot.names, tolower(colnames(df)))
  char.cols <- which(col.type == "character" | is.xannot)
  char.cols
  last.charcol <- tail(char.cols, 1)
  last.charcol

  ## if the rownames are not unique, and some more character columns
  ## exists, then search for best column and paste after rownames
  ## (first column).
  if (sum(duplicated(df[, 1])) && length(char.cols) >= 2) {
    ndup <- sapply(char.cols[-1], function(k) {
      sum(duplicated(paste0(df[, 1], "_", df[, k])))
    })
    ndup
    sel <- names(which.min(ndup))[1]
    rownames(df) <- paste0(rownames(df), "_", df[, sel])
  }

  ## drop numerical columns (these can be intensities). We check if we
  ## have equal or more than two columns. First column are rownames.
  if (length(char.cols) && last.charcol >= 2) {
    df <- df[, 1:last.charcol, drop = FALSE]
  } else {
    df <- NULL
  }
  return(df)
}

#' @export
first_feature <- function(x) {
  unname(sapply(strsplit(x, split = "[;,\\|]"), "[[", 1))
}
