## --------------------------------------------------------------------
## --------------------------------------------------------------------
## --------------------------------------------------------------------
#' Read data from a GMT file
#'
#' This function reads data from a GMT (Gene Matrix Transposed) file format.
#' The GMT format is commonly used to store gene sets or gene annotations.
#'
#' @param gmt.file The path to the GMT file.
#' @param dir (Optional) The directory where the GMT file is located.
#' @param add.source (Optional) Specifies whether to include the source information in the gene sets' names.
#' @param nrows (Optional) The number of rows to read from the GMT file.
#'
#' @export
#'
#' @return A list of gene sets, where each gene set is represented as a character vector of gene names.
#'
read.gmt <- function(gmt.file, dir = NULL, add.source = FALSE, nrows = -1) {
  f0 <- gmt.file
  if (strtrim(gmt.file, 1) == "/") dir <- NULL
  if (!is.null(dir)) f0 <- paste(sub("/$", "", dir), "/", gmt.file, sep = "")
  gmt <- utils::read.csv(f0, sep = "!", header = FALSE, comment.char = "#", nrows = nrows)[, 1]
  gmt <- as.character(gmt)
  gmt <- sapply(gmt, strsplit, split = "\t")
  names(gmt) <- NULL
  gmt.name <- sapply(gmt, "[", 1)
  gmt.source <- sapply(gmt, "[", 2)
  gmt.genes <- sapply(gmt, function(x) {
    if (length(x) < 3) {
      return("")
    }
    paste(x[3:length(x)], collapse = " ")
  })
  gset <- strsplit(gmt.genes, split = "[ \t]")
  gset <- lapply(gset, function(x) setdiff(x, c("", "NA", NA)))
  names(gset) <- gmt.name

  if (add.source) {
    names(gset) <- paste0(names(gset), " (", gmt.source, ")")
  }
  gset
}
