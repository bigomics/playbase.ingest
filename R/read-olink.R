#' Read Olink NPX data and create a counts matrix
#' @param NPX_data Path to Olink file. Must be standard format as per OlinkAnalyze R package.
#' @return NPX data matrix (features on rows; samples on columns)
#' @return Sample metadata matrix (samples on rows; metadata on columns)
#' @export
read_Olink_NPX <- function(NPX_data) {
  if (!requireNamespace("OlinkAnalyze", quietly = TRUE)) {
    stop("read_Olink_NPX() requires the 'OlinkAnalyze' package (Suggests); install it to read Olink NPX files.")
  }
  NPX <- try(OlinkAnalyze::read_NPX(NPX_data), silent = TRUE)
  if (inherits(NPX, "try-error")) {
    message("[read_Olink_NPX]: Uploaded file does not adhere with standard Olink format.")
    return(NULL)
  }

  NPX <- data.table::as.data.table(NPX)
  cols <- tolower(colnames(NPX))

  npx.id <- colnames(NPX)[grep("npx", cols)[1]]
  ss.id <- colnames(NPX)[grep("sampleid", cols)[1]]
  ff.id <- colnames(NPX)[c(grep("uniprot", cols), grep("assay$", cols))[1]]

  if (is.na(npx.id)) message("[read_Olink_NPX]: 'NPX' is missing.")
  if (is.na(ss.id)) message("[read_Olink_NPX]: 'SampleID' is missing.")
  if (is.na(ff.id)) message("[read_Olink_NPX]: 'Uniprot' or 'Assay' is missing.")
  if (is.na(npx.id) | is.na(ss.id) | is.na(ff.id)) {
    return(NULL)
  }

  ## Counts
  fm <- as.formula(paste0(ff.id, "~", ss.id))
  counts.df <- data.table::dcast(NPX, fm, value.var = npx.id, fun.aggregate = mean)
  counts <- as.matrix(counts.df[, -1, with = FALSE])
  rownames(counts) <- counts.df[[1]]
  counts <- counts[!is.na(rownames(counts)), , drop = FALSE]

  ## Metadata
  NPX <- as.data.frame(NPX)
  hh <- grepl("uniprot|olinkid|assay|npx|freq|lod", cols)
  meta_cols <- colnames(NPX)[!hh]
  samples <- NPX[!duplicated(NPX[[ss.id]]), meta_cols, drop = FALSE]
  rownames(samples) <- samples[[ss.id]]
  samples <- samples[, setdiff(colnames(samples), ss.id), drop = FALSE]

  cm <- intersect(colnames(counts), rownames(samples))
  counts <- counts[, cm, drop = FALSE]
  samples <- samples[cm, , drop = FALSE]

  return(list(counts = counts, samples = samples))
}
