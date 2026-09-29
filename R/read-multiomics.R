#' Read and merge per-layer counts into one multi-omics matrix
#'
#' @param counts List of counts, one per omics layer. Each element is a file
#'   path or a matrix, as accepted by \code{read_counts()}.
#' @param datatypes Character vector, one data type per element of
#'   \code{counts}: "RNA-seq" (or microarray / transcriptomics), "proteomics"
#'   or "metabolomics".
#'
#' @details Each layer's feature names get the layer prefix (\code{gx:},
#' \code{px:}, \code{mx:}) that downstream preprocessing uses to infer layers.
#' Only samples present in every layer are kept.
#'
#' @return \code{list(counts = , dropped_samples = )}: the stacked numeric
#'   matrix (features x shared samples) and the sample names that were removed
#'   because they are missing from at least one layer.
#'
#' @examples
#' \dontrun{
#' mo <- read_multiomics_counts(
#'   list("rna.csv", "prot.csv"),
#'   c("RNA-seq", "proteomics")
#' )
#' }
#' @export
read_multiomics_counts <- function(counts, datatypes) {
  tag <- "[playbase.ingest::read_multiomics_counts]"
  if (length(counts) < 2) stop(tag, " at least two layers are required")
  if (length(counts) != length(datatypes)) {
    stop(tag, " counts and datatypes must have the same length")
  }
  prefix <- vapply(datatypes, multiomics_prefix, character(1), USE.NAMES = FALSE)
  if (anyDuplicated(prefix)) stop(tag, " each data type can be used only once")

  mats <- lapply(counts, read_counts)
  all_samples <- unique(unlist(lapply(mats, colnames)))
  shared <- Reduce(intersect, lapply(mats, colnames))
  if (length(shared) == 0) stop(tag, " the layers share no sample")
  dropped <- setdiff(all_samples, shared)
  if (length(dropped)) {
    message(tag, " dropping samples not present in every layer: ", paste(dropped, collapse = ", "))
  }

  mats <- Map(function(m, p) {
    m <- m[, shared, drop = FALSE]
    rownames(m) <- paste0(p, ":", rownames(m))
    m
  }, mats, prefix)
  list(counts = do.call(rbind, unname(mats)), dropped_samples = dropped)
}

## Layer prefix for a data type; same matching as the OmicsPlayground upload.
multiomics_prefix <- function(datatype) {
  dt <- tolower(datatype)
  if (grepl("proteomics", dt)) return("px")
  if (grepl("metabolomics", dt)) return("mx")
  if (grepl("rna|microarray|micro.array|transcriptomics", dt)) return("gx")
  stop("[playbase.ingest::read_multiomics_counts] unsupported data type: ", datatype)
}
