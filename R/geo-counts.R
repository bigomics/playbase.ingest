#' @describeIn pgx.getGEOcounts Download count data from GEO. First check
#' if the GEO ID is in archs5, then in recount. If not, try to get from GEO.
#' @return List of counts matrix and source.
#' @param id GEO accession ID.
#' @param archs.h5 Path to archs.h5 dataset.
#' @export
pgx.getGEOcounts <- function(accession, archs.h5) {
  id <- accession
  is.valid.id <- is.GEO.id.valid(id)
  if (!is.valid.id) stop("[pgx.getGEOcounts] FATAL: ID is invalid. Exiting.")
  id <- as.character(id)

  expr <- NULL
  meta <- NULL
  source <- ""

  if (!is.null(archs.h5) && is.null(expr)) {
    message("[pgx.getGEOcounts]: pgx.getGEOcounts.archs4...")
    expr <- pgx.getGEOcounts.archs4(id, archs.h5)
    if (!is.null(expr)) source <- "ARCHS4"
  }

  if (is.null(expr)) {
    message("[pgx.getGEOcounts]: pgx.getGEOcounts.GEOquery...")
    geo <- pgx.getGEOcounts.GEOquery(accession = id)
    if (!is.null(geo)) {
      expr <- geo[["expr"]]
      meta <- geo[["meta"]]
      source <- "GEO"
    }
  }

  if (is.null(expr)) {
    message("[pgx.getGEOcounts]: pgx.getGEOcounts.recount...")
    expr <- pgx.getGEOcounts.recount(accession = id)
    if (!is.null(expr)) source <- "recount"
  }

  if (is.null(expr)) {
    message("[pgx.getGEOcounts]: pgx.getGEOcounts.arrayexpress...")
    ae.data <- pgx.getArrayExpress.data(accession = id)
    if (!is.null(ae.data)) {
      expr <- ae.data[["expr"]]
      meta <- ae.data[["samples"]]
      source <- "ArrayExpress"
    }
  }

  if (is.null(expr)) {
    cat("WARNING:: Could not retrieve dataset. Please download manually.\n")
    return(NULL)
  }

  LL <- list(expr = expr, samples = meta, source = source)
  return(LL)
}

## -------------------------------------------------------------------------------------
## Query GEO expression
## -------------------------------------------------------------------------------------

#' @describeIn pgx.getGEOcounts.archs4 Downloads and extracts gene expression count
#' data from a GEO series stored in an HDF5 file (if available). It searches the
#' HDF5 file metadata to find samples matching the input GEO series ID, and returns
#' the count matrix for those samples. It detects log2-scale and convert to linear.
#' It also removes duplicated genes by summing in the linear scale
#' @export
pgx.getGEOcounts.archs4 <- function(id, h5.file) {
  is.valid.id <- is.GEO.id.valid(id)
  if (!is.valid.id) message("[pgx.getGEOcounts.archs4] Dataset ID is invalid. Please use a valid ID.")
  id <- as.character(id)

  if (is.null(h5.file) || h5.file == "") {
    stop("[pgx.getGEOcounts.archs4] FATAL: invalid path to h5.file ID. Exiting.\n")
  }

  sample.series <- rhdf5::h5read(h5.file, "meta/Sample_series_id")
  sample.series <- strsplit(as.character(sample.series), split = "Xx-xX")
  idx <- which(sapply(sample.series, function(s) id %in% s))
  if (!id %in% sample.series) {
    message("[pgx.getGEOcounts.archs4] WARNING: series ", id, " not in ARCHS4. Exiting.\n")
    return(NULL)
  }
  message("[pgx.getGEOcounts.archs4] Series ", id, " found in ARCHS4.")

  ## get matrix
  counts <- rhdf5::h5read(h5.file, "data/expression", index = list(NULL, idx))
  sample.acc <- rhdf5::h5read(h5.file, "meta/Sample_geo_accession")
  gene_name <- rhdf5::h5read(h5.file, "meta/genes")
  colnames(counts) <- sample.acc[idx]
  rownames(counts) <- gene_name

  ## ensure counts
  qx <- as.numeric(stats::quantile(counts, c(0., 0.25, 0.5, 0.75, 0.99, 1.0), na.rm = T))
  is.count <- (qx[5] > 100) || (qx[6] - qx[1] > 50 && qx[2] > 0) ||
    (qx[2] > 0 && qx[2] < 1 && qx[4] > 1 && qx[4] < 2) ## from GEO2R script
  if (!is.count) counts <- 2**counts

  ## rm missing genes and sum linear intensities
  jj <- which(!is.na(gene_name) & gene_name != "")
  counts <- counts[jj, ]
  gene_name <- gene_name[jj]
  counts <- tapply(1:nrow(counts), gene_name, function(ii) {
    Matrix::colSums(counts[ii, , drop = FALSE], na.rm = TRUE)
  })
  counts <- do.call(rbind, counts)
  message("[pgx.getGEOcounts.archs4] Success!")

  return(counts)
}

#' @describeIn pgx.getGEOcounts.recount Downloads and processes gene-level count data
#' for a GEO series from the recount database. It takes a GEO ID, searches recount,
#' downloads the RangedSummarizedExperiment object, and returns the count matrix.
#' It detects log2-scale and convert to linear. It also removes duplicated genes
#' by summing in the linear scale.
#' Vignette recount-quickstart.html
#' @export
pgx.getGEOcounts.recount <- function(accession) {
  id <- accession
  is.valid.id <- is.GEO.id.valid(id)
  if (!is.valid.id) stop("[pgx.getGEOcounts.recount] FATAL: ID is invalid. Exiting.")
  id <- as.character(id)

  project_info <- recount::abstract_search(id)
  pid <- project_info$project
  if (length(pid) == 0) {
    message("[pgx.getGEOcounts.recount] WARNING: series ", id, " not in recount. Exiting.\n")
    return(NULL)
  }
  message("[pgx.getGEOcounts.recount] Series ", id, " found in recount database.")

  ## Download the gene-level RangedSummarizedExperiment data
  outdir <- file.path(tempdir(), pid)
  cc <- try(recount::download_study(pid, outdir = outdir), silent = TRUE)
  if (inherits(cc, "try-error")) {
    message("[pgx.getGEOcounts.recount] Error: could not retrieve ", id, ". Exiting.\n")
    return(NULL)
  }

  ## Load the data
  load(file.path(outdir, "rse_gene.Rdata"))

  ## Scale counts by taking into account the total coverage per sample
  rse <- recount::scale_counts(rse_gene)
  counts <- MultiAssayExperiment::assay(rse)

  ## ensure counts
  qx <- as.numeric(stats::quantile(counts, c(0., 0.25, 0.5, 0.75, 0.99, 1.0), na.rm = T))
  is.count <- (qx[5] > 100) || (qx[6] - qx[1] > 50 && qx[2] > 0) ||
    (qx[2] > 0 && qx[2] < 1 && qx[4] > 1 && qx[4] < 2) ## from GEO2R script
  if (!is.count) counts <- 2**counts

  ## rm missing genes and sum linear intensities
  jj <- which(!is.na(rownames(counts)) & rownames(counts) != "")
  counts <- counts[jj, ]
  counts <- tapply(1:nrow(counts), rownames(counts), function(ii) {
    Matrix::colSums(counts[ii, , drop = FALSE], na.rm = TRUE)
  })
  counts <- do.call(rbind, counts)
  message("[pgx.getGEOcounts.recount] Success!")

  return(counts)
}

#' @describeIn pgx.getArrayExpress.data retrieves expression count data for a
#' accession ID using the arrayExpress R package. It downloads the counts and
#' metadata. It detects log2-scale and convert to linear.
#' \name {pgx.getArrayExpress.data}
#' \title {Get ArrayExpress Data}
#' @export
pgx.getArrayExpress.data <- function(accession) {
  id <- accession
  valid.ID <- is.GEO.id.valid(id)
  if (!valid.ID) {
    message("[pgx.getArrayExpress.data]: No valid ArrayExpress accession ID")
    return(NULL)
  }

  ae.data <- try(ArrayExpress::ArrayExpress(id), silent = TRUE)
  if (inherits(ae.data, "try-error")) {
    message("[pgx.getArrayExpress.data] Could not retrieve ", id, " from ArrayExpress\n")
    return(NULL)
  }

  counts <- Biobase::exprs(ae.data)
  meta <- Biobase::pData(ae.data)

  features <- rownames(Biobase::fData(ae.data))
  if (length(features) == nrow(counts)) rownames(counts) <- features

  ## has.fdata <- !is.null(Biobase::fData(eset)) && NCOL(Biobase::fData(eset)) > 0
  ## if (has.fdata) {
  ##   fdata <- Biobase::fData(eset)
  ## } else {
  ##   gpl.annot <- GEOquery::getGEO(eset@annotation)
  ##   fdata <- GEOquery::Table(gpl.annot)
  ## }
  ## if ("ID" %in% colnames(fdata)) {
  ##   cm <- intersect(rownames(ex), as.character(fdata$ID))
  ##   jj <- match(cm, as.character(fdata$ID))
  ##   fdata <- fdata[jj, , drop = FALSE]
  ##   rownames(fdata) <- cm
  ##   fdata <- fdata[, colnames(fdata) != "ID"]
  ##   ex <- ex[cm, , drop = FALSE]
  ## } else {
  ##   cm <- intersect(rownames(ex), rownames(fdata))
  ##   fdata <- fdata[cm, , drop = FALSE]
  ##   ex <- ex[cm, , drop = FALSE]
  ## }

  ## perform linear transformation if appropriate
  qq <- c(0., 0.25, 0.5, 0.75, 0.99, 1.0)
  qx <- as.numeric(stats::quantile(counts, qq, na.rm = T))
  is.count <- (qx[5] > 100) || (qx[6] - qx[1] > 50 && qx[2] > 0) ||
    (qx[2] > 0 && qx[2] < 1 && qx[4] > 1 && qx[4] < 2)
  if (!is.count) counts <- 2**counts

  LL <- list(expr = counts, samples = meta, source = "ArrayExpress")
  rm(ae.data, counts, features, meta, qq, qx)

  return(LL)
}
