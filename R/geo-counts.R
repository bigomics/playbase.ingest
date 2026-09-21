#' @title Download counts for a GEO series
#' @description Tries each source in turn until one returns data: ARCHS4
#' (when \code{archs.h5} is given), GEOquery, recount, then ArrayExpress.
#' @param accession GEO accession ID.
#' @param archs.h5 Path to the ARCHS4 HDF5 file, or NULL to skip ARCHS4.
#' @return List with \code{expr} (counts), \code{samples} (metadata or NULL)
#'   and \code{source}; NULL if no source had the series.
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

#' @title Read GEO series counts from an ARCHS4 HDF5 file
#' @description Finds the samples of a GEO series in the ARCHS4 HDF5 file and
#' returns their count matrix. Log2 data is converted to linear, and
#' duplicated genes are summed on the linear scale.
#' @param id GEO accession ID.
#' @param h5.file Path to the ARCHS4 HDF5 file.
#' @return Counts matrix (genes x samples), or NULL if the series is not in
#'   the file.
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

#' @title Download GEO series counts from recount
#' @description Searches recount for the GEO series, downloads its
#' RangedSummarizedExperiment and returns the gene counts. Log2 data is
#' converted to linear, and duplicated genes are summed on the linear scale.
#' See the recount-quickstart vignette.
#' @param accession GEO accession ID.
#' @return Counts matrix (genes x samples), or NULL if recount has no match.
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

#' @title Download data from ArrayExpress
#' @description Downloads counts and sample metadata for an accession with the
#' ArrayExpress package. Log2 data is converted to linear.
#' @param accession ArrayExpress (or GEO) accession ID.
#' @return List with \code{expr}, \code{samples} and \code{source}, or NULL.
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
