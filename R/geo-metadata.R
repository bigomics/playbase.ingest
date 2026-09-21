#' @title Download sample metadata for a GEO series
#' @description Tries the per-sample GSM records first, then the GSE Series
#' Matrix.
#' @param accession GEO accession ID.
#' @return Data frame of sample metadata (samples in rows), or NULL.
#' @export
pgx.getGEOmetadata <- function(accession) {
  id <- accession
  is.valid.id <- is.GEO.id.valid(id)
  if (!is.valid.id) stop("[pgx.getGEOmetadata] FATAL: ID is invalid. Exiting.")
  id <- as.character(id)

  meta <- pgx.getGEOmetadata.fromGSM(id) ## no GSEMatrix
  if (is.null(meta)) meta <- pgx.getGEOmetadata.fromEset(id) ## with GSEMatrix

  ## Sometimes the phenotype is coded in the title string
  # if ("title" %in% colnames(meta) && NCOL(meta) == 0) {
  #  px <- title2pheno(meta$title, split = NULL, trim = TRUE, summarize = TRUE)
  #  if (!is.null(px) && NCOL(px) > 0 && is.null(meta))
  #    pheno <- px
  #  if (!is.null(px) && NCOL(px) > 0 && !is.null(pheno))
  #     pheno <- cbind(pheno, px)
  # }

  return(meta)
}

## -------------------------------------------------------------------------------------
## Query GEO metadata
## -------------------------------------------------------------------------------------

#' @title Download the description of a GEO series
#' @param id GEO accession ID.
#' @return The GSE header (a named list: title, summary, design, ...), or NULL.
#' @export
pgx.getGEOexperimentInfo <- function(id) {
  is.valid.id <- is.GEO.id.valid(id)
  if (!is.valid.id) stop("[pgx.getGEOexperimentInfo] FATAL: ID is invalid. Exiting.")
  id <- as.character(id)

  suppressMessages(
    gse <- try(GEOquery::getGEO(id, GSEMatrix = FALSE, getGPL = FALSE), silent = TRUE)
  )
  if (inherits(gse, "try-error")) {
    message("[pgx.getGEOexperimentInfo] Error: GEOquery::getGEO failed to get", id, ".\n")
    return(NULL)
  }

  return(gse@header) ## can be a big list!
}

#' @title Sample metadata from GEO GSM records
#' @description Builds the sample table from the individual GSM records of a
#' series, without the GSE Series Matrix files.
#' @param id GEO accession ID.
#' @return Data frame of sample metadata, or NULL.
#' @export
pgx.getGEOmetadata.fromGSM <- function(id) {
  is.valid.id <- is.GEO.id.valid(id)
  if (!is.valid.id) stop("[pgx.getGEOmetadata.fromGSM] FATAL: ID is invalid. Exiting.")
  id <- as.character(id)

  message("[pgx.getGEOmetadata.fromGSM] Attempt to download metadata without GSEMatrix...")
  suppressMessages(
    gse <- try(GEOquery::getGEO(id, GSEMatrix = FALSE, getGPL = FALSE), silent = TRUE)
  )
  if (inherits(gse, "try-error")) {
    message("[pgx.getGEOmetadata.fromGSM] Error: getGEO failed to retrieve metadata for ", id, "\n")
    return(NULL)
  }

  if (length(gse@gsms) == 0) {
    message("[pgx.getGEOmetadata.fromGSM] WARNING: no GSM information in object. Exiting. \n")
    return(NULL)
  }

  ## get metadata
  # summary <- gse@header$summary
  gsm.title <- sapply(gse@gsms, function(g) g@header$title)
  gsm.source <- sapply(gse@gsms, function(g) g@header$source_name_ch1)
  gsm.gpl <- sapply(gse@gsms, function(g) g@header$platform_id)
  gsm.samples <- gse@header$sample_id
  # geo.accession <- gse@header$geo_accession
  meta <- data.frame(
    # geo_accession = geo.accession,
    GPL = gsm.gpl,
    GSM = gsm.samples,
    title = gsm.title,
    source = gsm.source,
    stringsAsFactors = FALSE,
    check.names = FALSE
  )
  ch1_info <- lapply(gse@gsms, function(g) g@header$characteristics_ch1)
  cm <- intersect(names(ch1_info), gsm.samples)
  if (!is.null(ch1_info) && length(cm)) {
    gsm.samples <- gsm.samples[match(cm, gsm.samples)]
    ch1_info <- ch1_info[match(cm, names(ch1_info))]
    meta <- meta[match(cm, meta$GSM), ]
    ch1_info <- lapply(ch1_info, function(x) sub("^Clinical info: ", "", x))
    ch1_vars <- unique(unlist(lapply(ch1_info, function(x) trimws(sub("[:=].*", "", x)))))
    ch1_info <- lapply(ch1_info, function(x) {
      xvar <- trimws(sub("[:=].*", "", x))
      x <- trimws(sub(".*[:=] ", "", x))
      names(x) <- xvar
      x <- x[match(ch1_vars, names(x))]
      return(x)
    })
    ch1_info <- do.call(rbind, ch1_info)
    colnames(ch1_info) <- ch1_vars
    meta <- data.frame(cbind(meta, ch1_info), stringsAsFactors = FALSE, check.names = FALSE)
    colnames(meta) <- gsub("[ ]", "_", colnames(meta))

    message("[pgx.getGEOmetadata.fromGSM] Success!")
    return(meta)
  } else {
    message("[pgx.getGEOmetadata.fromGSM] WARNING: no shared samples between GSM & ch1_info. Exiting.\n")
    return(NULL)
  }

  # is.underscored <- length(gsm.title) && all(grepl("_", gsm.title))
  # title_info <- NULL
  ## NEED RETHINK!!!!!!!!!!!!!!!!!!
  # if (FALSE && is.underscored) {
  #  title2 <- trimws(gsm.title)
  #  title_info <- eset.parsePhenoFromTitle(title2, split = "_")
  # }
  # if (!is.null(title_info)) sample_info <- cbind(sample_info, title_info)
}

#' @title Sample metadata from the GEO Series Matrix
#' @description Downloads the GSE Series Matrix files and extracts the sample
#' metadata from their ExpressionSets.
#' @param id GEO accession ID.
#' @return Data frame of sample metadata, or NULL.
#' @export
pgx.getGEOmetadata.fromEset <- function(id) {
  is.valid.id <- is.GEO.id.valid(id)
  if (!is.valid.id) stop("[pgx.getGEOmetadata.fromEset] FATAL: ID is invalid. Exiting.")
  id <- as.character(id)

  message("[pgx.getGEOmetadata.fromEset] Attempt to download metadata with GSEMatrix...")
  suppressMessages(
    gse <- try(GEOquery::getGEO(id, GSEMatrix = TRUE, getGPL = FALSE), silent = TRUE)
  )
  if (inherits(gse, "try-error")) {
    message("[pgx.getGEOmetadata.fromEset] Error: getGEO failed to retrieve metadata for ", id, "\n")
    return(NULL)
  }

  nsamples <- sapply(gse, function(s) nrow(Biobase::pData(Biobase::phenoData(s))))
  gse <- gse[which(nsamples >= 3)]
  meta.list <- lapply(gse, function(x) pgx.getGEOmetadata.fromEset.helper(x))
  meta <- do.call(rbind, meta.list)
  meta <- data.frame(meta, stringsAsFactors = FALSE, check.names = FALSE)
  colnames(meta) <- gsub("[ ]", "_", colnames(meta))
  message("[pgx.getGEOmetadata.fromEset] Success!")
  rm(meta.list)

  return(meta)
}

#' @title Sample metadata from one ExpressionSet
#' @param eset ExpressionSet from GEOquery.
#' @return Data frame of sample metadata.
#' @export
pgx.getGEOmetadata.fromEset.helper <- function(eset) {
  if (!class(eset) %in% "ExpressionSet") {
    message("[pgx.getGEOmetadata.fromEset.helper] Error: eset must be of class 'ExpressionSet'")
    return(NULL)
  }

  ## get metadata
  meta0 <- Biobase::pData(Biobase::phenoData(eset))
  gsm.gpl <- as.character(meta0$platform_id)
  gsm.samples <- as.character(meta0$geo_accession)
  gsm.title <- as.character(meta0$title)
  gsm.source <- as.character(meta0$source_name_ch1)
  meta <- data.frame(
    GPL = gsm.gpl,
    GSM = gsm.samples,
    title = gsm.title,
    source = gsm.source,
    stringsAsFactors = FALSE,
    check.names = FALSE
  )

  ch1_info <- NULL
  ch1_sel <- grepl("characteristics_ch1", colnames(meta0))
  if (any(ch1_sel)) {
    ch1_info <- meta0[, ch1_sel, drop = FALSE]
    colnames(ch1_info) <- paste0("characteristics_", 1:ncol(ch1_info))
    ch1_info <- apply(ch1_info, 2, function(x) sub("^Clinical info: ", "", x))
    ch_vars <- apply(ch1_info, 2, function(x) {
      return(unique(sub("[:=].*", "", x)))
    })
    colnames(ch1_info) <- unname(ch_vars)
    ch1_info <- apply(ch1_info, 2, function(x) trimws(sub(".*[:=] ", "", x)))
    meta <- cbind(meta, ch1_info)
    meta <- data.frame(meta, stringsAsFactors = FALSE, check.names = FALSE)
  }
  rm(meta0)

  ## Base sample_info from characteristics (ch1) column
  # ch1_info <- eset.getCH1(eset)
  # We can get extra information from title
  # is.underscored <- length(gsm.title) && all(grepl("_", gsm.title))
  # title_info <- NULL
  # if (FALSE && is.underscored) {
  #  title2 <- trimws(gsm.title)
  #  title_info <- eset.parsePhenoFromTitle(title2, split = "_")
  # }
  ## All sample_info: from characterisctis_ch1 and title
  # if (!is.null(ch1_info)) sample_info <- cbind(sample_info, ch1_info)
  # if (!is.null(title_info)) sample_info <- cbind(sample_info, title_info)

  return(meta)
}
