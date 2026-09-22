##
## This file is part of the Omics Playground project.
## Copyright (c) 2018-2026 BigOmics Analytics SA. All rights reserved.
##


## -------------------------------------------------------------------------------------
## Query GEO
## -------------------------------------------------------------------------------------

#' @title Download a GEO series
#' @description Downloads counts and sample metadata for a GEO series.
#' @param accession GEO accession ID.
#' @param archs.h5 Path to an ARCHS4 HDF5 file. NULL (default) skips ARCHS4.
#' @param get.info Logical. Also download the series description with
#'   \code{pgx.getGEOexperimentInfo()}.
#' @return List with \code{counts}, \code{samples}, \code{info} and
#'   \code{source}.
#' @details Downloads GEO accession ID data. It tries ARCHS4 only when
#' \code{archs.h5} is given, then GEO, then recount. Counts and sample
#' matrices are aligned.
#' @export
pgx.getGEOseries <- function(accession,
                             archs.h5 = NULL,
                             get.info = TRUE) {
  id <- accession
  is.valid.id <- is.GEO.id.valid(id)
  if (!is.valid.id) stop("[pgx.getGEOseries] FATAL: ID is invalid. Exiting.")
  id <- as.character(id)

  meta <- NULL
  geo <- pgx.getGEOcounts(id, archs.h5 = archs.h5)
  source <- geo[["source"]]
  counts <- geo[["expr"]]
  meta <- geo[["samples"]]
  if (is.null(counts)) {
    message("[pgx.getGEOseries] WARNING:", id, " not found in GEO, recount, ArrayExpress.\n")
    return(NULL)
  }

  if (is.null(meta)) meta <- pgx.getGEOmetadata(id)
  if (!is.null(meta)) {
    nn <- apply(meta, 2, function(x) length(unique(x)))
    ex.vars <- names(nn)[which(nn == 1)]
    if (length(ex.vars) > 0) {
      meta <- meta[, !colnames(meta) %in% ex.vars, drop = FALSE]
    }
    hh <- grep("characteristics_ch", colnames(meta))
    if (any(hh)) meta <- meta[, -hh, drop = FALSE]
  } else {
    message("[pgx.getGEOseries] WARNING: Metadata not retrieved.")
  }

  ## conform matrices
  if (!is.null(meta)) {
    samples <- intersect(rownames(meta), colnames(counts))
    if (length(samples) == 0) {
      i <- 1
      for (i in 1:ncol(counts)) {
        hh1 <- grep(colnames(counts)[i], meta)
        if (length(hh1) == 0) next
        hh2 <- grep(colnames(counts)[i], meta[, hh1])
        if (length(hh2) == 0) next
        colnames(counts)[i] <- unique(rownames(meta)[hh2])[1]
      }
      samples <- intersect(rownames(meta), colnames(counts))
      if (length(samples) == 0) {
        message("[pgx.getGEOseries] WARNING: No shared samples between counts and metadata.")
      } else {
        meta <- meta[samples, , drop = FALSE]
        counts <- counts[, samples, drop = FALSE]
      }
    } else {
      meta <- meta[samples, , drop = FALSE]
      counts <- counts[, samples, drop = FALSE]
    }
  }

  ## get ID experiment info
  info <- NULL
  if (get.info) info <- pgx.getGEOexperimentInfo(id)

  ## get categorical phenotypes
  # meta1 <- apply(meta, 2, trimsame)
  # rownames(meta1) <- rownames(meta)
  # sampleinfo <- pgx.discretizePhenotypeMatrix(meta1, min.ncat = 2,
  # max.ncat = 20, remove.dup = TRUE)
  # sampleinfo <- data.frame(sampleinfo, stringsAsFactors = FALSE, check.names = FALSE)
  ## automagically create contrast matrix
  # contrasts <- NULL
  # if (NCOL(sampleinfo) > 0) {
  #  mingrp <- 3
  #  slen <- 15
  #  ref <- NA
  #  ct <- pgx.makeAutoContrasts(sampleinfo, mingrp = 3, slen = 20, ref = NA)
  #  if (is.null(ct)) {
  #    ct <- pgx.makeAutoContrasts(sampleinfo, mingrp = 2, slen = 20, ref = NA)
  #  }
  #  if (!is.null(ct$exp.matrix)) {
  #    contrasts <- ct$exp.matrix
  #  } else {
  #    contrasts <- ct$design %*% ct$contr.matrix
  #  }
  # }

  LL <- list(counts = counts, samples = meta, info = info, source = source)
  return(LL)
}

## -------------------------------------------------------------------------------------
## HELPER functions
## -------------------------------------------------------------------------------------

#' @title Check a GEO accession ID
#' @param accession GEO accession ID.
#' @return TRUE if the ID is non-empty and contains both letters and digits.
#' @export
is.GEO.id.valid <- function(accession) {
  id <- accession
  is.valid <- TRUE
  if (is.null(id) || id == "" || !grepl("[A-Za-z]", id) || !grepl("[0-9]", id)) {
    is.valid <- FALSE
  }
  return(is.valid)
}
