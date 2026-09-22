#' @title Download GEO series counts with GEOquery
#' @description Downloads the series matrix (or RNA-seq counts / supplementary
#' files), platform metadata and probe annotation from GEO, and maps probes to
#' gene symbols. Log2 data is converted to linear, and duplicated genes are
#' summed on the linear scale.
#' @param accession GEO accession ID.
#' @return List with \code{expr} (counts) and \code{meta} (sample metadata),
#'   or NULL.
#' @export
pgx.getGEOcounts.GEOquery <- function(accession) {
  id <- accession
  is.valid.id <- is.GEO.id.valid(id)
  if (!is.valid.id) stop("[pgx.getGEOcounts.GEOquery] FATAL: ID is invalid. Exiting.")
  id <- as.character(id)

  meta <- NULL

  gse <- try(GEOquery::getGEO(GEO = id, GSEMatrix = TRUE, getGPL = TRUE), silent = TRUE)
  has.expr <- FALSE
  if (!inherits(gse, "try-error")) {
    has.expr <- sapply(gse, function(x) nrow(Biobase::exprs(x)) > 0)
  }

  if (inherits(gse, "try-error") | !any(has.expr)) {
    gse <- try(GEOquery::getRNASeqData(accession = id), silent = TRUE)
    if (inherits(gse, "try-error")) {
      message("[pgx.getGEOcounts.GEOquery] getGEO failed to retrieve ", id, "\n")
      return(NULL)
    }

    if (class(gse) %in% "SummarizedExperiment") {
      counts <- try(SummarizedExperiment::assay(gse), silent = TRUE)
      meta <- try(as.data.frame(SummarizedExperiment::colData(gse)), silent = TRUE)
      if (inherits(counts, "try-error")) {
        message("[pgx.getGEOcounts.GEOquery] getGEO failed to retrieve ", id, "\n")
        return(NULL)
      }

      features <- try(as.data.frame(SummarizedExperiment::rowData(gse)), silent = TRUE)
      if (!inherits(features, "try-error")) {
        jj <- which(!is.na(rownames(features)) & rownames(features) != "")
        features <- features[jj, , drop = FALSE]
        cm <- intersect(rownames(features), rownames(counts))
        if (length(cm) > 0) {
          features <- features[cm, , drop = FALSE]
          counts <- counts[cm, , drop = FALSE]
        }
        hh1 <- grepl("symbol", colnames(features), ignore.case = TRUE)
        hh2 <- grepl("Ensembl", colnames(features), ignore.case = TRUE)
        ff <- NULL
        if (any(hh1)) ff <- features[, which(hh1)[1]]
        if (!any(hh1) && any(hh2)) ff <- features[, which(hh2)[1]]
        if (!is.null(ff)) {
          jj <- which(is.na(ff) | ff == "")
          if (length(jj) > 0) ff[jj] <- rownames(features)[jj]
          rownames(counts) <- as.character(ff)
        }
        rm(ff, features)
        gc()
      }

      LL <- list(expr = counts, meta = meta)
      rm(counts, meta)
      gc()
      return(LL)
    }
  }

  supp_file <- NULL
  has.expr <- sapply(gse, function(x) nrow(Biobase::exprs(x)) > 0)
  if (any(has.expr)) {
    gse <- gse[which(has.expr)]
  } else {
    message("[pgx.getGEOcounts.GEOquery] WARNING: no data found in ", id, " from GEO.\n")
    supp_file <- sapply(gse, function(g) g@experimentData@other$supplementary_file)
    supp_file <- unname(supp_file[[1]])
    if (!is.null(supp_file)) {
      sfiles <- strsplit(supp_file, "\n")[[1]]
      csvfile <- head(grep(".csv", sfiles, fixed = TRUE), 1)
      if (length(csvfile)) {
        ext <- ".csv"
        if (tools::file_ext(sfiles[csvfile]) == "gz") ext <- ".csv.gz"
        destfile <- tempfile(fileext = ext)
        dd <- try(utils::download.file(url = sfiles[csvfile], destfile = destfile), silent = TRUE)
        if (inherits(dd, "try-error")) {
          message("[pgx.getGEOcounts.GEOquery] 1st attempt in downloading supp file failed. Trying again.\n")
          dd <- try(utils::download.file(url = sfiles[csvfile], destfile = destfile), silent = TRUE)
          if (inherits(dd, "try-error")) {
            message("[pgx.getGEOcounts.GEOquery] Error in downloading supp file: ", sfiles[csvfile], ". Exiting \n")
            return(NULL)
          }
        }
        file <- destfile
        if (ext == ".csv.gz") {
          R.utils::gunzip(destfile, remove = TRUE)
          file <- gsub(".gz", "", destfile)
        }
        counts <- read_counts(file)
        base::file.remove(file)
        return(counts)
      }
    } else {
      message("[pgx.getGEOcounts.GEOquery] getGEO failed to retrieve ", id, "\n")
      return(NULL)
    }
  }

  ## select preferred platform is multiple exists
  k <- 1
  expr.list <- list()
  for (k in 1:length(gse)) {
    eset <- gse[[k]]
    ex <- Biobase::exprs(eset)
    if (ncol(ex) <= 3) {
      message("[pgx.getGEOcounts.GEOquery] WARNING: ", id, " contains <= 3 samples. Skipping.\n")
      next()
    }

    ## perform linear transformation (unlog) if required
    qx <- as.numeric(stats::quantile(ex, c(0., 0.25, 0.5, 0.75, 0.99, 1.0), na.rm = T))
    is.count <- (qx[5] > 100) || (qx[6] - qx[1] > 50 && qx[2] > 0) ||
      (qx[2] > 0 && qx[2] < 1 && qx[4] > 1 && qx[4] < 2) ## from GEO2R script
    if (!is.count) ex <- 2**ex

    ## featuredata
    has.fdata <- !is.null(Biobase::fData(eset)) && NCOL(Biobase::fData(eset)) > 0
    if (has.fdata) {
      fdata <- Biobase::fData(eset)
    } else {
      gpl.annot <- GEOquery::getGEO(eset@annotation)
      fdata <- GEOquery::Table(gpl.annot)
    }
    if ("ID" %in% colnames(fdata)) {
      cm <- intersect(rownames(ex), as.character(fdata$ID))
      jj <- match(cm, as.character(fdata$ID))
      fdata <- fdata[jj, , drop = FALSE]
      rownames(fdata) <- cm
      fdata <- fdata[, colnames(fdata) != "ID"]
      ex <- ex[cm, , drop = FALSE]
    } else {
      cm <- intersect(rownames(ex), rownames(fdata))
      fdata <- fdata[cm, , drop = FALSE]
      ex <- ex[cm, , drop = FALSE]
    }

    ## get symbol from featuredata; clean and sum linear intensities
    fsymbol <- pgx.getSymbolFromFeatureData(fdata)
    jj <- which(!is.na(fsymbol) & fsymbol != "")
    ex <- ex[jj, ]
    fsymbol <- fsymbol[jj]
    fsymbol <- gsub(" /// ", ";", fsymbol)
    ex2 <- tapply(1:nrow(ex), fsymbol, function(ii) {
      Matrix::colSums(ex[ii, , drop = FALSE], na.rm = TRUE) ## not log!!
    })
    ex2 <- do.call(rbind, ex2)
    expr.list[[names(gse)[k]]] <- ex2
  }

  if (length(expr.list) == 0) {
    return(NULL)
  }

  if (length(expr.list) > 1) {
    ## merge/join all expressions
    probes <- sort(unique(unlist(lapply(expr.list, rownames))))
    samples <- sort(unique(unlist(lapply(expr.list, colnames))))
    expr.list2 <- lapply(expr.list, function(x) {
      x[match(probes, rownames(x)), match(samples, colnames(x))]
    })
    expr.list2 <- lapply(expr.list2, function(x) {
      x[is.na(x)] <- 0
      x
    })
    expr <- Reduce("+", expr.list2)
    colnames(expr) <- samples
    rownames(expr) <- probes
  } else {
    expr <- expr.list[[1]]
  }

  LL <- list(expr = expr, meta = meta)
  rm(expr, meta)
  gc()
  return(LL)
  # return(expr) ## linear intensities
}

#' @title Extract gene symbols from GEO feature data
#' @param fdata The featureData table from a GEOquery GEO dataset object.
#' @description Extracts official gene symbols from feature datatable of a GEO dataset downloaded with GEOquery.
#' It first looks for a column containing gene symbols by matching against the org.Hs.egSYMBOL database.
#' If no direct symbol column is found, it looks for an ENTREZ identifier and maps to symbols using org.Hs.egSYMBOL.
#' Then it looks at REFSEQ identifiers. Then it looks at ENSEMBLE IDs. If no approach works, it returns NULL.
#' @return A character vector of gene symbols, or NULL if symbols could not be extracted.
#' @export
pgx.getSymbolFromFeatureData <- function(fdata) {
  symbol <- NULL

  ## SYMBOL column
  SYMBOL <- as.character(unlist(as.list(org.Hs.eg.db::org.Hs.egSYMBOL)))
  symbol.col <- grep("symbol|gene|hugo", colnames(fdata), ignore.case = TRUE)
  if (any(symbol.col)) {
    ok.symbol <- apply(
      fdata[, symbol.col, drop = FALSE], 2,
      function(g) mean(toupper(g[!is.na(g)]) %in% SYMBOL)
    )
    if (any(ok.symbol > 0.5)) {
      k <- which.max(ok.symbol)
      symbol <- fdata[, symbol.col[k]]
      message("[pgx.getSymbolFromFeatureData] SYMBOL column found. Returning gene symbols...")
      return(symbol)
    }
  }

  ## ENTREZ column
  ENTREZ <- AnnotationDbi::keys(org.Hs.eg.db::org.Hs.egSYMBOL)
  entrez.col <- grep("entrez", colnames(fdata), ignore.case = TRUE)
  if (any(entrez.col)) {
    entrez.match <- apply(
      fdata[, entrez.col, drop = FALSE], 2,
      function(g) mean(g[!is.na(g)] %in% ENTREZ)
    )
    entrez.ok <- length(entrez.col) && entrez.match > 0.5
    if (entrez.ok) {
      k <- entrez.col[which.max(entrez.match)]
      probes <- as.character(fdata[, k])
      symbol <- AnnotationDbi::mapIds(org.Hs.eg.db::org.Hs.eg.db, probes, "SYMBOL", "ENTREZID")
      message("[pgx.getSymbolFromFeatureData] ENTREZ column found. Returning gene symbols...")
      return(symbol)
    }
  }

  ## REFSEQ column
  REFSEQ <- unlist(as.list(org.Hs.eg.db::org.Hs.egREFSEQ))
  refseq.col <- grep("refseq", colnames(fdata), ignore.case = TRUE)
  if (any(refseq.col)) {
    refseq.match <- apply(
      fdata[, refseq.col, drop = FALSE], 2,
      function(g) mean(sub("[.].*", "", g[!is.na(g)]) %in% REFSEQ)
    )
    refseq.ok <- length(refseq.col) && refseq.match > 0.5
    if (refseq.ok) {
      k <- refseq.col[which.max(refseq.match)]
      probes <- sub("[.].*", "", as.character(fdata[, k]))
      symbol <- AnnotationDbi::mapIds(org.Hs.eg.db::org.Hs.eg.db, probes, "SYMBOL", "REFSEQ")
      message("[pgx.getSymbolFromFeatureData] REFSEQ column found. Returning gene symbols...")
      return(symbol)
    }
  }

  ## EnsembleID column
  gene.column <- grep("gene|mrna|transcript", colnames(fdata), ignore.case = TRUE)
  has.ens <- apply(fdata[, gene.column, drop = FALSE], 2, function(s) mean(grepl("ENS", s)))
  if (any(has.ens > 0.3)) {
    ens.col <- ifelse(max(has.ens) > 0, names(which.max(has.ens)), NA)
    ens.ann <- lapply(fdata[, ens.col], function(a) trimws(strsplit(a, split = "//|///")[[1]]))
    ens.probes <- sapply(ens.ann, function(s) Matrix::head(grep("^ENS", s, value = TRUE), 1))
    ens.probes[sapply(ens.probes, length) == 0] <- NA
    ens.probes <- sub("[.].*", "", unlist(ens.probes)) ## drop version suffix
    symbol <- rep(NA_character_, length(ens.probes))
    ok <- !is.na(ens.probes)
    symbol[ok] <- AnnotationDbi::mapIds(org.Hs.eg.db::org.Hs.eg.db, ens.probes[ok], "SYMBOL", "ENSEMBL")
    message("[pgx.getSymbolFromFeatureData] ENSEMBLE ID column found. Returning gene symbols...")
    return(symbol)
  }

  message("WARNING:: could not parse symbol information from featureData!")
  return(NULL)
}
