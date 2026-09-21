#' Read Spectronaut output abundance global proteome data file.
#' @param file Path to Spectronaut output abundance file.
#' @return abundance data matrix (features on rows; samples on columns)
#' @export
read_spectronaut <- function(file) {
  msg <- function(...) message("[playbase.ingest::read_spectronaut] ", ...)

  counts <- suppressMessages(suppressWarnings(
    try(read_counts(file), silent = TRUE)
  ))
  if (inherits(counts, "try-error")) {
    counts <- try(read.csv(file, sep = "\t"), silent = TRUE)
    if (inherits(counts, "try-error")) {
      msg("FATAL: could not read abundance file")
      return(NULL)
    }
  }

  counts <- counts[!is.na(rownames(counts)) & rownames(counts) != "", , drop = FALSE]
  if (!any(is.na(suppressWarnings(as.numeric(rownames(counts)))))) {
    rownames(counts) <- paste0("Feature", 1:nrow(counts))
  }
  msg("Initial matrix size: ", nrow(counts), " x ", ncol(counts))

  pg.cols <- grep("^PG\\.", colnames(counts), value = TRUE)
  if (length(pg.cols) > 0) msg("pg.cols = ", paste(pg.cols, collapse = "; "))

  ## Columns required: PG.ProteinGroups; PG.Genes; PG.Quantity;
  pg.idx <- grep("proteingroups|accession", pg.cols, ignore.case = TRUE)
  genes.idx <- grep("\\.genes", pg.cols, ignore.case = TRUE)
  if (length(pg.idx) == 0 & length(genes.idx) == 0) {
    msg("FATAL: No protein groups, accession or genes found. Exiting")
    return(NULL)
  }

  protgr <- genes <- NULL
  if (length(pg.idx) > 0) protgr <- counts[, pg.cols[pg.idx], drop = FALSE]
  if (length(genes.idx) > 0) genes <- counts[, pg.cols[genes.idx], drop = FALSE]

  ## Remove contaminants
  if (!is.null(protgr)) {
    is.contam <- rowSums(sapply(protgr, grepl, pattern = "Cont_", ignore.case = TRUE)) > 0
    if (any(is.contam)) {
      msg("Identified ", sum(is.contam), " contaminants. Removing...")
      counts <- counts[!is.contam, , drop = FALSE]
      protgr <- protgr[!is.contam, , drop = FALSE]
      if (!is.null(genes)) genes <- genes[!is.contam, , drop = FALSE]
    }
  }

  ## Abundance data
  quant.idx <- grep("\\.quantity", colnames(counts), ignore.case = TRUE)
  if (length(quant.idx) == 0) {
    msg("FATAL: No abundances data found. Exiting")
    return(NULL)
  }
  msg(length(quant.idx), " protein abundances columns found...")
  counts <- counts[, quant.idx, drop = FALSE]

  ## Combine annot and abundance
  counts <- cbind(protgr, genes, counts)
  hh <- grep(".proteingroups|accession", colnames(counts), ignore.case = TRUE)
  if (length(hh) > 0) {
    rownames(counts) <- make.unique(counts[, hh[1]])
    counts <- counts[, -hh[1], drop = FALSE]
  }

  ## Clean colnames
  colnames(counts) <- gsub("^PG\\.|.PG.Quantity", "", colnames(counts))
  colnames(counts) <- gsub("^X\\.[0-9]+|\\.{2,}", "", colnames(counts))
  colnames(counts) <- gsub("\\.d$", "", colnames(counts))
  colnames(counts) <- gsub("-", "_", gsub("\\.", "_", colnames(counts)))

  msg("Completed. Final matrix size: ", nrow(counts), " x ", ncol(counts))
  rm(pg.cols, protgr, genes)
  gc()

  return(counts)
}

#' Read Spectronaut output hPTM abundance data file.
#' @param file Path to Spectronaut output hPTM abundance file.
#' @param use_ptm_norm Boolean. Use PTM abudances normalized to the global proteome. Default TRUE.
#' @return Abundance data matrix + annotation (features on rows; samples on columns)
#' @export
read_spectronaut_hPTM <- function(file, use_ptm_norm = TRUE) {

  msg <- function(...) message("[playbase.ingest::read_spectronaut_hPTM] ", ...)

  counts <- try(read.csv(file, sep = "\t"), silent = TRUE)
  if (inherits(counts, "try-error")) {
    counts <- try(data.table::fread(file, data.table = FALSE), silent = TRUE)
    if (inherits(counts, "try-error")) {
      msg("FATAL: could not read abundance file")
      return(NULL)
    }
  }
  
  counts <- counts[!is.na(rownames(counts)) & rownames(counts) != "", , drop = FALSE]
  if (!any(is.na(suppressWarnings(as.numeric(rownames(counts)))))) {
    rownames(counts) <- paste0("Feature", 1:nrow(counts))
  }
  msg("Initial matrix size: ", nrow(counts), " x ", ncol(counts))

  ## PTM annotation
  ss <- "proteinid|^pg\\.|key|modification|position|location|sitelocation|siteAA"
  ss <- paste0(ss, "|multiplicity|flankingregion")
  ann_idx <- unique(c(grep(ss, colnames(counts), ignore.case = TRUE)))
  if (length(ann_idx) == 0) {
    msg("FATAL: No PTM annotation columns found. Exiting")
    return(NULL)
  }
  ann <- counts[, ann_idx, drop = FALSE]
  counts <- counts[, -ann_idx, drop = FALSE]

  ## Retain only PTM phospho STY & ""
  mod_idx <- grep("modification", colnames(ann), ignore.case = TRUE)[1]
  if (!is.na(mod_idx)) {
    ptm_data <- as.character(ann[, mod_idx])
    keep <- grepl("sty", ptm_data, ignore.case = TRUE) | ptm_data == "" | is.na(ptm_data)
    ann <- ann[keep, , drop = FALSE]
    counts <- counts[keep, , drop = FALSE]
  }

  ## Feature names
  pid_idx <- grep("proteinid", colnames(ann), ignore.case = TRUE)[1]
  aa_idx <- grep("siteAA", colnames(ann), ignore.case = TRUE)[1]
  pos_idx <- grep("position|location|sitelocation", colnames(ann), ignore.case = TRUE)[1]
  ff <- rownames(ann)
  if (!is.na(pid_idx)) {
    ff <- as.character(ann[, pid_idx])
    jj <- which(is.na(ff) | ff == "")
    if (length(jj) > 0) ff[jj] <- rownames(ann)[jj]
    if (!is.na(aa_idx)) {
      jj <- which(!is.na(ann[, aa_idx]) & ann[, aa_idx] != "")
      ff[jj] <- paste0(ff[jj], "_", ann[jj, aa_idx])
      if (!is.na(pos_idx)) {
        jj <- which(!is.na(ann[, pos_idx]) & ann[, pos_idx] != "")
        ff[jj] <- paste0(ff[jj], ann[jj, pos_idx])
        ann[, pos_idx] <- as.character(ann[, pos_idx]) ## else interpreted as 'abundance'
      }
    }
  }
  rownames(ann) <- rownames(counts) <- make.unique(ff)
  
  mult_idx <- grep("multiplicity", colnames(ann), ignore.case = TRUE)[1]
  if (!is.na(mult_idx)) ann[, mult_idx] <- as.character(ann[, mult_idx]) ## else interpreted as 'abundance'
  
  ## Remove contaminants
  pid <- grep("proteinid", colnames(ann), ignore.case = TRUE)[1]
  if (!is.na(pid)) {
    is.contam <- grep("cont_", ann[, pid], ignore.case = TRUE)
    if (any(is.contam)) {
      msg("Identified ", length(is.contam), " contaminants. Removing...")
      ann <- ann[-is.contam, , drop = FALSE]
      counts <- counts[rownames(ann), , drop = FALSE]
    }
  }
  
  ## PTM.Quantity & PTM.QuantityPerProtein abundances (both are Spectronaut outputs).
  ## PTM.Quantity: PTM abundance not-normalized to the global proteome.
  ## PTM.QuantityPerProtein: PTM abundance normalized to the global proteome.
  if (use_ptm_norm) {
    ptm_qty <- grep("ptm\\..*quantityperprotein$", colnames(counts), ignore.case = TRUE)
    if (length(ptm_qty) == 0) {
      msg("No PTM.QuantityPerProtein columns found. Checking for PTM.Quantity columns...")
      ptm_qty <- grep("ptm\\..*quantity$", colnames(counts), ignore.case = TRUE)
      if (length(ptm_qty) == 0) {
        msg("FATAL: No PTM.Quantity columns found either. Exiting")
        return(NULL)
      }
    }
  } else {
    ptm_qty <- grep("ptm\\..*quantity$", colnames(counts), ignore.case = TRUE)
    if (length(ptm_qty) == 0) {
      msg("FATAL: No PTM quantity columns found. Exiting")
      return(NULL)
    }
  }

  counts <- counts[, ptm_qty, drop = FALSE]
  
  ## Spectronaut seems to write "Filtered" for some values.
  ## For now i coerce to numeric so they become NA. Real NaN/NA are preserved
  ## Subject to verification with client.
  counts[] <- lapply(counts, function(x) suppressWarnings(as.numeric(x)))  
  res <- cbind(ann, as.matrix(counts))  
  
  ## Clean colnames
  ss <- "^PG\\.|[.]PTM[.]Quantity|[.]PTM[.]QuantityPerProtein|PTM[.]|[.]PTM[.]|[.]raw"
  colnames(res) <- gsub(ss, "", colnames(res))
  colnames(res) <- gsub("^X\\.[0-9]+|\\.{2,}", "", colnames(res))
  colnames(res) <- gsub("\\.d$", "", colnames(res))

  msg("Completed. Final matrix size: ", nrow(res), " x ", ncol(res))
  rm(counts, ann)
  
  return(res)

}
