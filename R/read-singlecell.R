## This file is part of the Omics Playground project.
## Copyright (c) 2018-2026 BigOmics Analytics SA. All rights reserved.


#' @title Convert Seurat to PGX
#' @param obj Seurat object to convert
#' @return List with the pgx slots a Seurat object carries: counts, X,
#'   samples, and any tsne/umap embeddings and seurat_clusters.
#' @description Converts a Seurat single-cell RNA-seq object into a pgx list.
#' @details Only extracts what the Seurat object already holds. Gene
#' annotation and sample clustering are playbase steps (e.g.
#' \code{playbase::pgx.clusterSamples}); run them on the result if needed.
#' @export
seurat2pgx <- function(obj) {
  message("[seurat2pgx] creating PGX object...")

  pgx <- list()
  pgx$name <- "SeuratProject"
  pgx$description <- "Seurat object converted using seurat2pgx"
  pgx$date <- Sys.Date()
  pgx$datatype <- "scRNA-seq"
  pgx$counts <- Seurat::GetAssayData(obj, assay = "RNA", layer = "counts")
  pgx$X <- Seurat::GetAssayData(obj, assay = "RNA", layer = "data")
  pgx$samples <- obj@meta.data

  ## copy clustering from Seurat
  if ("tsne" %in% names(obj@reductions)) {
    pos <- obj@reductions[["tsne"]]@cell.embeddings[, 1:2]
    pgx$cluster$pos[["tsne2d"]] <- pos
    pgx$tsne2d <- pos
    pgx$tsne3d <- cbind(pos, 0)
  }
  if ("umap" %in% names(obj@reductions)) {
    pgx$cluster$pos[["umap2d"]] <- obj@reductions[["umap"]]@cell.embeddings[, 1:2]
  }
  if ("seurat_clusters" %in% colnames(obj@meta.data)) {
    pgx$samples$cluster <- obj@meta.data[, "seurat_clusters"]
  }

  return(pgx)
}

#' @export
pgx.read_singlecell_counts <- function(filename) {
  counts <- NULL

  if (grepl("[.]csv$", filename)) {
    counts <- as.matrix(data.table::fread(filename, header = TRUE), rownames = 1)
  }

  if (grepl("[.]mtx$", filename)) {
    dir <- dirname(filename)
    barcode.file <- file.path(dir, "barcodes.tsv")
    genes.file <- file.path(dir, "genes.tsv")
    if (!file.exists(filename)) stop("could not find counts matrix: ", filename)
    if (!file.exists(barcode.file)) stop("could not find barcode file: ", barcode.file)
    if (!file.exists(genes.file)) stop("could not find genes file: ", genes.file)
    counts <- Matrix::readMM(filename)
    bc <- read.csv(barcode.file, header = FALSE, sep = "\t")
    gn <- read.csv(genes.file, header = FALSE, sep = "\t")
    rownames(counts) <- gn[, 2] ## gene names
    colnames(counts) <- bc[, 1]
  }

  if (grepl("[.]h5$", filename)) {
    counts <- Seurat::Read10X_h5(filename, use.names = TRUE, unique.features = TRUE)
  }

  counts
}

#' Read 10X Cell Ranger Software output (version V3 onwards).
#' @param file .tar.gz or .zip compressed directory
#' @return Count gene expression data matrix (sparse dgCMatrix)
#' @export
read_cellranger_output <- function(file) {
  if (!requireNamespace("Seurat", quietly = TRUE)) {
    stop("read_cellranger_output() requires the 'Seurat' package (Suggests); install it to read 10X Cell Ranger output.")
  }
  msg <- function(...) message("[playbase.ingest::read_cellranger_output] ", ...)

  msg("Reading 10X Cell Ranger output...")
  tmp <- tempfile()
  dir.create(tmp)
  on.exit(unlink(tmp, recursive = TRUE), add = TRUE) # deletes tmp always.

  if (grepl("\\.tar\\.gz$|\\.gz$", file)) {
    msg(".tar or .gz compressed file detected...")
    utils::untar(file, exdir = tmp)
    dir <- tmp
  } else if (grepl("\\.zip$", file)) {
    msg(".zip compressed file detected...")
    utils::unzip(file, exdir = tmp)
    dir <- tmp
  }

  ff1 <- c("barcodes.tsv.gz", "features.tsv.gz", "matrix.mtx.gz")
  ff2 <- c("barcodes.tsv", "genes.tsv", "matrix.mtx")
  mex_dir <- NULL
  Data <- list.dirs(dir, recursive = TRUE)
  for (i in 1:length(Data)) {
    files <- list.files(Data[i])
    has_mex <- all(ff1 %in% files) || all(ff2 %in% files)
    if (has_mex) {
      mex_dir <- Data[i]
      break
    }
  }

  if (is.null(mex_dir)) {
    msg("Could not find MEX directory (barcodes/features/matrix files) in: ", dir)
    return(NULL)
  }

  require(Seurat) # do not remove
  counts <- Seurat::Read10X(data.dir = mex_dir)
  if (is.list(counts)) counts <- counts[["Gene Expression"]]

  msg("Completed. Expression matrix: ", nrow(counts), " x ", ncol(counts), ".\n")
  gc()
  return(counts)
}
