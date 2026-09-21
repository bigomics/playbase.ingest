#' Read scRNA-seq counts matrix in h5 format.
#' Automatically 'infer' counts;features;cells from the h5 file structure.
#' Attempts multiple ways.
#' @export
read_h5_counts <- function(h5.file) {
  if (!requireNamespace("rhdf5", quietly = TRUE)) {
    stop("read_h5_counts() requires the 'rhdf5' package (Suggests); install it to read HDF5/h5ad files.")
  }
  message("[playbase.ingest::read_h5_counts] Reading h5 file: ", h5.file)
  df <- NULL

  FF <- tryCatch(
    {
      rhdf5::h5ls(h5.file)
    },
    error = function(w) {
      NULL
    }
  )

  if (!is.null(FF) && all(c("group", "name") %in% colnames(FF))) {
    h5.ems <- paste0(FF[, "group"], "/", FF[, "name"])
    LL <- lapply(h5.ems, function(ems) rhdf5::h5read(file = h5.file, name = ems))
    dims <- unlist(lapply(LL, function(x) length(dim(x))))
    ll <- lapply(LL, length)
    if (any(dims == 2)) df <- LL[[which(dims == 2)[1]]]
    if (!is.null(df)) {
      if (is.null(rownames(df)) && any(ll == nrow(df))) {
        rownames(df) <- as.character(LL[[which(ll == nrow(df))[1]]])
      }
      if (is.null(colnames(df)) && any(ll == ncol(df))) {
        colnames(df) <- as.character(LL[[which(ll == ncol(df))[1]]])
      }
    }
  }

  # AnnData h5ad: /X stored as sparse CSR (data + indices + indptr)
  if (is.null(df) && !is.null(FF)) {
    h5.paths <- paste0(FF[, "group"], "/", FF[, "name"])
    is_anndata <- all(c("/X/data", "/X/indices", "/X/indptr") %in% h5.paths)
    if (is_anndata) {
      df <- tryCatch(
        {
          x_data <- rhdf5::h5read(h5.file, "/X/data")
          x_indices <- rhdf5::h5read(h5.file, "/X/indices")
          x_indptr <- rhdf5::h5read(h5.file, "/X/indptr")
          ## rhdf5::h5read returns 1D arrays (with a dim attribute), not plain
          ## character vectors. Seurat's LogMap[[<-]] dispatch fails when dimnames
          ## are arrays rather than vectors. Strip with as.character().
          obs_names <- as.character(rhdf5::h5read(h5.file, "/obs/_index"))
          var_names <- as.character(rhdf5::h5read(h5.file, "/var/_index"))
          n_obs <- length(obs_names)
          n_vars <- length(var_names)
          ## CSR indptr length = n_rows + 1; some tools write X as vars x obs (transposed)
          n_csr_rows <- length(x_indptr) - 1L
          transposed <- (n_csr_rows == n_vars && n_csr_rows != n_obs)
          if (transposed) {
            row_names <- var_names
            col_names <- obs_names
            n_rows <- n_vars
            n_cols <- n_obs
          } else {
            row_names <- obs_names
            col_names <- var_names
            n_rows <- n_obs
            n_cols <- n_vars
          }
          message("[playbase.ingest::read_h5_counts] AnnData sparse CSR detected")
          message("[playbase.ingest::read_h5_counts] Size: ", n_obs, " cells; ", n_vars, " features")
          row_idx <- rep(seq_len(n_rows), diff(as.integer(x_indptr)))
          mat <- Matrix::sparseMatrix(
            i = row_idx,
            j = as.integer(x_indices) + 1L,
            x = as.numeric(x_data),
            dims = c(n_rows, n_cols),
            dimnames = list(row_names, col_names)
          )
          ## ensure final result is genes x cells
          counts_mat <- if (transposed) mat else Matrix::t(mat)

          ## read /obs cell metadata
          obs_meta <- tryCatch(
            {
              FF_obs <- FF[FF[, "group"] == "/obs", , drop = FALSE]
              col_names_obs <- FF_obs[, "name"]
              ## skip internal AnnData fields and sub-group entries
              skip <- c("_index", "__categories")
              is_subgroup <- col_names_obs %in% FF[FF[, "group"] != "/obs", "name"]
              top_cols <- col_names_obs[!col_names_obs %in% skip & !is_subgroup]
              ## also skip columns that have sub-paths (categorical stored as categories/codes)
              has_subpath <- vapply(top_cols, function(cn) {
                any(FF[, "group"] == paste0("/obs/", cn))
              }, logical(1))
              cat_cols <- top_cols[has_subpath]
              scal_cols <- top_cols[!has_subpath]

              meta <- data.frame(row.names = obs_names)
              for (cn in scal_cols) {
                v <- tryCatch(as.vector(rhdf5::h5read(h5.file, paste0("/obs/", cn))), error = function(e) NULL)
                if (!is.null(v) && length(v) == n_obs) meta[[cn]] <- v
              }
              for (cn in cat_cols) {
                cats <- tryCatch(as.character(rhdf5::h5read(h5.file, paste0("/obs/", cn, "/categories"))), error = function(e) NULL)
                codes <- tryCatch(as.integer(rhdf5::h5read(h5.file, paste0("/obs/", cn, "/codes"))), error = function(e) NULL)
                if (!is.null(cats) && !is.null(codes) && length(codes) == n_obs) {
                  ## AnnData encodes missing categoricals as -1; map to NA
                  idx <- codes + 1L
                  meta[[cn]] <- ifelse(idx > 0L, cats[ifelse(idx > 0L, idx, 1L)], NA_character_)
                }
              }
              if (ncol(meta) > 0) meta else NULL
            },
            error = function(e) NULL
          )

          list(counts = counts_mat, samples = obs_meta)
        },
        error = function(w) {
          message("[playbase.ingest::read_h5_counts] AnnData read failed: ", conditionMessage(w))
          NULL
        }
      )
    }
  }

  ## For AnnData path df is already list(counts, samples); unwrap for fallbacks
  if (is.list(df) && all(c("counts", "samples") %in% names(df))) {
    counts_mat <- df$counts
    samples_df <- df$samples
  } else {
    counts_mat <- df
    samples_df <- NULL
  }

  if (is.null(counts_mat)) {
    counts_mat <- tryCatch(
      {
        Seurat::Read10X_h5(h5.file)
      },
      error = function(w) {
        NULL
      }
    )
  }

  if (is.null(counts_mat)) {
    counts_mat <- tryCatch(
      {
        h5.readMatrix(h5.file)
      },
      error = function(w) {
        NULL
      }
    )
  }

  if (!is.null(counts_mat) & (all(class(counts_mat) %in% c("matrix", "array")) || is(counts_mat, "sparseMatrix"))) {
    message("[playbase.ingest::read_h5_counts] Reading operation successfully completed. \n")
  } else {
    message("[playbase.ingest::read_h5_counts] df is null!")
  }

  return(list(counts = counts_mat, samples = samples_df))
}

#' @title Read (sub-)matrix from HDF5 file
#'
#' @export
h5.readMatrix <- function(h5.file, rows = NULL, cols = NULL,
                          matrixid = "data/matrix", rowid = "data/rownames",
                          colid = "data/colnames") {
  cn <- rhdf5::h5read(h5.file, colid)
  rn <- rhdf5::h5read(h5.file, rowid)
  rowidx <- 1:length(rn)
  colidx <- 1:length(cn)
  if (!is.null(rows)) {
    if (all(is.integer(rows))) rows <- rn[rows]
    rowidx <- match(intersect(rows, rn), rn)
  }
  if (!is.null(cols)) {
    if (all(is.integer(cols))) cols <- cn[cols]
    colidx <- match(intersect(cols, cn), cn)
  }
  nr <- length(rowidx)
  nc <- length(colidx)
  message("[sigdb.getConnectivityMatrix] reading large H5 file: ", nr, "x", nc, "")
  X <- rhdf5::h5read(h5.file, matrixid, index = list(rowidx, colidx))
  rownames(X) <- rn[rowidx]
  colnames(X) <- cn[colidx]
  return(X)
}
