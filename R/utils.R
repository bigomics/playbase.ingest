## ponytail: private copies of playbase helpers that playbase still uses
## elsewhere. Kept unexported so the two packages share no API; see
## playbase.preprocess for the same leaf-owns-its-helpers convention.

#' @noRd
make_unique <- function(s, sep = ".") {
  s[is.na(s)] <- "NA"
  make.unique(s, sep = sep)
}

#' Check if values are logarithm
#'
#' @noRd
is_logged <- function(x, verbose = 0) {
  ## force as matrix
  if (any(class(x) == "data.frame")) x <- as.matrix(x)

  ## if all values are 'small' it may be log
  xmax <- max(x, na.rm = TRUE)
  all.lt60 <- xmax < 60
  has.bigx <- xmax > 1000

  ## if we have negative  values it may be log. The -1 because of
  ## possible RNAseq prior = 1 in log2(x+1)
  xmin <- min(x, na.rm = TRUE)
  minx.neg <- xmin < -1
  all.pos <- xmin >= 0

  ## check for ratio/fraction data: rows of ratio data matrices (like
  ## TMT) sum to some constant value.
  is.ratio <- FALSE
  if (NCOL(x) > 1 && all.pos) {
    rowx.mean <- Matrix::rowMeans(x, na.rm = TRUE)
    rowx.mean <- rowx.mean[!is.na(rowx.mean)]
    is.ratio <- (sd(rowx.mean) / mean(rowx.mean)) < 0.01
  }

  ## check raw count data: if all values are integer.
  ## For sparse matrices avoid mean(sparseMatrix) which dispatches to mean.default,
  ## returns NA, and then (x - round(x)) / NA densifies the full matrix by forcing
  ## every structural zero to an explicit NA. Instead operate on x@x directly.
  if (inherits(x, "sparseMatrix")) {
    xv <- x@x
    xv_mean <- if (length(xv) > 0) mean(abs(xv)) else 1
    all.integer <- length(xv) == 0 ||
      max(abs(xv - round(xv)), na.rm = TRUE) / xv_mean < 1e-8
  } else {
    fraction.diff <- (x - round(x)) / mean(x, na.rm = TRUE)
    all.integer <- max(fraction.diff, na.rm = TRUE) < 1e-8
  }
  is.counts <- all.pos && all.integer

  ## check for low-count single-cell data: they may be all <60
  ## so they look like log.
  ## For sparse matrices use nnzero(): avoids is.na(x)|x==0 which creates a
  ## dense-like logical matrix and then mean() densifies it again.
  if (inherits(x, "sparseMatrix")) {
    zero.inflated <- Matrix::nnzero(x) / (as.numeric(nrow(x)) * ncol(x)) < 0.5
  } else {
    zero.inflated <- mean((is.na(x) | x == 0)) > 0.5
  }
  is.singlecell <- NCOL(x) > 1000 && all.pos && all.lt60 && zero.inflated

  if (verbose > 0) {
    message("[is_logged] all.lt60 = ", all.lt60)
    message("[is_logged] minx.neg = ", minx.neg)
    message("[is_logged] has.bigx = ", has.bigx)
    message("[is_logged] all.pos = ", all.pos)
    message("[is_logged] is.ratio = ", is.ratio)
    message("[is_logged] all.integer = ", all.integer)
    message("[is_logged] is.counts = ", is.counts)
    message("[is_logged] zero.inflated = ", zero.inflated)
    message("[is_logged] is.singlecell = ", is.singlecell)
  }

  possible.log <- (all.lt60 || minx.neg)
  possible.linear <- (is.ratio || has.bigx || is.counts || is.singlecell)
  is.log <- possible.log && !possible.linear
  is.log
}

#' Convert contrast matrix to group labels
#'
#' @title Convert contrast matrix to group labels
#'
#' @description Converts a contrast matrix to a data frame of group labels.
#'
#' @param contr.matrix The contrast matrix.
#' @param as.factor Whether to return as factor. Default is FALSE.
#'
#' @details This function takes a contrast matrix and returns a data frame
#' with a column of group labels for each contrast. The labels are generated
#' from the contrast names.
#'
#' @return Data frame of group labels.
#'
#' @examples
#' \dontrun{
#' contrast <- playbase::CONTRASTS
#' z <- playbase::contrastAsLabels(contrast)
#' z
#' }
#' @noRd
contrastAsLabels <- function(contr.matrix, as.factor = FALSE) {
  contrastAsLabels.col <- function(contr, contr.name) {
    grp1 <- gsub(".*[:]|_vs_.*", "", contr.name)
    grp0 <- gsub(".*_vs_|@.*", "", contr.name)
    x <- rep(NA, length(contr))
    x[which(contr < 0)] <- grp0
    x[which(contr > 0)] <- grp1
    if (as.factor) x <- factor(x, levels = c(grp0, grp1))
    x
  }

  ## any data.frame to matrix
  contr.matrix <- as.matrix(contr.matrix)

  # convert values <0 to -1 and values >0 to 1 (for old design matrix with float values)
  # sign function will fail when non-numeric values are present, use try catch to preserve input
  contr.matrix <- tryCatch(sign(contr.matrix), error = function(e) contr.matrix)

  num.values <- c(-1, 0, 1, NA, "NA", "na", "", " ")
  is.num <- all(apply(contr.matrix, 2, function(x) all(x %in% num.values)))
  is.num <- is.num && all(c(-1, 1) %in% contr.matrix) ## must have at least -1 and 1!!
  if (!is.num) {
    ## message("[contrastAsLabels] already as label!")
    contr.matrix[which(contr.matrix == "")] <- NA
    return(contr.matrix)
  }
  K <- contr.matrix[, 0]
  i <- 1
  for (i in 1:ncol(contr.matrix)) {
    contr <- as.numeric(contr.matrix[, i])
    contr.name <- colnames(contr.matrix)[i]
    k1 <- contrastAsLabels.col(contr, contr.name)
    K <- cbind(K, k1)
  }
  K[which(K == "")] <- NA
  colnames(K) <- colnames(contr.matrix)
  rownames(K) <- rownames(contr.matrix)
  return(K)
}

#' Converts old-style contrast matrix to sample-wise labeled contrast matrix
#'
#' @title Convert old-style contrast matrix to sample-wise labeled matrix
#'
#' @param contrasts Matrix of contrasts
#' @param samples   Sample information dataframe
#'
#' @return Contrast matrix
#'
#' @description Convert old-style contrast to new-style contrast of sample-wise labels.
#'
#' @details This function takes a matrix of sample labels as input, where the column
#' names indicate the sample groups being compared (e.g. "Group1_vs_Group2"). It parses
#' the column names to extract the two groups, then constructs a labels contrast matrix by
#' assigning the proper condition name to samples belonging to each group.
#'
#' The resulting contrast matrix has rows corresponding to samples, and columns
#' corresponding to the label matrix column names. This encodes the contrasts between
#' each pair of groups.
#' @noRd
contrasts.convertToLabelMatrix <- function(contrasts, samples) {
  if (NCOL(contrasts) == 1 && is.null(dim(contrasts))) {
    stop("contrasts must be a matrix with column names")
  }
  if (is.null(colnames(contrasts))) {
    stop("contrasts must have column names")
  }
  contrasts <- type.convert(contrasts, as.is = TRUE)
  is.numeric.matrix <- inherits(contrasts, "matrix") &&
    all(apply(contrasts, 2, class) %in% c("integer", "numeric"))
  is.numeric.df <- inherits(contrasts, "data.frame") &&
    all(apply(contrasts, 2, class) %in% c("integer", "numeric"))
  is.numeric.contrast <- is.numeric.matrix | is.numeric.df
  ## Preserve numeric group labels (e.g. "1","2","3") as character
  if (is.numeric.contrast) {
    vals <- unique(as.vector(as.matrix(contrasts)))
    vals <- vals[!is.na(vals)]
    if (!all(vals %in% c(-1, 0, 1))) {
      rn <- rownames(contrasts)
      contrasts <- apply(contrasts, 2, as.character)
      contrasts[contrasts == "NA"] <- NA
      rownames(contrasts) <- rn
      is.numeric.contrast <- FALSE
    }
  }
  if (is.numeric.contrast) {
    has.negpos <- any(contrasts < 0, na.rm = TRUE) && any(contrasts > 0, na.rm = TRUE)
    if (has.negpos) {
      contrasts[contrasts %in% c(NA, "NA", "na", "")] <- 0
      contrasts <- sign(contrasts)
    }
  }

  ## first match of group (or condition) in colum names, regard as group column
  group.col <- head(grep("group|condition", tolower(colnames(samples))), 1)
  group.col
  if (length(group.col) == 0) {
    ## try to detect automatically
    group.col <- head(which(apply(samples, 2, function(x) all(rownames(contrasts) %in% x))), 1)
    group.col
  }
  has.group.col <- length(group.col) > 0
  is.group.contrast <- has.group.col && all(samples[, group.col] %in% rownames(contrasts))
  is.sample.contrast <- nrow(contrasts) > 0 && all(rownames(contrasts) %in% rownames(samples))

  if (!is.sample.contrast && !has.group.col) {
    message("[contrasts.convertToLabelMatrix] ERROR: Invalid group-wise contrast. could not find 'group' column.")
    return(NULL)
  }
  if (!is.group.contrast && !is.sample.contrast) {
    message("[contrasts.convertToLabelMatrix] ERROR: Invalid contrast.")
    return(NULL)
  }

  new.contrasts <- contrasts
  ## old1: group-wise -1/0/1 matrix
  if (is.group.contrast && is.numeric.contrast) {
    message("[contrasts.convertToLabelMatrix] WARNING: converting old1 style contrast to new format")
    new.contrasts <- samples[, 0]
    if (NCOL(contrasts) > 0) {
      contrasts2 <- apply(contrasts, 2, as.numeric)
      contrasts2 <- contrastAsLabels(contrasts2)
      rownames(contrasts2) <- rownames(contrasts)
      grp <- as.character(samples[, group.col])
      new.contrasts <- contrasts2[grp, , drop = FALSE]
      rownames(new.contrasts) <- rownames(samples)
    }
  }
  ## old2: sample-wise -1/0/1 matrix
  if (is.sample.contrast && is.numeric.contrast) {
    message("[contrasts.convertToLabelMatrix] WARNING: converting old2 style contrast to new format")
    new.contrasts <- samples[, 0]

    if (NCOL(contrasts) > 0) {
      contrasts2 <- apply(contrasts, 2, as.numeric)
      rownames(contrasts2) <- rownames(contrasts)
      contrasts2 <- contrastAsLabels(contrasts2)
      new.contrasts <- matrix(NA, nrow(samples), ncol(contrasts))
      new.contrasts <- contrasts2[match(rownames(samples), rownames(contrasts2)), , drop = FALSE]
      rownames(new.contrasts) <- rownames(samples)
    }
  }
  ## old3: group-wise label matrix
  if (is.group.contrast && !is.numeric.contrast) {
    message("[contrasts.convertToLabelMatrix] WARNING: converting group-wise label contrast to new format")
    new.contrasts <- samples[, 0]
    if (NCOL(contrasts) > 0) {
      grp <- as.character(samples[, group.col])
      new.contrasts <- contrasts[grp, , drop = FALSE]
      rownames(new.contrasts) <- rownames(samples)
    }
  }

  ## always clean up
  new.contrasts <- as.matrix(new.contrasts)
  new.contrasts <- apply(new.contrasts, 2, as.character)
  new.contrasts[trimws(new.contrasts) %in% c("", "NA", "na", " ")] <- NA
  rownames(new.contrasts) <- rownames(samples)
  new.contrasts
}

#' @noRd
iconv2utf8 <- function(s) {
  iconv(s, to = "UTF-8//TRANSLIT", sub = "")
}
