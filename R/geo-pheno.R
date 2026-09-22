#' Extract phenotype data field from ExpressionSet
#' @param eset ExpressionSet object
#' @param field Character string specifying phenotype data field name
#' @return Vector of values for the specified field
#' @details Extracts a specific phenotype data field from an ExpressionSet object.
#' The \code{field} parameter specifies the phenotype data column name to extract.
#' The phenotype data pData slot is extracted from the ExpressionSet; the specified
#' field is returned as a vector.
#' @export
eset.getPhenoData <- function(eset, field) {
  Biobase::pData(Biobase::phenoData(eset))[, field]
}

#' Parse phenotype columns from sample titles
#' @param title character vector with the titles
#' @param split delimiter character to split on `c(",", ";", "\\|", "_", " ")`.
#' @param trim trim leading and trailing whitespace from each term
#' @param summarize summarize the terms by counting the number of occurrences of each term
#' @export
title2pheno <- function(title, split = NULL, trim = TRUE, summarize = TRUE) {
  if (is.null(split)) {
    split.chars <- c(",", ";", "\\|", "_", " ")
    i <- 1
    ss <- c()
    for (i in 1:length(title)) {
      ns <- sapply(split.chars, function(s) sum(gregexpr(s, title[i])[[1]] > 0))
      split0 <- names(ns)[which.max(ns)]
      ns1 <- ns[setdiff(names(ns), " ")]
      if (split0 == " " && any(ns1 > 0)) {
        split0 <- names(ns1)[which.max(ns1)]
      }
      ss[i] <- split0
    }
    split <- names(which.max(table(ss)))
  }

  ## Check if all titles have equal splitted parts (e.g. nicely formatted)
  nsplit <- sapply(title, function(tt) sum(gregexpr(split, tt)[[1]] > 0))
  nsplit.equal <- all(nsplit == nsplit[1])
  nsplit.equal
  if (!nsplit.equal) {
    cat("splitted title terms not equal lengths\n")
    return(NULL)
  }

  ## split
  tt <- as.character(sapply(as.character(title), function(s) trimws(s)))
  ff <- sapply(as.character(tt), strsplit, split = split)

  ## cleanup
  ff <- lapply(ff, trimws) ## trim whitespace
  ff <- lapply(ff, function(s) gsub("hours$|hour$|hrs$|hr$", "h", s)) ## hours
  ff <- lapply(ff, function(s) gsub("[ ][ ]*", " ", s)) ## double space

  ## make dataframe
  F1 <- do.call(rbind, ff)
  F1[is.na(F1)] <- NA
  rownames(F1) <- NULL

  ## Guess names
  getmax.term <- function(s) {
    tt <- table(unlist(strsplit(s, split = "[ _]")))
    vip.tt <- grep("hour|repl|hr|time|treat|infec|pati|sampl", names(tt))
    if (length(vip.tt)) tt[vip.tt] <- 1.1 * tt[vip.tt] ## boost known keywords
    names(which.max(tt))
  }
  maxterm <- apply(F1, 2, function(s) getmax.term(s))
  maxterm <- paste0("_", maxterm)
  colnames(F1) <- maxterm

  ## trims same words/characters on both ends
  if (trim) F1 <- apply(F1, 2, trimsame, summarize = summarize)
  F1
}

#' @describeIn eset.getPhenoData Phenotype data from the title of an
#' ExpressionSet object by splitting the title on a specified delimiter and guessing column names.
#' @param title Character vector of sample titles.
#' @param split Delimiter to split titles on; guessed when NULL.
#' @export
eset.parsePhenoFromTitle <- function(title, split = NULL) {
  if (!all(grepl(split, title))) {
    return(NULL)
  }

  tt <- as.character(sapply(as.character(title), function(s) trimws(s)))
  tt <- sapply(tt, function(s) gsub("[ ]*hours|[ ]*hour|[ ]*hrs|[ ]*hr", "h", s)) ## hours
  tt <- sapply(tt, function(s) gsub("([0-9]*)[ _]h", "\\1h", s))
  tt <- trimsame(tt, split = split, ends = TRUE)
  tt <- gsub(paste0(split, split, split), split, tt)
  tt <- gsub(paste0(split, split), split, tt)
  ff <- sapply(as.character(tt), strsplit, split = split)
  nf <- max(sapply(ff, length))
  ff <- lapply(ff, function(x) Matrix::head(c(x, rep(NA, nf)), nf))

  ## cleanup
  ff <- lapply(ff, trimws) ## trim whitespace
  ff <- lapply(ff, function(s) gsub("hours$|hour$|hrs$|hr$", "h", s)) ## hours
  ff <- lapply(ff, function(s) gsub("[ ][ ]*", " ", s)) ## double space

  F1 <- do.call(rbind, ff)
  F1[is.na(F1)] <- NA
  lapply(apply(F1, 2, table), sort, decreasing = TRUE)

  AA <- setdiff(unique(Biostrings::GENETIC_CODE), "*")

  i <- 1
  G <- list()
  for (i in 1:(ncol(F1) - 1)) {
    k <- min(ncol(F1), (i + 1))
    a2 <- factor(as.vector(F1[, i:k]))
    aa.dict <- levels(a2)
    names(aa.dict) <- AA[1:length(levels(a2))]
    levels(a2) <- AA[1:length(levels(a2))]
    F2 <- matrix(a2, nrow(F1))
    F2[is.na(F2)] <- "-"

    ff <- apply(F2, 1, paste, collapse = "")
    names(ff) <- paste0("tt", 1:nrow(F2))
    aln <- as.character(msa::msa(ff, type = "protein"))
    aln <- aln[names(ff)]
    F.aln <- do.call(rbind, sapply(aln, strsplit, split = ""))
    F.aln2 <- apply(F.aln, 2, function(x) aa.dict[x])

    if (ncol(F.aln2) > ncol(F2) && i < (ncol(F1) - 1)) {
      G[[i]] <- F.aln2[, 1:(ncol(F.aln) - 1)]
    } else if (i == (ncol(F1) - 1)) {
      G[[i]] <- F.aln2
    } else {
      G[[i]] <- F.aln2[, 1]
    }
  }
  G <- do.call(cbind, G)
  G <- G[, colMeans(is.na(G)) < 1, drop = FALSE]
  rownames(G) <- NULL
  colnames(G) <- paste0("V", 1:ncol(G))

  return(G)
}

#' @describeIn trimsame0 trimsame is a function that trims common prefixes and/or
#' suffixes from a character vector by applying trimsame0 forwards and/or backwards.
#' @param ends Logical. Trim common words at both ends (TRUE) or only the
#'   prefix (FALSE).
#' @export
trimsame <- function(s, split = " ", ends = TRUE, summarize = FALSE) {
  if (all(is.na(s)) || all(s == "")) {
    return(s)
  }
  if (ends) {
    return(trimsame.ends(s, split = split, summarize = summarize))
  }
  return(trimsame0(s, split = split, summarize = summarize))
}

#' @describeIn trimsame0 trimsame.ends is a function that trims common prefixes and
#' suffixes from a character vector by applying trimsame0 forwards and backwards.
#' @export
trimsame.ends <- function(s, split = " ", summarize = FALSE) {
  s1 <- trimsame0(s, split = split, summarize = summarize)
  s2 <- sapply(strsplit(s1, split = split), function(x) paste(rev(x), collapse = split))
  s3 <- trimsame0(s2, split = split, summarize = summarize, rev = TRUE)
  s4 <- sapply(strsplit(s3, split = split), function(x) paste(rev(x), collapse = split))
  s4
}

#' @title Trim Common Prefix from Strings
#'
#' @description This function trims the common prefix from a character vector of strings.
#'
#' @param s A character vector of strings to be trimmed.
#' @param split An optional character string specifying the delimiter used to split the strings into words.
#' The default value is a space (" ").
#' @param summarize An optional logical value indicating whether to summarize the common prefix using the first letter of each word.
#' The default value is `FALSE`.
#' @param rev An optional logical value indicating whether to reverse the order of the words in the summarized prefix.
#' The default value is `FALSE`.
#'
#' @details The function takes a character vector of strings `s` as input and searches for a common prefix among the strings.
#' The common prefix is defined as the longest sequence of characters that is shared by all strings and ends with the specified delimiter.
#' If a common prefix is found, it is trimmed from all strings.
#' If the `summarize` parameter is `TRUE`, the common prefix is summarized using the first letter of each word and appended to the beginning of each string.
#' If the `rev` parameter is `TRUE`, the order of the words in the summarized prefix is reversed.
#'
#' @return A character vector of the same length as `s`, containing the trimmed strings.
#'
#' @export
trimsame0 <- function(s, split = " ", summarize = FALSE, rev = FALSE) {
  if (all(is.na(s)) || all(s == "")) {
    return(s)
  }
  s <- strsplit(s, split)
  s.orig <- s

  ##
  i <- 1
  done <- FALSE
  while (i < 1000 && !done) {
    s1 <- sapply(s, "[", 1)
    slen <- sapply(s, length)
    ss <- setdiff(s1, NA)
    if (all(ss == ss[1])) {
      s <- lapply(s, "[", -1)
    } else {
      done <- TRUE
    }
    i <- i + 1
  }
  sapply(s, length)

  if (all(sapply(s, length) == 0) || all(s == "")) {
    sx <- sapply(s.orig, "[", 1)
    sx
    return(sx)
  }

  i <- 1
  done <- FALSE
  while (i < 1000 && !done) {
    slen <- sapply(s, length)
    ## tail() on an emptied element returns character(0), which would make
    ## sapply return a list and the comparison below error. Pad with NA so
    ## setdiff() drops it, matching how "[" handles the leading-token loop.
    s2 <- sapply(s, function(x) if (length(x)) utils::tail(x, 1) else NA_character_)
    ss <- setdiff(s2, NA)
    if (all(ss == ss[1])) {
      ## SIMPLIFY = FALSE: when every element has the same remaining length
      ## mapply would collapse the list into a matrix, and the final
      ## sapply(s, paste) would then iterate cells instead of strings,
      ## returning more elements than were passed in.
      s <- mapply(head, s, slen - 1, SIMPLIFY = FALSE)
    } else {
      done <- TRUE
    }
    i <- i + 1
  }

  if (all(sapply(s, length) == 0) || all(s == "")) {
    s <- sapply(s.orig, "[", 1)
  }
  s <- sapply(s, paste, collapse = split)
  s
}
