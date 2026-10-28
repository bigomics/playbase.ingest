#' Test for title2pheno
test_that("title2pheno extracts phenotype terms", {
  titles <- c(
    "GSE1234-tissue_liver-disease-cancer",
    "GSE4321-tissue_lung-disease-cancer"
  )

  expected <- matrix(
    c(
      "GSE1234", "tissue_liver", "disease", "cancer",
      "GSE4321", "tissue_lung", "disease", "cancer"
    ),
    nrow = 2,
    dimnames = list(NULL, c("_GSE1234", "_tissue", "_disease", "_cancer"))
  )

  result <- playbase.ingest::title2pheno(titles, split = "-", trim = FALSE)

  # Check class, dim, and first column
  expect_equal(class(result), c("matrix", "array"))
  expect_equal(dim(result), c(length(titles), 4))
  expect_equal(result[, 1], c("GSE1234", "GSE4321"))
})

#' Test for trimsame
test_that("trimsame dispatches to trimsame0 or trimsame.ends", {
  # ends = FALSE delegates to trimsame0 (prefix trim only)
  expect_equal(
    playbase.ingest::trimsame(c("a b c", "a b d"), ends = FALSE),
    c("c", "d")
  )
  # ends = TRUE (default) trims both ends, non-space split
  expect_equal(
    playbase.ingest::trimsame(c("grp_a_b_c_end", "grp_z_end"), split = "_", ends = TRUE),
    c("a_b_c", "z")
  )
  # all-NA and all-empty short-circuit and return input unchanged
  expect_equal(playbase.ingest::trimsame(c(NA_character_, NA_character_)), c(NA_character_, NA_character_))
  expect_equal(playbase.ingest::trimsame(c("", "")), c("", ""))
})

#' Test for trimsame.ends
test_that("trimsame.ends trims common prefix and suffix", {
  expect_equal(
    unname(playbase.ingest::trimsame.ends(c("pre mid post", "pre other post"))),
    c("mid", "other")
  )
  # non-space split must round-trip correctly (collapse uses split, not " ")
  expect_equal(
    playbase.ingest::trimsame.ends(c("grp_a_b_c_end", "grp_z_end"), split = "_"),
    c("a_b_c", "z")
  )
})

#' Test for trimsame0
test_that("trimsame0 trims common leading/trailing tokens without erroring", {
  # length >= 2 input used to hard-error with
  # "'length = 2' in coercion to 'logical(1)'"
  expect_equal(playbase.ingest::trimsame0(c("a b c", "a b d")), c("c", "d"))

  # common prefix only
  expect_equal(playbase.ingest::trimsame0(c("a b c", "a b d e")), c("c", "d e"))

  # common suffix only
  expect_equal(playbase.ingest::trimsame0(c("a b c x", "z x")), c("a b c", "z"))

  # common prefix and suffix
  expect_equal(
    unname(playbase.ingest::trimsame0(c("pre mid post", "pre other post"))),
    c("mid", "other")
  )

  # no common leading/trailing token at all
  expect_equal(playbase.ingest::trimsame0(c("x y", "z w")), c("x y", "z w"))

  # non-space split
  expect_equal(playbase.ingest::trimsame0(c("a_b_c", "a_b_d"), split = "_"), c("c", "d"))

  # all-NA and all-empty short-circuit and return input unchanged
  expect_equal(playbase.ingest::trimsame0(c(NA_character_, NA_character_)), c(NA_character_, NA_character_))
  expect_equal(playbase.ingest::trimsame0(c("", "")), c("", ""))

  # trimming consumes every token (strings identical) -> falls back to
  # the first original token for each element
  expect_equal(playbase.ingest::trimsame0(c("a b", "a b")), c("a", "a"))
})

test_that("trimsame0 preserves input length and never errors on ragged input", {
  # equal remaining lengths after the suffix trim used to collapse the
  # intermediate list into a matrix, returning 4 elements for a 2-element input
  expect_equal(playbase.ingest::trimsame0(c("x y c", "z w c")), c("x y", "z w"))
  expect_length(playbase.ingest::trimsame0(c("x y c", "z w c")), 2L)

  # one element emptying before the others used to error with
  # "comparison of these types is not implemented"
  expect_equal(playbase.ingest::trimsame0(c("p a x", "p b x", "p x")), c("a", "b", ""))

  # length of the result always matches the length of the input
  for (input in list(
    c("a x", "b x", "c x"),
    c("p a x", "p b x", "p c x"),
    c("pre a post", "pre post"),
    c("one", "two", "three")
  )) {
    expect_length(playbase.ingest::trimsame0(input), length(input))
    expect_length(playbase.ingest::trimsame.ends(input), length(input))
  }
})

test_that("eset.getPhenoData works without Biobase attached", {
  pheno <- data.frame(group = c("a", "b"), row.names = c("s1", "s2"))
  eset <- Biobase::ExpressionSet(
    matrix(1:4, 2, dimnames = list(c("g1", "g2"), c("s1", "s2"))),
    phenoData = Biobase::AnnotatedDataFrame(pheno)
  )
  expect_equal(playbase.ingest::eset.getPhenoData(eset, "group"), c("a", "b"))
})
