test_that("read_multiomics_counts prefixes layers and keeps shared samples", {
  rna <- matrix(1:6, 2, 3, dimnames = list(c("A", "B"), c("s1", "s2", "s3")))
  prot <- matrix(7:10, 2, 2, dimnames = list(c("A", "C"), c("s2", "s1")))

  expect_message(
    mo <- read_multiomics_counts(list(rna, prot), c("RNA-seq", "proteomics")),
    "s3"
  )
  expect_equal(rownames(mo$counts), c("gx:A", "gx:B", "px:A", "px:C"))
  expect_equal(colnames(mo$counts), c("s1", "s2"))
  expect_equal(unname(mo$counts["px:A", ]), c(9, 7))
  expect_equal(mo$dropped_samples, "s3")
})

test_that("read_multiomics_counts rejects bad input", {
  m <- matrix(1:4, 2, 2, dimnames = list(c("A", "B"), c("s1", "s2")))
  n <- matrix(1:4, 2, 2, dimnames = list(c("A", "B"), c("x1", "x2")))
  expect_error(read_multiomics_counts(list(m), "RNA-seq"), "two layers")
  expect_error(read_multiomics_counts(list(m, m), c("RNA-seq", "microarray")), "only once")
  expect_error(read_multiomics_counts(list(m, m), c("RNA-seq", "lipidomics")), "unsupported")
  expect_error(read_multiomics_counts(list(m, n), c("RNA-seq", "proteomics")), "share no sample")
})
