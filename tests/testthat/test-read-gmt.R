test_that("read.gmt reads gene sets", {
  tmp <- tempfile(fileext = ".gmt")
  writeLines(c("SET1\tsrc1\tA\tB\tC", "SET2\tsrc2\tD"), tmp)
  gmt <- playbase.ingest::read.gmt(tmp)
  expect_equal(gmt, list(SET1 = c("A", "B", "C"), SET2 = "D"))
  expect_equal(names(playbase.ingest::read.gmt(tmp, add.source = TRUE)), c("SET1 (src1)", "SET2 (src2)"))
  unlink(tmp)
})
