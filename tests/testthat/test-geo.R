test_that("is.GEO.id.valid accepts GSE ids only", {
  expect_true(playbase.ingest::is.GEO.id.valid("GSE10846"))
  expect_false(playbase.ingest::is.GEO.id.valid("foo"))
})
