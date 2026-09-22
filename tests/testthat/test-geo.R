test_that("is.GEO.id.valid accepts GSE ids only", {
  expect_true(playbase.ingest::is.GEO.id.valid("GSE10846"))
  expect_false(playbase.ingest::is.GEO.id.valid("foo"))
})

test_that("pgx.getGEOseries skips ARCHS4 unless a file is given", {
  expect_null(formals(playbase.ingest::pgx.getGEOseries)$archs.h5)
})
