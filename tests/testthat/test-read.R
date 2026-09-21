testthat::test_that("read_counts works", {
  counts <- playbase.ingest::read_counts(
    playbase.ingest::example_file("counts.csv")
  )
  testthat::expect_equal(nrow(counts), 7439)
})

testthat::test_that("read_samples works", {
  samples <- playbase.ingest::read_samples(playbase.ingest::example_file("samples.csv"))
  testthat::expect_equal(nrow(samples), 18)
})

testthat::test_that("read_contrasts works", {
  contrasts <- playbase.ingest::read_contrasts(playbase.ingest::example_file("contrasts.csv"))
  testthat::expect_equal(nrow(contrasts), 6)
})
