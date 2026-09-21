testthat::test_that("OlinkAnalyze::read_NPX works with parquet file", {
  testthat::skip_if_not_installed("OlinkAnalyze")
  testthat::skip_if_not_installed("arrow")
  parquet_file <- system.file("npx_data_ext.parquet", package = "playbase.ingest")
  testthat::expect_true(file.exists(parquet_file))
  npx <- OlinkAnalyze::read_NPX(parquet_file)
  testthat::expect_s3_class(npx, "data.frame")
  testthat::expect_gt(nrow(npx), 0)
})

testthat::test_that("OlinkAnalyze::read_NPX works with CSV file", {
  testthat::skip_if_not_installed("OlinkAnalyze")
  csv_file <- system.file("abundance_NPX_Data_3K.csv", package = "playbase.ingest")
  testthat::expect_true(file.exists(csv_file))
  npx <- OlinkAnalyze::read_NPX(csv_file)
  testthat::expect_s3_class(npx, "data.frame")
  testthat::expect_gt(nrow(npx), 0)
})

testthat::test_that("read_Olink_NPX returns counts and samples", {
  testthat::skip_if_not_installed("OlinkAnalyze")
  res <- playbase.ingest::read_Olink_NPX(system.file("abundance_NPX_Data_3K.csv", package = "playbase.ingest"))
  testthat::expect_true(is.matrix(res$counts))
  testthat::expect_equal(colnames(res$counts), rownames(res$samples))
})
