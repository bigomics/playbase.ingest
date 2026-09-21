#' Test for pgx.crosscheckINPUT
# Define the test data
COUNTS <- data.frame(matrix(rnorm(40), nrow = 10))

colnames(COUNTS) <- c("sample1", "sample2", "sample3", "sample4")

rownames(COUNTS) <- paste0("gene", 1:10)

SAMPLES <- data.frame(
  type = c("A", "A", "B", "B"),
  group = c("high", "low", "high", "low")
)

rownames(SAMPLES) <- colnames(COUNTS)

CONTRASTS <- data.frame(
  `type:A_vs_B` = c("A", "A", "B", "B"),
  `group:high_vs_low` = c("high", "low", "high", "low")
)

# TODO: make test work with samples, counts and contrasts above

# Test 2: Check if the function returns a list
test_that("crosscheckINPUT function returns a list", {
  result <- playbase.ingest::pgx.crosscheckINPUT(
    playbase.ingest::read_samples(playbase.ingest::example_file("samples.csv")),
    playbase.ingest::read_counts(playbase.ingest::example_file("counts.csv")),
    playbase.ingest::read_contrasts(playbase.ingest::example_file("contrasts.csv"))
  )
  expect_type(result, "list")
  expect_equal(result$checks, list())
  expect_true(result$PASS)
})

# # Test 3: Check if the function handles non-matching samples
# test_that("crosscheckINPUT function handles non-matching sample and count names", {
#   SAMPLES_mismatch <- SAMPLES
#   rownames(SAMPLES_mismatch)[1] <- "mismatch"
#   suppressWarnings(result <- playbase.ingest::pgx.crosscheckINPUT(SAMPLES_mismatch, COUNTS, CONTRASTS))
#   expect_true(result$PASS)
#   expect_equal(result$checks$e19, c("mismatch", "sample1"))
#   expect_equal(result$checks$e17, c("sample2", "sample3", "sample4", "1", "2", "3", "4"))
# })

# Test 6: Check if the function handles non-matching order of sample and count names
# test_that("crosscheckINPUT function handles non-matching order of sample and count names", {
#   COUNTS_mismatch_order <- COUNTS[, 2:1]
#   result <- playbase.ingest::pgx.crosscheckINPUT(SAMPLES, COUNTS_mismatch_order, CONTRASTS)
#   expect_true(result$PASS)
#   expect_equal(result$checks$e19, c("sample3", "sample4"))
#   expect_equal(result$checks$e17, c("sample1", "sample2", "1", "2", "3", "4"))
# })

#' Test for contrasts_conversion_check
