#' Test for pgx.checkINPUT
# Load the necessary libraries

# Define the common test data
df <- data.frame(matrix(rnorm(20), nrow = 10))

colnames(df) <- c("sample1", "sample2")

rownames(df) <- paste0("gene", 1:10)

# Test 1:
test_that("checkINPUT function returns a list", {
  expect_type(playbase.ingest::pgx.checkINPUT(df, "SAMPLES"), "list")
})

# Test 2:
test_that("checkINPUT function handles non-existent datatype", {
  expect_error(playbase.ingest::pgx.checkINPUT(df, "NON_EXISTENT"))
})

# Test 3:
test_that("checkINPUT function handles duplicated column names", {
  df_dup <- cbind(df, df)
  check_result <- playbase.ingest::pgx.checkINPUT(df_dup, "COUNTS")
  expect_false(check_result$PASS)
  expect_equal(check_result$checks$e6, c("sample1", "sample2"))
})

# Test 4:
test_that("checkINPUT function handles duplicated row names", {
  df_dup <- as.matrix(df)
  rownames(df_dup)[1] <- rownames(df_dup)[2]

  # No pass for samples
  samples_check_result <- playbase.ingest::pgx.checkINPUT(df_dup, "SAMPLES")
  expect_false(samples_check_result$PASS)
  expect_equal(samples_check_result$checks$e1, "gene2")

  # Pass for counts
  counts_check_result <- playbase.ingest::pgx.checkINPUT(df_dup, "COUNTS")
  expect_true(counts_check_result$PASS)
  expect_equal(counts_check_result$checks$e7, "gene2")
})

# Test 5:
test_that("checkINPUT function handles zero count rows", {
  df_zero <- df
  df_zero[c(1, 5), ] <- 0
  check_result <- playbase.ingest::pgx.checkINPUT(df_zero, "COUNTS")
  expect_true(check_result$PASS)
  expect_equal(check_result$checks$e9, c("gene1", "gene5"))
})

# Test 6:
test_that("checkINPUT function handles zero count columns", {
  df_zero <- df
  df_zero[, 1] <- 0
  counts_zero_col_check <- playbase.ingest::pgx.checkINPUT(df_zero, "COUNTS")
  expect_true(counts_zero_col_check$PASS)
  expect_equal(counts_zero_col_check$checks$e10, "sample1")
})

# Test 7:
test_that("checkINPUT function handles valid contrast names", {
  df_contrasts <- df
  colnames(df_contrasts) <- c("sample1_sample2", "sample2_sample3")
  expect_equal(pgx.checkINPUT(df_contrasts, "CONTRASTS")$PASS, FALSE)
})

# Test 8:
test_that("checkINPUT function handles invalid contrast names", {
  df_contrasts <- df
  colnames(df_contrasts) <- c("sample1_sample2", "sample2_sample3")
  expect_equal(playbase.ingest::pgx.checkINPUT(df_contrasts, "CONTRASTS")$PASS, FALSE)
})
