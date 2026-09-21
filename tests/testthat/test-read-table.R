test_that("read.as_matrix works correctly", {
  # Test that read.as_matrix correctly reads a simple matrix
  # write.table(matrix(1:4, nrow = 2), file = ".//tests/data/test1.csv", sep = ";", row.names = c("a", "b"), col.names = c("c","d") )
  temp_file1 <- testthat::test_path("..", "data", "test1.csv")
  expect_equal(as.numeric(read.as_matrix(temp_file1)), c(1, 2, 3, 4))

  # Test that read.as_matrix correctly reads a matrix with extreme values
  # write.table(matrix(c(1e100, 1e-100, -1e100, -1e-100), nrow = 2), file = ".//tests/data/test2.csv", sep = ";", row.names = c("a", "b"), col.names = c("c","d") )
  temp_file2 <- testthat::test_path("..", "data", "test2.csv")
  expect_equal(as.numeric(read.as_matrix(temp_file2)), c(1e100, 1e-100, -1e100, -1e-100))

  temp_file3 <- testthat::test_path("..", "data", "large_integers.csv")
  expect_equal(
    as.numeric(read.as_matrix(temp_file3)),
    c(395000000, 895050000, 84760000000, 4760700000, 2390000000, 1290000000, 4680000000, 4680000000)
  )
})

test_that("read.as_matrix reads file as matrix", {
  tmp <- tempfile()
  writeLines(c(",1,2,3", "gene1,1,2,3", "gene2,4,5,6"), tmp)
  expected <- matrix(c(1, 2, 3, 4, 5, 6),
    ncol = 3, byrow = TRUE,
    dimnames = list(c("gene1", "gene2"), NULL)
  )
  colnames(expected) <- c(1, 2, 3)
  result <- playbase.ingest::read.as_matrix(tmp, as.char = FALSE)
  expect_equal(class(result), c("matrix", "array"))
  expect_equal(result, expected)
  unlink(tmp)
})
