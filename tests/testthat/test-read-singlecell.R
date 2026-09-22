counts <- matrix(c(0, 3, 1, 0, 2, 5), 2,
  dimnames = list(c("g1", "g2"), c("c1", "c2", "c3"))
)

test_that("pgx.read_singlecell_counts reads csv with gene rownames", {
  f <- tempfile(fileext = ".csv")
  write.csv(counts, f)
  x <- playbase.ingest::pgx.read_singlecell_counts(f)
  expect_true(is.numeric(x))
  expect_equal(dimnames(x), dimnames(counts))
})

test_that("pgx.read_singlecell_counts reads 10X mtx with gene names", {
  d <- tempfile()
  dir.create(d)
  Matrix::writeMM(Matrix::Matrix(counts, sparse = TRUE), file.path(d, "matrix.mtx"))
  write.table(data.frame(c("ENSG1", "ENSG2"), rownames(counts)), file.path(d, "genes.tsv"),
    sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE
  )
  writeLines(colnames(counts), file.path(d, "barcodes.tsv"))
  x <- playbase.ingest::pgx.read_singlecell_counts(file.path(d, "matrix.mtx"))
  expect_equal(dimnames(x), dimnames(counts))
  expect_equal(as.matrix(x), counts)
})

test_that("seurat2pgx extracts counts and meta from a Seurat v5 object", {
  skip_if_not_installed("Seurat")
  m <- Matrix::Matrix(matrix(rpois(200, 5), 20,
    dimnames = list(paste0("G", 1:20), paste0("C", 1:10))
  ), sparse = TRUE)
  obj <- suppressWarnings(Seurat::CreateSeuratObject(m))
  pgx <- suppressWarnings(suppressMessages(playbase.ingest::seurat2pgx(obj)))
  expect_equal(dim(pgx$counts), c(20, 10))
  expect_equal(rownames(pgx$samples), colnames(m))
})

test_that("read_cellranger_output rejects unsupported file types", {
  expect_error(
    suppressMessages(playbase.ingest::read_cellranger_output("counts.csv")),
    "expected a .tar.gz, .gz or .zip file"
  )
})
