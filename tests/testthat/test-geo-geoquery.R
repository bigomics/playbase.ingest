test_that("pgx.getSymbolFromFeatureData maps an Ensembl column", {
  skip_if_not_installed("org.Hs.eg.db")
  skip_if_not_installed("AnnotationDbi")
  fdata <- data.frame(gene_assignment = c(
    "ENSG00000141510.18 // TP53", "NM_000546 // no ens", "ENSG00000146648 // EGFR"
  ))
  sym <- suppressMessages(playbase.ingest::pgx.getSymbolFromFeatureData(fdata))
  expect_equal(unname(sym), c("TP53", NA, "EGFR"))
})
