# playbase.ingest

Data ingestion for OmicsPlayground, split out of
[playbase](https://github.com/bigomics/playbase): delimited table readers
(`read.as_matrix`, `read_counts`, `read_samples`, `read_contrasts`,
`read_annot`, `fread.csv`), platform readers (Olink NPX, Spectronaut,
10X Cell Ranger, h5/h5ad), `read.gmt`, and the upload checks
(`pgx.checkINPUT`, `pgx.crosscheckINPUT`, `validate_*`, `PGX_CHECKS`).

A leaf package: it depends on nothing from playbase. Heavy readers
(OlinkAnalyze, rhdf5, Seurat) are Suggests and checked at call time.

```r
remotes::install_github("bigomics/playbase.ingest")
counts <- playbase.ingest::read_counts(playbase.ingest::example_file("counts.csv"))
```
