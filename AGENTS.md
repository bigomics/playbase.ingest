# playbase.ingest: guide for agents and contributors

`playbase.ingest` is the data-ingestion code of OmicsPlayground, split out of
`playbase`. It does two things:

- **Read** user files and public repositories into plain R objects (matrices,
  data frames).
- **Validate** uploads before a pgx object is built.

Anything that computes on the data belongs in another package: normalization,
annotation, clustering, statistics.

## Hard rules

1. **This is a leaf package.** Never depend on `playbase`,
   `playbase.preprocess` or `playbase.epigenetics`. That means no Imports and
   no `playbase::` calls. `playbase` depends on this package, so the reverse
   is a cycle.
2. **Need a playbase helper? Copy it.** Put the copy in `R/utils.R`, marked
   `@noRd` and not exported. The five helpers already there are copies of
   this kind.
3. **Exported names are API.** `playbase` and downstream apps call them as
   `playbase.ingest::fn()`. Don't rename exports or change what they return
   without checking those callers.
4. **playbase never redefines these names.** It does not re-export this
   package, and its `tests/testthat/test-reexports.R` fails if playbase
   defines a function with the same name as one exported here. When you move
   a function from playbase to this package, delete it from playbase in the
   same change.
5. **Every package an exported function calls goes in Imports.** Call it as
   `pkg::fn()`. playbase installs this package with its Imports only, so a
   package that is only in Suggests is missing at runtime and the feature
   breaks. Suggests is for test-only packages (`testthat`, `arrow`), which
   tests guard with `skip_if_not_installed()`.
6. **Epigenomics (IDAT) ingestion is left out on purpose for now.** It lives
   in `playbase.epigenetics`.

## Where code goes

File names follow the function family: `read-*`, `check*`, `validate`,
`geo-*`. Every source file has one test file, named
`tests/testthat/test-<file>.R`.

| File | Holds | Put new code here when… |
|---|---|---|
| `R/read.R` | `read_counts`, `read_samples`, `read_contrasts`, `read_annot`, `read_files`, `first_feature` | it reads one of the standard upload tables |
| `R/read-table.R` | `read.as_matrix`, `detect_delim`, `detect_decimal`, `fread.csv` | it's low-level delimited-file parsing (separators, decimals, headers) |
| `R/read-olink.R`, `R/read-spectronaut.R`, `R/read-h5.R`, `R/read-singlecell.R` | one platform or file format each | it reads a vendor or platform export. **New platform → new `R/read-<platform>.R`** |
| `R/read-gmt.R` | `read.gmt` | it reads gene-set files |
| `R/check.R` | `pgx.checkINPUT` | it's a check on a single table |
| `R/check-cross.R` | `pgx.crosscheckINPUT`, `contrasts_conversion_check` | it's a check across samples, counts and contrasts |
| `R/validate.R` | `validate_*`, `getError` | it's a user-facing validation wrapper |
| `R/geo.R` | `pgx.getGEOseries`, `is.GEO.id.valid` | it's the GEO entry point |
| `R/geo-counts.R` | `pgx.getGEOcounts` and its ARCHS4 / recount / ArrayExpress sources | it's another counts source in the `pgx.getGEOcounts()` fallback chain |
| `R/geo-geoquery.R` | the GEOquery path, `pgx.getSymbolFromFeatureData` | it's GEOquery-specific |
| `R/geo-metadata.R` | `pgx.getGEOmetadata*`, `pgx.getGEOexperimentInfo` | it fetches sample or series metadata |
| `R/geo-pheno.R` | `title2pheno`, `eset.*`, `trimsame*` | it parses phenotypes out of GEO titles |
| `R/data.R` | `example_file`, docs for data objects | it documents a data object |
| `R/utils.R` | private helper copies only | never put exported code here |
| `R/playbase.ingest-package.R` | package doc, `@import` / `@importFrom` | you add a package-wide import |

- **A new public repository** (not GEO or ArrayExpress) gets its own file,
  `R/<repo>.R` or `R/<repo>-*.R`.
- **Aim for files under ~300 lines.** When a file grows past that, split it
  by topic, not into arbitrary halves.

## Adding a reader

1. **Function:** `read_<format>(file, ...)`. Return a matrix (features ×
   samples), or `list(counts = , samples = )` when the file also carries
   sample metadata. Prefix messages with `[playbase.ingest::read_<format>]`.
2. **Roxygen:** a title, `@param` for every argument, `@return` and
   `@export`. `R CMD check` warns about any mismatch.
3. **Dependencies:** add any new package to Imports in `DESCRIPTION`.
4. **Fixture:**
   - small, shareable files go in `inst/extdata/` and are reached with
     `example_file()`;
   - otherwise use `tests/data/`, or write a temp file inside the test.
5. **Test:** add `tests/testthat/test-read-<format>.R`.
6. **Run** `devtools::document()`, `devtools::test()` and `R CMD check` (see
   below).

## Upload error codes (`PGX_CHECKS`)

- The error codes returned by `pgx.checkINPUT()` / `pgx.crosscheckINPUT()`
  (e.g. `e25`) are described in `data/PGX_CHECKS.rda`. `getError()` looks
  messages up there, and the app shows them to users.
- That file is built from `data-raw/PGX_CHECKS.csv` by
  `data-raw/pgx_checks_builder.R`.
- To change a code: edit the CSV, run the builder from the package root,
  then commit both the CSV and the `.rda`. Never edit the `.rda` by hand.

## Workflow

- **After roxygen changes,** run `devtools::document()`. Never edit
  `NAMESPACE` or `man/` by hand.
- **Before a PR,** `devtools::test()` must pass, and
  `rcmdcheck::rcmdcheck(args = "--no-manual")` must give 0 errors and no new
  warnings. Known warnings:
  - the non-standard license string;
  - `require(Seurat)` in `read_cellranger_output()`;
  - `GEOquery::getRNASeqData` missing, with older GEOquery versions.
- **Moving code between files:**
  - first commit a pure `git mv` of the largest piece, then do the split in
    a second commit;
  - merge with a merge commit or rebase, **not squash**, so
    `git log --follow` keeps the history.
- **Don't fix unrelated bugs while moving code.** Open a separate PR.
- **playbase installs this package from GitHub `HEAD`** (its Remotes and
  unittest image). So `main` reaches playbase CI right away: keep it green
  and keep exports stable.
