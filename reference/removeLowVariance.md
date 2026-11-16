# Filter low variance features

Filter low variance features

## Usage

``` r
removeLowVariance(sce, threshold = 0, robust = FALSE, assay_type = "features")
```

## Arguments

- sce:

  [`SingleCellExperiment`](https://rdrr.io/pkg/SingleCellExperiment/man/SingleCellExperiment.html)
  object

- threshold:

  numeric(1), keep features whose spread exceeds it

- robust:

  logical(1), if `TRUE` the spread is the median absolute deviation,
  otherwise the variance

- assay_type:

  character(1), name of the assay with the features

## Value

[`SingleCellExperiment`](https://rdrr.io/pkg/SingleCellExperiment/man/SingleCellExperiment.html)
object

## Examples

``` r
set.seed(23)
cell_file <- generateData()
sce <- loadData(cell_file)
sce <- removeLowVariance(sce)
```
