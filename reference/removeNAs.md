# Remove cells with missing features

Remove cells with missing features

## Usage

``` r
removeNAs(sce, assay_type = "features")
```

## Arguments

- sce:

  [`SingleCellExperiment`](https://rdrr.io/pkg/SingleCellExperiment/man/SingleCellExperiment.html)
  object

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
sce <- removeNAs(sce)
```
