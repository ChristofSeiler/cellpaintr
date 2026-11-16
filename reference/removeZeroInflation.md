# Remove zero-inflated features

Remove zero-inflated features

## Usage

``` r
removeZeroInflation(sce, proportion = 0.2, assay_type = "features")
```

## Arguments

- sce:

  [`SingleCellExperiment`](https://rdrr.io/pkg/SingleCellExperiment/man/SingleCellExperiment.html)
  object

- proportion:

  numeric(1) between 0 and 1, remove features with a larger proportion
  of zeros

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
sce <- removeZeroInflation(sce)
```
