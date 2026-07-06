# Remove zero-inflated features

Remove zero-inflated features

## Usage

``` r
removeZeroInflation(sce, proportion = 0.2)
```

## Arguments

- sce:

  [`SingleCellExperiment`](https://rdrr.io/pkg/SingleCellExperiment/man/SingleCellExperiment.html)
  object

- proportion:

  Remove features exceeding this zero-inflation proportion

## Value

[`SingleCellExperiment`](https://rdrr.io/pkg/SingleCellExperiment/man/SingleCellExperiment.html)
object

## Examples

``` r
set.seed(23)
cell_file <- generate_data()
sce <- loadData(cell_file)
sce <- removeZeroInflation(sce)
```
