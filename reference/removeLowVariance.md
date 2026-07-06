# Filter low variance features

Filter low variance features

## Usage

``` r
removeLowVariance(sce, threshold = 0, robust = FALSE)
```

## Arguments

- sce:

  [`SingleCellExperiment`](https://rdrr.io/pkg/SingleCellExperiment/man/SingleCellExperiment.html)
  object

- threshold:

  Keep features that have larger variance than this threshold

- robust:

  If true use median absolute deviation, otherwise use variance

## Value

[`SingleCellExperiment`](https://rdrr.io/pkg/SingleCellExperiment/man/SingleCellExperiment.html)
object

## Examples

``` r
set.seed(23)
cell_file <- generate_data()
sce <- loadData(cell_file)
sce <- removeLowVariance(sce)
```
