# Remove cells if not enough or too many in one image

Remove cells if not enough or too many in one image

## Usage

``` r
removeOutliers(sce, min, max)
```

## Arguments

- sce:

  [`SingleCellExperiment`](https://rdrr.io/pkg/SingleCellExperiment/man/SingleCellExperiment.html)
  object

- min:

  Remove cells below that number

- max:

  Remove cells above that number

## Value

[`SingleCellExperiment`](https://rdrr.io/pkg/SingleCellExperiment/man/SingleCellExperiment.html)
object

## Examples

``` r
set.seed(23)
cell_file <- generate_data()
sce <- loadData(cell_file)
sce <- removeOutliers(sce, min = 0, max = 300)
```
