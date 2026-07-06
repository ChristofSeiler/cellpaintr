# Plot number of cells per image

Plot number of cells per image

## Usage

``` r
plotCellsPerImage(sce, bins = 100)
```

## Arguments

- sce:

  [`SingleCellExperiment`](https://rdrr.io/pkg/SingleCellExperiment/man/SingleCellExperiment.html)
  object

- bins:

  Number of histogram bins

## Value

[`ggplot2`](https://ggplot2.tidyverse.org/reference/ggplot2-package.html)
object

## Examples

``` r
set.seed(23)
cell_file <- generate_data()
sce <- loadData(cell_file)
plotCellsPerImage(sce)

```
