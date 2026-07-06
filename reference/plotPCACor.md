# Plot number of cells per image

Plot number of cells per image

## Usage

``` r
plotPCACor(sce, filter_by = 1, top = 20, pcs = seq(5))
```

## Arguments

- sce:

  [`SingleCellExperiment`](https://rdrr.io/pkg/SingleCellExperiment/man/SingleCellExperiment.html)
  object

- filter_by:

  PC to use for feature selection

- top:

  Number of top features to select

- pcs:

  Number of PCs to plot

## Value

[`ggplot2`](https://ggplot2.tidyverse.org/reference/ggplot2-package.html)
object

## Examples

``` r
set.seed(23)
cell_file <- generate_data()
sce <- loadData(cell_file)
sce <- transformLogScale(sce, robust = TRUE)
sce <- scater::runPCA(sce, exprs_values = "tfmfeatures", ncomponents = 10)
plotPCACor(sce, filter_by = 1)

```
