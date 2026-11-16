# Plot correlations between features and principal components

Computes the Pearson correlation of every feature with the first
principal components in `reducedDim(sce, "PCA")` and shows the features
most correlated with one of them.

## Usage

``` r
plotPCACor(
  sce,
  filter_by = 1,
  top = 20,
  pcs = seq(5),
  assay_type = "tfmfeatures"
)
```

## Arguments

- sce:

  [`SingleCellExperiment`](https://rdrr.io/pkg/SingleCellExperiment/man/SingleCellExperiment.html)
  object

- filter_by:

  numeric(1), principal component used to select features

- top:

  numeric(1), number of features to show

- pcs:

  integer vector, principal components to show

- assay_type:

  character(1), name of the assay with the features

## Value

[`ggplot2`](https://ggplot2.tidyverse.org/reference/ggplot2-package.html)
object

## Examples

``` r
set.seed(23)
cell_file <- generateData()
sce <- loadData(cell_file)
sce <- transformScale(sce, robust = TRUE)
sce <- scater::runPCA(sce, exprs_values = "tfmfeatures", ncomponents = 10)
plotPCACor(sce, filter_by = 1)

```
