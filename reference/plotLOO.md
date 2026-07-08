# Plot predicted leave-one-out probabilities

Plot predicted leave-one-out probabilities

## Usage

``` r
plotLOO(sce, target, group, assay_type = "tfmfeatures")
```

## Arguments

- sce:

  [`SingleCellExperiment`](https://rdrr.io/pkg/SingleCellExperiment/man/SingleCellExperiment.html)
  object

- target:

  Name of target variable for prediction

- group:

  Grouping variable for cross-validation, e.g., patient

- assay_type:

  A string specifying the assay

## Value

[`ggplot2`](https://ggplot2.tidyverse.org/reference/ggplot2-package.html)
object

## Examples

``` r
set.seed(23)
cell_file <- generate_data()
sce <- loadData(cell_file)
sce <- transformLogScale(sce)

sce$Drug <- as.factor(sce$Drug)
sce$Drug <- relevel(sce$Drug, ref = "D1")
types <- c("AreaShape", "Intensity", "Texture")

sce_single <- predictLOO(
    sce,
    target = "Drug", group = "Patient",
    interest_level = "D7", reference_level = "D1",
    types = types,
    n_threads = 1
)

plotLOO(sce_single, target = "Drug", group = "Patient")

```
