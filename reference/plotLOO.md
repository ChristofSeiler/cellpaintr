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

  character(1), name of the column of `colData(sce)` with the condition
  of each cell, such as the drug

- group:

  character(1), name of the column of `colData(sce)` with the replicate
  of each cell, such as the patient; the replicates are the folds of the
  leave-one-out cross-validation

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
sce <- transformScale(sce)

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
#>  ■■■■■■■■■■■■■■■■                  50% |  ETA:  1s

plotLOO(sce_single, target = "Drug", group = "Patient")

```
