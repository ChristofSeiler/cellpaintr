# Aggregate predicted leave-one-out probabilities over meta variables over a list of SingleCellExperiment objects

Aggregate predicted leave-one-out probabilities over meta variables over
a list of SingleCellExperiment objects

## Usage

``` r
calculateStats(sce, target, group, assay_type = "tfmfeatures")
```

## Arguments

- sce:

  A
  [`SingleCellExperiment`](https://rdrr.io/pkg/SingleCellExperiment/man/SingleCellExperiment.html)
  object

- target:

  Name of target variable for prediction

- group:

  Grouping variable for cross-validation, e.g., patient

- assay_type:

  A string specifying the assay

## Value

`data.frame`

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

calculateStats(sce_single, target = "Drug", group = "Patient")
#>   Target   Feature       pvalue log2FoldChange
#> 1     D7       all 0.0008479434     0.36380804
#> 2     D7 AreaShape 0.8285342602    -0.04039513
#> 3     D7 Intensity 0.6854724293    -0.03239716
#> 4     D7   Texture 0.0002283086     0.54088688
```
