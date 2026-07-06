# Aggregate predicted leave-one-out probabilities over meta variables

Aggregate predicted leave-one-out probabilities over meta variables

## Usage

``` r
aggregateYhat(
  sce,
  assay_type = "tfmfeatures",
  meta_vars = c("Patient", "Treatment")
)
```

## Arguments

- sce:

  [`SingleCellExperiment`](https://rdrr.io/pkg/SingleCellExperiment/man/SingleCellExperiment.html)
  object

- assay_type:

  A string specifying the assay

- meta_vars:

  a vector of variables from \`colData\`

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

aggregateYhat(sce_single, meta_vars = c("Patient", "Drug"))
#> Error in aggregateYhat(sce_single, meta_vars = c("Patient", "Drug")): could not find function "aggregateYhat"
```
