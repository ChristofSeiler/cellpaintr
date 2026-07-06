# Predict target from features

Predict target from features

## Usage

``` r
predictLOO(
  sce,
  assay_type = "tfmfeatures",
  target = "Treatment",
  interest_level,
  reference_level,
  types = NULL,
  channels = NULL,
  group = "Patient",
  weights = NULL,
  n_threads = 1
)
```

## Arguments

- sce:

  [`SingleCellExperiment`](https://rdrr.io/pkg/SingleCellExperiment/man/SingleCellExperiment.html)
  object

- assay_type:

  A string specifying the assay

- target:

  Name of target variable for prediction

- interest_level:

  Factor interest level in \`target\` variable

- reference_level:

  Factor reference level in \`target\` variable

- types:

  Vector of strings of feature types

- channels:

  Vector of strings of feature channels

- group:

  Grouping variable for cross-validation, e.g., patient

- weights:

  Weights variable when features are aggregated

- n_threads:

  Number of parallel threads for fitting of models

## Value

[`tibble`](https://tibble.tidyverse.org/reference/tibble.html) data
frame

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
```
