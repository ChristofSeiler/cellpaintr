# Internal function to compute y_hat on a subset of the features

Internal function to compute y_hat on a subset of the features

## Usage

``` r
compute_y_hat(
  feature_name,
  sce_feature,
  starts,
  assay_type,
  target,
  interest_level,
  reference_level,
  group,
  weights,
  n_threads
)
```

## Arguments

- feature_name:

  String with feature name

- sce_feature:

  [`SingleCellExperiment`](https://rdrr.io/pkg/SingleCellExperiment/man/SingleCellExperiment.html)

- starts:

  Starting string

- assay_type:

  A string specifying the assay

- target:

  Name of target variable for prediction

- interest_level:

  Factor interest level in \`target\` variable

- reference_level:

  Factor reference level in \`target\` variable

- group:

  Grouping variable for cross-validation, e.g., patient

- weights:

  Weights variable when features are aggregated

- n_threads:

  Number of parallel threads for fitting of models

## Value

[`tibble`](https://tibble.tidyverse.org/reference/tibble.html) data
frame
