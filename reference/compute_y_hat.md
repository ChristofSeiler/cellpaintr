# Internal function to compute y_hat on a subset of the features

Internal function to compute y_hat on a subset of the features

## Usage

``` r
compute_y_hat(
  feature_name,
  sce_feature,
  starts,
  target,
  group,
  interest_level,
  reference_level,
  weights,
  n_threads,
  assay_type
)
```

## Arguments

- feature_name:

  String with feature name

- sce_feature:

  [`SingleCellExperiment`](https://rdrr.io/pkg/SingleCellExperiment/man/SingleCellExperiment.html)

- starts:

  Starting string

- target:

  Name of target variable for prediction

- group:

  Grouping variable for cross-validation, e.g., patient

- interest_level:

  Factor interest level in \`target\` variable

- reference_level:

  Factor reference level in \`target\` variable

- weights:

  Weights variable when features are aggregated

- n_threads:

  Number of parallel threads for fitting of models

- assay_type:

  A string specifying the assay

## Value

[`tibble`](https://tibble.tidyverse.org/reference/tibble.html) data
frame
