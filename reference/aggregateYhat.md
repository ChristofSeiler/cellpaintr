# Aggregate predicted leave-one-out probabilities over meta variables

Aggregate predicted leave-one-out probabilities over meta variables

## Usage

``` r
aggregateYhat(sce, target, group, assay_type = "tfmfeatures")
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

`data.frame`
