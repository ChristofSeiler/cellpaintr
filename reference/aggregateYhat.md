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
