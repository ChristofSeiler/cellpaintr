# Transform, center, and scale features

Transforms each family of CellProfiler features to an approximately
unbounded, symmetric scale (logarithm for intensities and
dispersion-type texture measures, Fisher's z for correlations, logit for
proportions, and a rank-based inverse normal transformation for
thresholded colocalization coefficients), caps infinite values at the
finite extremes of each feature, and standardizes every feature.

## Usage

``` r
transformScale(sce, robust = FALSE, assay_type = "features")
```

## Arguments

- sce:

  [`SingleCellExperiment`](https://rdrr.io/pkg/SingleCellExperiment/man/SingleCellExperiment.html)
  object

- robust:

  logical(1), if `TRUE` robust z-scores (median and median absolute
  deviation), otherwise standard z-scores

- assay_type:

  character(1), name of the assay with the features

## Value

[`SingleCellExperiment`](https://rdrr.io/pkg/SingleCellExperiment/man/SingleCellExperiment.html)
object with the transformed features in the new assay `"tfmfeatures"`

## Examples

``` r
set.seed(23)
cell_file <- generateData()
sce <- loadData(cell_file)
sce <- transformScale(sce)
```
