# Predict target from features

Predict target from features

## Usage

``` r
predictLOO(
  sce,
  target,
  group,
  interest_level,
  reference_level,
  types = NULL,
  channels = NULL,
  weights = NULL,
  n_threads = 1,
  assay_type = "tfmfeatures",
  progress = TRUE
)
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

- interest_level:

  character(1), level of \`target\` compared with the reference, such as
  a drug

- reference_level:

  character(1), level of \`target\` used as the reference, such as the
  control

- types:

  `NULL` or character vector of CellProfiler measurement types, such as
  `"Intensity"`; each type is a feature group of the features whose
  names start with it

- channels:

  `NULL` or character vector of imaging channels, such as `"DNA"` or
  `"Mito"`; each channel is a feature group of the features whose names
  contain it

- weights:

  `NULL` or character(1), name of the column of `colData(sce)` with case
  weights, for example cell counts when cells are aggregated to images

- n_threads:

  numeric(1), number of threads of each random forest fit

- assay_type:

  character(1), name of the assay with the features

- progress:

  logical(1), show a progress bar while fitting the models of feature
  types and channels

## Value

[`SingleCellExperiment`](https://rdrr.io/pkg/SingleCellExperiment/man/SingleCellExperiment.html)
object with the cells of the two levels, the leave-one-out prediction
scores in `reducedDim(sce, "prevalidated")`, and the settings of this
call in `metadata(sce)$predictLOO`, which
[`calculateStats`](https://christofseiler.github.io/cellpaintr/reference/calculateStats.md)
uses to refit the models in its permutation test

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
#>  ■■■■■■■■■■■■■■■■■■■■■■■           75% |  ETA:  1s
```
