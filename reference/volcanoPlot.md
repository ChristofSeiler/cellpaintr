# Plot predicted leave-one-out probabilities

Plot predicted leave-one-out probabilities

## Usage

``` r
volcanoPlot(
  sce,
  target,
  group,
  p_cutoff = NULL,
  fc_cutoff = 1,
  assay_type = "tfmfeatures",
  stats = NULL,
  ...
)
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

- p_cutoff:

  Cut-off for statistical significance. A horizontal line will be drawn
  at -log10(p_cutoff). Defaults to 0.01 divided by the number of feature
  groups. With few patients, permutation p-values cannot fall below this
  default, because the smallest attainable p-value is one over the
  number of relabelings.

- fc_cutoff:

  Cut-off for absolute log2 fold-change. A vertical lines will be drawn
  at fc_cutoff.

- assay_type:

  A string specifying the assay

- stats:

  Output of
  [`calculateStats`](https://christofseiler.github.io/cellpaintr/reference/calculateStats.md);
  if \`NULL\`, it is computed from \`sce\`

- ...:

  Further arguments to
  [`calculateStats`](https://christofseiler.github.io/cellpaintr/reference/calculateStats.md),
  such as \`method\` or \`n_permutations\`

## Value

[`ggplot2`](https://ggplot2.tidyverse.org/reference/ggplot2-package.html)
object

## Examples

``` r
set.seed(23)
cell_file <- generate_data()
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

# p-values from the t-distribution keep the example fast; the default
# permutation test refits all models for every relabeling of the patients
volcanoPlot(sce_single,
    target = "Drug", group = "Patient",
    p_cutoff = 0.05, fc_cutoff = 0.5, method = "t"
)

```
