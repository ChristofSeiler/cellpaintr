# Changelog

## cellpaintr 0.1.0

- Experimental version still under development

## cellpaintr 0.2.0

- Version used internally in our lab

## cellpaintr 0.3.0

- Getting it ready for submission to Bioconductor

## cellpaintr 0.99.0

- Submitted to Bioconductor

## cellpaintr 0.99.4

- [`calculateStats()`](https://christofseiler.github.io/cellpaintr/reference/calculateStats.md)
  calibrates the t-statistic with a permutation test that relabels whole
  patients and refits all leave-one-out models
  (`method = "permutation"`, the new default). It enumerates all
  relabelings if there are at most `n_permutations + 1` of them and
  samples them otherwise. `method = "t"` keeps the p-values from the
  t-distribution, which are anti-conservative because the scores of
  different held-out patients are correlated.
- [`predictLOO()`](https://christofseiler.github.io/cellpaintr/reference/predictLOO.md)
  records its settings in `metadata()`, and
  [`calculateStats()`](https://christofseiler.github.io/cellpaintr/reference/calculateStats.md)
  takes the interest and reference levels from there instead of from the
  order of the factor levels.
- [`volcanoPlot()`](https://christofseiler.github.io/cellpaintr/reference/volcanoPlot.md)
  accepts precomputed statistics (`stats`) and passes further arguments,
  such as `method`, to
  [`calculateStats()`](https://christofseiler.github.io/cellpaintr/reference/calculateStats.md).

## cellpaintr 0.99.5

- The permutation test of
  [`calculateStats()`](https://christofseiler.github.io/cellpaintr/reference/calculateStats.md)
  handles mixed designs, in which only some patients have cells of both
  levels, for example after quality control removed all cells of one
  level of a patient. It swaps the levels of the patients with both
  levels, permutes the levels among the other patients, and reports
  which patients have cells of one level only.

## cellpaintr 0.99.6

- [`calculateStats()`](https://christofseiler.github.io/cellpaintr/reference/calculateStats.md)
  offers a corrected paired t-test (`method = "corrected"`). For every
  pair of patients, it refits the models without both of them and adds
  the products of their influences on each other’s score differences to
  the sample variance, which makes the variance unbiased under the null
  hypothesis. It needs n(n-1)/2 fits per feature group for n patients
  instead of refitting all models for every relabeling, and it can reach
  small p-values with few patients. The attribute `variance_ratio` holds
  the corrected over the uncorrected variance. The permutation test
  stays the default.
- Internal: model fitting of `compute_y_hat()` moved into helpers shared
  with the corrected t-test; results are unchanged.

## cellpaintr 0.99.7

Changes in response to the first Bioconductor review:

- `generate_data()` is renamed
  [`generateData()`](https://christofseiler.github.io/cellpaintr/reference/generateData.md),
  in line with the rest of the API.
- New example data: `inst/extdata/jump_subset.csv.gz`, 240 cells of the
  JUMP Cell Painting dataset (DMSO, AMG900, and dexamethasone on four
  plates) as exported by CellProfiler; `inst/scripts/make-jump-subset.R`
  documents how it was made, and `inst/scripts/README.md` the provenance
  of all files in `inst/extdata`.
- The vignette is rewritten around the JUMP data: it describes the
  CellProfiler input and its `Metadata_` columns, explains every step,
  and shows the messages of
  [`calculateStats()`](https://christofseiler.github.io/cellpaintr/reference/calculateStats.md).
- New package man page
  ([`?cellpaintr`](https://christofseiler.github.io/cellpaintr/reference/cellpaintr-package.md))
  with an overview of the workflow.
- Exported functions validate their inputs and stop with informative
  messages;
  [`predictLOO()`](https://christofseiler.github.io/cellpaintr/reference/predictLOO.md)
  accepts a character `target`, needs at least two levels of `group`,
  and keeps the cell names of `sce` instead of replacing them with
  integers.
- [`loadData()`](https://christofseiler.github.io/cellpaintr/reference/loadData.md)
  passes further arguments to
  [`readr::read_csv()`](https://readr.tidyverse.org/reference/read_delim.html),
  for example `num_threads`.
- [`removeNAs()`](https://christofseiler.github.io/cellpaintr/reference/removeNAs.md),
  [`removeLowVariance()`](https://christofseiler.github.io/cellpaintr/reference/removeLowVariance.md),
  [`removeZeroInflation()`](https://christofseiler.github.io/cellpaintr/reference/removeZeroInflation.md),
  and
  [`transformScale()`](https://christofseiler.github.io/cellpaintr/reference/transformScale.md)
  have an `assay_type` argument, and use vectorized
  [`colSums()`](https://rdrr.io/r/base/colSums.html),
  [`rowMeans()`](https://rdrr.io/r/base/colSums.html), and
  `MatrixGenerics` instead of
  [`apply()`](https://rdrr.io/r/base/apply.html).
- The relabelings of the permutation test run through
  [`future.apply::future_lapply()`](https://future.apply.futureverse.org/reference/future_lapply.html),
  sequentially by default and in parallel after
  [`future::plan()`](https://future.futureverse.org/reference/plan.html),
  with reproducible random numbers.
- Functions are imported instead of declared as global variables;
  `rlang`, `stats`, `methods`, `MatrixGenerics`, and `future.apply` are
  imported.
- Long functions are split into helpers, and the code follows the
  Bioconductor style (80 characters, 4-space indents).
- `compute_y_hat()` and `aggregateYhat()` are internal and have no man
  pages.
- More unit tests check values, not only classes, and cover the input
  checks.

Changes to improve the package but unrelated to the Bioconductor review:

- The corrected t-test of
  [`calculateStats()`](https://christofseiler.github.io/cellpaintr/reference/calculateStats.md)
  (`method = "corrected"`) now divides the paired t-statistic by
  sqrt(2(n-1)/(n-2)), the factor by which it is too large if each score
  difference is a sum of pairwise interactions with the other patients.
  It needs no refits and no longer needs the settings of
  [`predictLOO()`](https://christofseiler.github.io/cellpaintr/reference/predictLOO.md).
  It replaces the leave-two-out variance correction of 0.99.6 and its
  attribute `variance_ratio`.
- [`calculateStats()`](https://christofseiler.github.io/cellpaintr/reference/calculateStats.md)
  offers the correction of the corrected resampled t-test of Nadeau and
  Bengio (2003) (`method = "nadeau_bengio"`), which divides the paired
  t-statistic by sqrt((2n-1)/(n-1)). It corrects less than
  `method = "corrected"` and needs a paired design with at least two
  patients and no refits.
- The permutation test of
  [`calculateStats()`](https://christofseiler.github.io/cellpaintr/reference/calculateStats.md)
  no longer refits relabelings that only repeat another one. Exchanging
  the two levels in all patients at once turns every score into one
  minus the score and leaves the statistic unchanged, so one patient now
  keeps its labels whenever this exchange is itself a relabeling (always
  in a paired design). This halves the number of refits, and the test
  stays exact. In a paired design with n patients, there are 2^(n-1)
  relabelings, and the smallest attainable p-value is 2^(1-n) instead of
  a nominal 2^-n that depended on the forest seeds.
- The vignette explains that ROC curves and AUCs use the mean scores of
  each replicate and condition.
- The vignette scales the well means before the PCA. Robust z-scores are
  extreme for features with a median absolute deviation close to zero,
  and one well dominated the first principal component.
- The DESCRIPTION and the vignette describe when cellpaintr is more
  powerful than a location test on patient means (most cells shift, or
  the spread of features changes) and when it is less powerful (a strong
  response of a small subpopulation), following a comparison on effects
  spiked into JUMP control wells.

Bug fixes:

- [`plotROC()`](https://christofseiler.github.io/cellpaintr/reference/plotROC.md)
  and
  [`plotAUC()`](https://christofseiler.github.io/cellpaintr/reference/plotAUC.md)
  take the interest and reference levels from the settings of
  [`predictLOO()`](https://christofseiler.github.io/cellpaintr/reference/predictLOO.md)
  instead of from the order of the factor levels of `target`. Before, a
  character `target` whose interest level sorted first alphabetically,
  such as AMG900 against DMSO, swapped the two levels and turned every
  AUC into one minus the AUC.
- [`predictLOO()`](https://christofseiler.github.io/cellpaintr/reference/predictLOO.md)
  orders the levels of a character `target` with the reference level
  first.
- [`removeRowColEffect()`](https://christofseiler.github.io/cellpaintr/reference/removeRowColEffect.md)
  subtracted the row and column effects of the wrong rows or columns
  when the first column of a plate lacked a row or the first row lacked
  a column, because it matched the effects to cells by position. It now
  matches them by name, and its median polish leaves out empty wells
  instead of setting their medians to zero.

## cellpaintr 0.99.8

Changes in response to the second Bioconductor review:

- [`removeRowColEffect()`](https://christofseiler.github.io/cellpaintr/reference/removeRowColEffect.md)
  has a `plate` argument: with the column of the plate or batch of each
  cell, it fits the median polish to each plate separately. Without it,
  all cells are polished as one plate, and the function warns if
  `colData(sce)` has a `Plate` column with several plates.
- The vignette’s example of parallel permutations uses two workers.
- The vignette and
  [`?calculateStats`](https://christofseiler.github.io/cellpaintr/reference/calculateStats.md)
  cite the paper that derives the corrected t-tests and studies their
  type I error and power with 4 to 12 replicates; the preprint will be
  available soon.
- The biocViews include DataImport, CellBasedAssays, and
  StatisticalMethod.
