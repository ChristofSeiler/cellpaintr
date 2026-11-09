# Changelog

## cellpaintr 0.1.0

- Experimental version still under development

## cellpaintr 0.2.0

- Version used internally in our lab

## cellpaintr 0.3.0

- Getting it ready for submission to Bioconductor

## cellpaintr 0.99.0

- Submitted to Bioconductor

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
- Internal: model fitting of
  [`compute_y_hat()`](https://christofseiler.github.io/cellpaintr/reference/compute_y_hat.md)
  moved into helpers shared with the corrected t-test; results are
  unchanged.

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
