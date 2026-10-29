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
