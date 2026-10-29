# Test each feature group for a difference between two levels

Average the leave-one-out prediction scores of
[`predictLOO`](https://christofseiler.github.io/cellpaintr/reference/predictLOO.md)
per patient and level, and test in each feature group whether the
interest level has higher scores than the reference level.

## Usage

``` r
calculateStats(
  sce,
  target,
  group,
  assay_type = "tfmfeatures",
  method = c("permutation", "t"),
  n_permutations = 999
)
```

## Arguments

- sce:

  Output of
  [`predictLOO`](https://christofseiler.github.io/cellpaintr/reference/predictLOO.md)

- target:

  Name of target variable for prediction

- group:

  Grouping variable for cross-validation, e.g., patient

- assay_type:

  A string specifying the assay

- method:

  Either \`"permutation"\` (default) for a permutation test that
  relabels whole patients, or \`"t"\` for p-values from the
  t-distribution

- n_permutations:

  Largest number of relabelings besides the observed labels; if there
  are more, this many are drawn at random

## Value

`data.frame` with one row per feature group and the columns \`Target\`
(interest level), \`Feature\` (feature group), \`pvalue\`, and
\`log2FoldChange\` (log2 ratio of the mean scores of the interest and
the reference level). With the permutation test, the attribute
\`relabelings\` holds the number of relabelings, the observed one
included.

## Details

The test statistic is the one-sided t-statistic for higher mean scores
in the interest level. It is paired if every patient (\`group\`) has
cells of both levels, and uses a pooled variance if every patient has
cells of one level.

With \`method = "permutation"\` (the default), the p-value comes from
relabeling whole patients: in a paired design, the two levels of a
patient are swapped; in an unpaired design, the levels are permuted
across patients. For each relabeling, all leave-one-out models are
refitted with the settings of the call to
[`predictLOO`](https://christofseiler.github.io/cellpaintr/reference/predictLOO.md),
and the p-value is the proportion of relabelings, the observed one
included, with a statistic at least as large. If there are at most
\`n_permutations + 1\` relabelings (\\2^n\\ for \\n\\ patients in a
paired design), all of them are enumerated and the test is exact;
otherwise, \`n_permutations\` relabelings are drawn at random. The
smallest attainable p-value is one over the number of relabelings.
Because all models are refitted for every relabeling, the test takes
that many times longer than
[`predictLOO`](https://christofseiler.github.io/cellpaintr/reference/predictLOO.md).
Use [`set.seed`](https://rdrr.io/r/base/Random.html) for reproducible
results.

With \`method = "t"\`, the p-value comes from the t-distribution. This
is fast but anti-conservative, because the models of different held-out
patients share most of their training data, which makes the scores of
different patients positively correlated.

## Examples

``` r
set.seed(23)
cell_file <- generate_data()
sce <- loadData(cell_file)
sce <- transformScale(sce)

sce$Drug <- as.factor(sce$Drug)
sce$Drug <- relevel(sce$Drug, ref = "D1")

sce_single <- predictLOO(
    sce,
    target = "Drug", group = "Patient",
    interest_level = "D7", reference_level = "D1",
    n_threads = 1
)

# permutation test; by default, it enumerates all 64 relabelings of the
# six patients, here 9 random relabelings keep the example fast
calculateStats(sce_single,
    target = "Drug", group = "Patient",
    n_permutations = 9
)
#> Permutation test with 10 relabelings of 6 levels of 'Patient' (paired design); smallest attainable p-value 1/10.
#> Refitting models ■■■■■■■■                          22% | ETA:  5s
#> Refitting models ■■■■■■■■■■■■■■■■■■                56% | ETA:  3s
#> Refitting models ■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■  100% | ETA:  0s
#>   Target Feature pvalue log2FoldChange
#> 1     D7     all    0.1      0.9304184

# p-values from the t-distribution: fast but anti-conservative
calculateStats(sce_single,
    target = "Drug", group = "Patient",
    method = "t"
)
#>   Target Feature       pvalue log2FoldChange
#> 1     D7     all 3.051405e-05      0.9304184
```
