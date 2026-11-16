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
  method = c("permutation", "corrected", "nadeau_bengio", "t"),
  n_permutations = 999
)
```

## Arguments

- sce:

  Output of
  [`predictLOO`](https://christofseiler.github.io/cellpaintr/reference/predictLOO.md)

- target:

  character(1), name of the column of `colData(sce)` with the condition
  of each cell, such as the drug

- group:

  character(1), name of the column of `colData(sce)` with the replicate
  of each cell, such as the patient; the replicates are the folds of the
  leave-one-out cross-validation

- assay_type:

  character(1), name of the assay with the features

- method:

  Either \`"permutation"\` (default) for a permutation test that
  relabels whole patients, \`"corrected"\` for the corrected paired
  t-test, \`"nadeau_bengio"\` for the paired t-test with the correction
  of Nadeau and Bengio, or \`"t"\` for p-values from the t-distribution

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
cells of both levels, and uses a pooled variance otherwise.

With \`method = "permutation"\` (the default), the p-value comes from
relabeling whole patients: in a paired design, the two levels of a
patient are swapped; in an unpaired design, the levels are permuted
across patients; and in a mixed design, where only some patients have
cells of both levels, the levels of those patients are swapped and the
levels of the other patients are permuted among them. For each
relabeling, all leave-one-out models are refitted with the settings of
the call to
[`predictLOO`](https://christofseiler.github.io/cellpaintr/reference/predictLOO.md),
and the p-value is the proportion of relabelings, the observed one
included, with a statistic at least as large. Exchanging the two levels
in all patients at once turns every score into one minus the score and
leaves the statistic unchanged, so relabelings come in pairs with the
same statistic. Whenever this exchange is itself a relabeling (always in
a paired design, and if the patients with one level are split evenly
between the two levels otherwise), one patient keeps its labels, which
halves the number of relabelings and refits without changing the test.
If there are at most \`n_permutations + 1\` relabelings (\\2^{n-1}\\ for
\\n\\ patients in a paired design), all of them are enumerated and the
test is exact; otherwise, \`n_permutations\` relabelings are drawn at
random. The smallest attainable p-value is one over the number of
relabelings, \\2^{1-n}\\ in a paired design. Because all models are
refitted for every relabeling, the test takes that many times longer
than
[`predictLOO`](https://christofseiler.github.io/cellpaintr/reference/predictLOO.md).
Use [`set.seed`](https://rdrr.io/r/base/Random.html) for reproducible
results.

With \`method = "corrected"\`, the p-value comes from a corrected paired
t-test, which needs a paired design with at least three patients. Each
patient is in the training set of every other patient's model, which
makes the score differences of different patients positively correlated,
so the t-statistic is too large. If each score difference is a sum of
pairwise interactions with the other patients, two score differences
have correlation \\1/(n-1)\\ for \\n\\ patients, and the t-statistic is
too large by the factor \\\sqrt{2(n-1)/(n-2)}\\. The corrected t-test
divides the t-statistic by this factor and refers it to the
t-distribution with \\n-1\\ degrees of freedom. It needs no refits, and
its p-values are not limited to multiples of \\2^{1-n}\\, but it is
approximate.

With \`method = "nadeau_bengio"\`, the t-statistic is divided instead by
\\\sqrt{(2n-1)/(n-1)}\\, the factor of the corrected resampled t-test of
Nadeau and Bengio (2003, Machine Learning 52:239-281) for \\n\\ splits
with one held-out patient each, which assumes a correlation of \\1/n\\
instead of \\1/(n-1)\\. It corrects less and therefore has more power,
which pays off if part of each score difference is independent of the
other patients, for example noise from the random seeds of the forests.
If the score differences consist mostly of pairwise interactions, it can
exceed the nominal level with four to six patients, where \`method =
"corrected"\` holds it. It needs a paired design with at least two
patients and no refits.

With \`method = "t"\`, the p-value comes from the t-distribution. This
is fast but anti-conservative, because the score differences of
different patients are positively correlated.

## Examples

``` r
set.seed(23)
cell_file <- generateData()
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
#> 'P1' keeps its labels, because exchanging the two levels of 'Drug' in all levels of 'Patient' at once gives the same statistic.
#>   Target Feature pvalue log2FoldChange
#> 1     D7     all    0.1      0.9328005

# corrected t-test: no refits
calculateStats(sce_single,
    target = "Drug", group = "Patient",
    method = "corrected"
)
#>   Target Feature       pvalue log2FoldChange
#> 1     D7     all 0.0001791293      0.9328005

# correction of Nadeau and Bengio: no refits, corrects less
calculateStats(sce_single,
    target = "Drug", group = "Patient",
    method = "nadeau_bengio"
)
#>   Target Feature       pvalue log2FoldChange
#> 1     D7     all 0.0001322961      0.9328005

# p-values from the t-distribution: fast but anti-conservative
calculateStats(sce_single,
    target = "Drug", group = "Patient",
    method = "t"
)
#>   Target Feature       pvalue log2FoldChange
#> 1     D7     all 1.971207e-05      0.9328005
```
