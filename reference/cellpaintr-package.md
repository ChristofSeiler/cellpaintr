# cellpaintr: perturbation analysis of Cell Painting data with few replicates

cellpaintr tests whether a perturbation, such as a drug, changes the
morphology of cells reproducibly across a handful of replicates, such as
patients or plates. It learns on single cells and tests on replicates: a
random forest trained on the cells of all other replicates scores the
cells of each held-out replicate (pre-validation), and a patient-level
test compares the scores between the two conditions, for all features
and for interpretable feature groups such as measurement types or
imaging channels.

## Workflow

1.  Load a per-cell CellProfiler table with
    [`loadData`](https://christofseiler.github.io/cellpaintr/reference/loadData.md).

2.  Remove poor cells and features with
    [`removeOutliers`](https://christofseiler.github.io/cellpaintr/reference/removeOutliers.md),
    [`removeNAs`](https://christofseiler.github.io/cellpaintr/reference/removeNAs.md),
    [`removeLowVariance`](https://christofseiler.github.io/cellpaintr/reference/removeLowVariance.md),
    and
    [`removeZeroInflation`](https://christofseiler.github.io/cellpaintr/reference/removeZeroInflation.md).

3.  Transform and scale the features with
    [`transformScale`](https://christofseiler.github.io/cellpaintr/reference/transformScale.md),
    and optionally remove plate row and column effects with
    [`removeRowColEffect`](https://christofseiler.github.io/cellpaintr/reference/removeRowColEffect.md).

4.  Score every cell with a classifier that never saw its replicate,
    with
    [`predictLOO`](https://christofseiler.github.io/cellpaintr/reference/predictLOO.md).

5.  Test each feature group with
    [`calculateStats`](https://christofseiler.github.io/cellpaintr/reference/calculateStats.md):
    an exact permutation test that refits all models for every
    relabeling of the replicates, or a corrected paired t-test that
    needs no refits.

6.  Summarize the results with
    [`volcanoPlot`](https://christofseiler.github.io/cellpaintr/reference/volcanoPlot.md),
    [`plotLOO`](https://christofseiler.github.io/cellpaintr/reference/plotLOO.md),
    [`plotROC`](https://christofseiler.github.io/cellpaintr/reference/plotROC.md),
    and
    [`plotAUC`](https://christofseiler.github.io/cellpaintr/reference/plotAUC.md).

[`generateData`](https://christofseiler.github.io/cellpaintr/reference/generateData.md)
simulates a small CellProfiler table, and
`system.file("extdata", "jump_subset.csv.gz", package = "cellpaintr")`
holds a small subset of the JUMP Cell Painting data. The vignette walks
through a complete analysis:
[`vignette("cellpaintr-vignette", package = "cellpaintr")`](https://christofseiler.github.io/cellpaintr/articles/cellpaintr-vignette.md).

## See also

Useful links:

- <https://github.com/ChristofSeiler/cellpaintr>

- <https://christofseiler.github.io/cellpaintr/>

- Report bugs at <https://github.com/ChristofSeiler/cellpaintr/issues>

## Author

**Maintainer**: Christof Seiler <christof.seiler@usz.ch>
([ORCID](https://orcid.org/0000-0001-8802-3642))

Other contributors:

- Phatthamon Laphanuwat <phatthamon.laphanuwat@usz.ch> \[contributor\]

- Caroline Ospelt <caroline.ospelt@usz.ch> \[contributor\]

- European League Against Rheumatism \[funder\]

- Swiss National Science Foundation \[funder\]
