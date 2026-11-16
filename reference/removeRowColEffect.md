# Remove row and column effect using the median polish procedure

Fits Tukey's median polish to the per-well medians of each feature and
subtracts the overall, row, and column effects from every cell, similar
to the B-score of high-throughput screening. The column data need the
factors `Row` and `Col` with the plate position of each cell.

## Usage

``` r
removeRowColEffect(sce, assay_type = "tfmfeatures")
```

## Arguments

- sce:

  [`SingleCellExperiment`](https://rdrr.io/pkg/SingleCellExperiment/man/SingleCellExperiment.html)
  object

- assay_type:

  character(1), name of the assay with the features

## Value

[`SingleCellExperiment`](https://rdrr.io/pkg/SingleCellExperiment/man/SingleCellExperiment.html)
object with the corrected features in the new assay `"corrected"`

## Examples

``` r
set.seed(23)
cell_file <- generateData()
sce <- loadData(cell_file)
sce <- transformScale(sce)
row <- substr(sce$Well, 1, 1)
row <- factor(row, levels = sort(unique(row)))
col <- substr(sce$Well, 2, 3)
col <- factor(as.numeric(col))
sce$Row <- row
sce$Col <- col
sce <- removeRowColEffect(sce)
```
