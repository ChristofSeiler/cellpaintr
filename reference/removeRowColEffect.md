# Remove row and column effect using the median polish procedure

Remove row and column effect using the median polish procedure

## Usage

``` r
removeRowColEffect(sce, assay_type = "tfmfeatures")
```

## Arguments

- sce:

  [`SingleCellExperiment`](https://rdrr.io/pkg/SingleCellExperiment/man/SingleCellExperiment.html)
  object

- robust:

  If true robust z-score, otherwise standard z-score

## Value

[`SingleCellExperiment`](https://rdrr.io/pkg/SingleCellExperiment/man/SingleCellExperiment.html)
object

## Examples

``` r
set.seed(23)
cell_file <- generate_data()
sce <- loadData(cell_file)
sce <- transformScale(sce)
row <- substr(sce$Well, 1, 1)
row <- factor(row, levels = LETTERS[1:length(unique(row))])
col <- substr(sce$Well, 2, 3)
col <- factor(as.numeric(col))
sce$Row <- row
sce$Col <- col
sce <- removeRowColEffect(sce)
#> Error in join_by(Row): could not find function "join_by"
```
