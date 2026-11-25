# Remove row and column effect using the median polish procedure

Fits Tukey's median polish to the per-well medians of each feature and
subtracts the overall, row, and column effects from every cell, similar
to the B-score of high-throughput screening. The column data need the
factors `Row` and `Col` with the plate position of each cell. Row and
column effects differ between plates, so with several plates, pass the
column with the plate (or batch) of each cell as \`plate\`; the median
polish is then fitted to each plate separately. Wells without cells are
left out of the fit. Because the median polish removes everything that
varies with rows or columns, it also removes the part of a treatment
effect that is confounded with them, so the conditions should be spread
over the rows and columns of each plate.

## Usage

``` r
removeRowColEffect(sce, assay_type = "tfmfeatures", plate = NULL)
```

## Arguments

- sce:

  [`SingleCellExperiment`](https://rdrr.io/pkg/SingleCellExperiment/man/SingleCellExperiment.html)
  object

- assay_type:

  character(1), name of the assay with the features

- plate:

  `NULL` or character(1), name of the column of `colData(sce)` with the
  plate or batch of each cell; if `NULL`, all cells come from one plate

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

# several plates: one median polish per plate
sce$Plate <- ifelse(sce$Patient %in% c("P1", "P2", "P3"),
    "plate1", "plate2"
)
sce <- removeRowColEffect(sce, plate = "Plate")
```
