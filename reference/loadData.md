# Load cell painting data from file and convert to a SingleCellExperiment

Reads a per-cell table exported by CellProfiler's ExportToSpreadsheet
module, for example `MyExpt_Cells.csv`, with one row per cell. The table
needs the column `ImageNumber`, at least one image metadata column whose
name starts with `Metadata_` (CellProfiler adds them with the option
"Add image metadata columns to your object data file"), and numeric
measurement columns such as `AreaShape_Area` or
`Intensity_MeanIntensity_DNA`. Metadata columns become the column data
of the cells, without the `Metadata_` prefix, together with
`ImageNumber`, `ObjectNumber`, and `Location_Center_X`, `_Y`, and `_Z`;
all other columns become the features in the assay `"features"`.

## Usage

``` r
loadData(cell_file, ...)
```

## Arguments

- cell_file:

  character(1), path to a csv file from CellProfiler, e.g.,
  `MyExpt_Cells.csv`; compressed files such as `.csv.gz` work

- ...:

  Further arguments to
  [`read_csv`](https://readr.tidyverse.org/reference/read_delim.html),
  for example `num_threads` to limit the threads used to read large
  files

## Value

[`SingleCellExperiment`](https://rdrr.io/pkg/SingleCellExperiment/man/SingleCellExperiment.html)
object with features as rows and cells as columns

## Examples

``` r
set.seed(23)
cell_file <- generateData()
loadData(cell_file)
#> class: SingleCellExperiment 
#> dim: 30 1000 
#> metadata(0):
#> assays(1): features
#> rownames(30): AreaShape_Zernike_2_2 AreaShape_Zernike_8_4 ...
#>   Texture_InverseDifferenceMoment_CorrBlue_5_00_256
#>   Texture_InfoMeas1_CorrBlue_3_01_256
#> rowData names(0):
#> colnames: NULL
#> colData names(8): ObjectNumber ImageNumber ... Location_Center_Y
#>   Location_Center_Z
#> reducedDimNames(0):
#> mainExpName: NULL
#> altExpNames(0):
```
