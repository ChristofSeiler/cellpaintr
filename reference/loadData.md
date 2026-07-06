# Load cell painting data from file and convert to a SingleCellExperiment

Load cell painting data from file and convert to a SingleCellExperiment

## Usage

``` r
loadData(cell_file)
```

## Arguments

- cell_file:

  path to csv file from CellProfiler, e.g., MyExpt_Cells.csv

## Value

[`SingleCellExperiment`](https://rdrr.io/pkg/SingleCellExperiment/man/SingleCellExperiment.html)
object

## Examples

``` r
set.seed(23)
cell_file <- generate_data()
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
