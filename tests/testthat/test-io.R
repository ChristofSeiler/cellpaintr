test_that("simulate data", {
    cell_file <- generate_data()
    expect_true(file.exists(cell_file))
})

test_that("simulate data", {
    cell_file <- generate_data()
    sce <- loadData(cell_file)
    expect_true(is(sce, "SingleCellExperiment"))
})
