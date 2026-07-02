test_that("remove missing values", {
    set.seed(23)
    cell_file <- generate_data()
    sce <- loadData(cell_file)

    assay(sce)[1, 1] <- NA
    sce <- removeNAs(sce)
    n_nas <- sum(is.na(assay(sce)))

    expect_equal(n_nas, 0)
})
