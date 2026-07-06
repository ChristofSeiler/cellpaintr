test_that("remove missing values", {
    set.seed(23)
    cell_file <- generate_data()
    sce <- loadData(cell_file)

    assay(sce)[1, 1] <- NA
    sce <- removeNAs(sce)
    n_nas <- sum(is.na(assay(sce)))

    expect_equal(n_nas, 0)
})

test_that("remove outliers", {
    set.seed(23)
    cell_file <- generate_data()
    sce <- loadData(cell_file)

    sce <- removeOutliers(sce, min = 0, max = 17)

    expect_equal(sum(table(sce$ImageNumber) > 17), 0)
})

test_that("remove low variance features (robust = FALSE)", {
    set.seed(23)
    cell_file <- generate_data()
    sce <- loadData(cell_file)

    assay(sce)["AreaShape_Zernike_2_2", ] <- 0
    sce <- removeLowVariance(sce, robust = FALSE)

    expect_disjoint("AreaShape_Zernike_2_2", rownames(sce))
})

test_that("remove low variance features (robust = TRUE)", {
    set.seed(23)
    cell_file <- generate_data()
    sce <- loadData(cell_file)

    assay(sce)["AreaShape_Zernike_2_2", ] <- 0
    sce <- removeLowVariance(sce, robust = TRUE)

    expect_disjoint("AreaShape_Zernike_2_2", rownames(sce))
})

test_that("remove zero inflation features", {
    set.seed(23)
    cell_file <- generate_data()
    sce <- loadData(cell_file)

    assay(sce)["AreaShape_Zernike_2_2", ] <- 0
    sce <- removeZeroInflation(sce)

    expect_disjoint("AreaShape_Zernike_2_2", rownames(sce))
})
