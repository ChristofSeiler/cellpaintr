test_that("train prediction model", {
    set.seed(23)
    cell_file <- generate_data()
    sce <- loadData(cell_file)
    sce <- transformLogScale(sce)

    sce$Drug <- as.factor(sce$Drug)
    sce$Drug <- relevel(sce$Drug, ref = "D1")
    types <- c("AreaShape", "Intensity", "Texture")

    # leave-one-out score
    sce_single <- predictLOO(
        sce,
        target = "Drug", group = "Patient",
        interest_level = "D7", reference_level = "D1",
        types = types,
        n_threads = 1
    )
    scores <- reducedDim(sce_single, "prevalidated")

    expect_setequal(colnames(scores), c("all", types))
    expect_gte(min(scores), 0)
    expect_lte(max(scores), 1)
})

test_that("extract prediction scores", {
    set.seed(23)
    cell_file <- generate_data()
    sce <- loadData(cell_file)
    sce <- transformLogScale(sce)

    sce$Drug <- as.factor(sce$Drug)
    sce$Drug <- relevel(sce$Drug, ref = "D1")
    types <- c("AreaShape", "Intensity", "Texture")

    # leave-one-out score
    sce_single <- predictLOO(
        sce,
        target = "Drug", group = "Patient",
        interest_level = "D7", reference_level = "D1",
        types = types,
        n_threads = 1
    )
    stats <- calculateStats(sce_single,
        meta_vars = c("Patient", "Drug"),
        target = "Drug"
    )

    expect_setequal(
        colnames(stats),
        c("Target", "Feature", "pvalue", "log2FoldChange")
    )
    expect_equal(nrow(stats), 4)
    expect_equal(ncol(stats), 4)
    expect_gte(min(stats$pvalue), 0)
    expect_lte(max(stats$pvalue), 1)
})
