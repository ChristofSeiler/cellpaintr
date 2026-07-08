test_that("train prediction model on types", {
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

test_that("train prediction model on types with weights", {
    set.seed(23)
    cell_file <- generate_data()
    sce <- loadData(cell_file)
    sce <- transformLogScale(sce)

    # aggregate
    sce_aggr <- aggregateAcrossCells(
        sce,
        id = colData(sce)[, c("ImageNumber", "Patient")],
        use.assay.type = "tfmfeatures",
        statistics = "mean"
    )

    sce_aggr$Drug <- as.factor(sce_aggr$Drug)
    sce_aggr$Drug <- relevel(sce_aggr$Drug, ref = "D1")
    types <- c("AreaShape", "Intensity", "Texture")

    # leave-one-out score
    sce_single <- predictLOO(
        sce_aggr,
        weights = "ncells",
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

test_that("train prediction model on channels", {
    set.seed(23)
    cell_file <- generate_data()
    sce <- loadData(cell_file)
    sce <- transformLogScale(sce)

    sce$Drug <- as.factor(sce$Drug)
    sce$Drug <- relevel(sce$Drug, ref = "D1")
    channels <- c("CorrBlue", "CorrTMRM", "CorrLyso", "CorrCytoID")

    # leave-one-out score
    sce_single <- predictLOO(
        sce,
        target = "Drug", group = "Patient",
        interest_level = "D7", reference_level = "D1",
        channels = channels,
        n_threads = 1
    )
    scores <- reducedDim(sce_single, "prevalidated")

    expect_setequal(colnames(scores), c("all", channels))
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
        target = "Drug", group = "Patient"
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
