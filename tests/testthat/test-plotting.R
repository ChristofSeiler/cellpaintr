test_that("plot cells per image", {
    set.seed(23)
    cell_file <- generate_data()
    sce <- loadData(cell_file)
    p <- plotCellsPerImage(sce)

    expect_s3_class(p, "ggplot")
    expect_equal(rlang::as_name(p$mapping$x), "n")
})

test_that("plot pca cor", {
    set.seed(23)
    cell_file <- generate_data()
    sce <- loadData(cell_file)
    sce <- transformLogScale(sce, robust = TRUE)
    sce <- scater::runPCA(sce, exprs_values = "tfmfeatures", ncomponents = 10)
    p <- plotPCACor(sce, filter_by = 1)

    expect_s3_class(p, "ggplot")
    expect_equal(rlang::as_name(p$mapping$x), "PC")
})

test_that("plot prediction scores", {
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
    p <- plotLOO(sce_single, target = "Drug", group = "Patient")
    expect_s3_class(p, "ggplot")
    expect_equal(rlang::as_name(p$mapping$y), "value")

    # volcano plot
    p <- volcanoPlot(sce_single,
        target = "Drug", group = "Patient",
        p_cutoff = 0.05, fc_cutoff = 0.5
    )
    expect_s3_class(p, "ggplot")
    expect_equal(rlang::as_name(p$mapping$x), "log2FoldChange")

    # roc curve
    p <- plotROC(sce_single, target = "Drug", group = "Patient")
    expect_s3_class(p, "ggplot")
    expect_equal(rlang::as_name(p$mapping$y), "sensitivity")

    # auc summary
    p <- plotAUC(sce_single, target = "Drug", group = "Patient")
    expect_s3_class(p, "ggplot")
    expect_equal(rlang::as_name(p$mapping$x), ".estimate")

    # no p_cutoff defined
    p <- volcanoPlot(sce_single,
        target = "Drug", group = "Patient",
        fc_cutoff = 0.5
    )
    expect_s3_class(p, "ggplot")
})
