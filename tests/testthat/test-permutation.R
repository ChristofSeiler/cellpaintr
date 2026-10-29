# small designs keep the number of refits low
simulate_sce <- function() {
    set.seed(23)
    cell_file <- generate_data()
    sce <- loadData(cell_file)
    sce <- transformScale(sce)
    sce$Drug <- relevel(as.factor(sce$Drug), ref = "D1")
    sce
}

fit_loo <- function(sce) {
    predictLOO(sce,
        target = "Drug", group = "Patient",
        interest_level = "D8", reference_level = "D1",
        n_threads = 1, progress = FALSE
    )
}

on_grid <- function(pvalue, n) all(abs(pvalue * n - round(pvalue * n)) < 1e-8)

test_that("predictLOO records its settings", {
    sce <- simulate_sce()
    sce_single <- fit_loo(sce[, sce$Patient %in% c("P1", "P2", "P3")])
    settings <- S4Vectors::metadata(sce_single)$predictLOO

    expect_equal(settings$target, "Drug")
    expect_equal(settings$group, "Patient")
    expect_equal(settings$interest_level, "D8")
    expect_equal(settings$reference_level, "D1")
    expect_equal(settings$assay_type, "tfmfeatures")
})

test_that("permutation test enumerates all swaps in a paired design", {
    sce <- simulate_sce()
    sce_single <- fit_loo(sce[, sce$Patient %in% c("P1", "P2", "P3")])
    expect_message(
        stats <- calculateStats(sce_single, target = "Drug", group = "Patient"),
        "all 8 relabelings"
    )

    expect_setequal(
        colnames(stats),
        c("Target", "Feature", "pvalue", "log2FoldChange")
    )
    expect_equal(attr(stats, "relabelings"), 8)
    expect_true(on_grid(stats$pvalue, 8))
    expect_gte(min(stats$pvalue), 1 / 8)
    expect_lte(max(stats$pvalue), 1)
})

test_that("permutation test samples swaps when there are too many", {
    sce <- simulate_sce()
    sce_single <- fit_loo(sce[, sce$Patient %in% c("P1", "P2", "P3")])
    stats <- suppressMessages(calculateStats(sce_single,
        target = "Drug", group = "Patient", n_permutations = 4
    ))

    expect_equal(attr(stats, "relabelings"), 5)
    expect_true(on_grid(stats$pvalue, 5))
})

test_that("permutation test permutes patients in an unpaired design", {
    sce <- simulate_sce()
    keep <- (sce$Patient %in% c("P1", "P2") & sce$Drug == "D1") |
        (sce$Patient %in% c("P3", "P4") & sce$Drug == "D8")
    sce_single <- fit_loo(sce[, keep])
    expect_message(
        stats <- calculateStats(sce_single, target = "Drug", group = "Patient"),
        "unpaired design"
    )

    expect_equal(attr(stats, "relabelings"), choose(4, 2))
    expect_true(on_grid(stats$pvalue, 6))
})

test_that("permutation test rejects mixed designs", {
    sce <- simulate_sce()
    keep <- sce$Patient %in% c("P1", "P2", "P3") &
        !(sce$Patient == "P1" & sce$Drug == "D8")
    sce_single <- fit_loo(sce[, keep])

    expect_error(
        calculateStats(sce_single, target = "Drug", group = "Patient"),
        "paired or an unpaired design"
    )
    stats <- calculateStats(sce_single,
        target = "Drug", group = "Patient", method = "t"
    )
    expect_equal(nrow(stats), 1)
})

test_that("permutation test needs the settings of predictLOO", {
    sce <- simulate_sce()
    sce_single <- fit_loo(sce[, sce$Patient %in% c("P1", "P2", "P3")])
    S4Vectors::metadata(sce_single)$predictLOO <- NULL

    expect_error(
        calculateStats(sce_single, target = "Drug", group = "Patient"),
        "predictLOO"
    )
    expect_no_error(calculateStats(sce_single,
        target = "Drug", group = "Patient", method = "t"
    ))
})

test_that("levels come from predictLOO, not from the factor order", {
    sce <- simulate_sce()
    sce <- sce[, sce$Patient %in% c("P1", "P2", "P3")]
    sce$Drug <- relevel(sce$Drug, ref = "D8")
    sce_single <- fit_loo(sce)
    stats <- calculateStats(sce_single,
        target = "Drug", group = "Patient", method = "t"
    )

    expect_equal(stats$Target, "D8")
    expect_gt(stats$log2FoldChange, 0)
})

test_that("volcano plot accepts precomputed statistics", {
    sce <- simulate_sce()
    sce_single <- fit_loo(sce[, sce$Patient %in% c("P1", "P2", "P3")])
    stats <- calculateStats(sce_single,
        target = "Drug", group = "Patient", method = "t"
    )
    p <- volcanoPlot(sce_single,
        target = "Drug", group = "Patient",
        stats = stats, p_cutoff = 0.05
    )

    expect_s3_class(p, "ggplot")
    expect_equal(rlang::as_name(p$mapping$x), "log2FoldChange")
})
