simulate_paired <- function() {
    set.seed(23)
    cell_file <- generate_data()
    sce <- loadData(cell_file)
    sce <- transformScale(sce)
    sce$Drug <- relevel(as.factor(sce$Drug), ref = "D1")
    sce
}

test_that("corrected t-test without cross term is the paired t-test", {
    d <- c(0.02, 0.05, -0.01, 0.04, 0.03)
    result <- corrected_t(d, matrix(0, 5, 5))
    expect_equal(result$pvalue, t.test(d, alternative = "greater")$p.value)
    expect_equal(result$variance_ratio, 1)
})

test_that("corrected t-test adds the cross term to the sample variance", {
    d <- c(0.02, 0.05, -0.01, 0.04)
    a <- matrix(c(
        0, 0.01, 0.02, 0.00,
        0.03, 0, 0.01, 0.02,
        0.01, 0.02, 0, 0.01,
        0.02, 0.01, 0.03, 0
    ), nrow = 4, byrow = TRUE)
    cross <- sum((a * t(a))[row(a) != col(a)])
    variance <- var(d) + cross / 3
    result <- corrected_t(d, a)
    expect_equal(result$pvalue, pt(mean(d) / sqrt(variance / 4), 3, lower.tail = FALSE))
    expect_equal(result$variance_ratio, variance / var(d))
})

test_that("a negative cross term leaves the sample variance unchanged", {
    d <- c(0.02, 0.05, -0.01)
    a <- matrix(c(0, 0.02, 0.01, -0.02, 0, 0.03, -0.01, -0.03, 0), nrow = 3)
    result <- corrected_t(d, a)
    expect_equal(result$pvalue, t.test(d, alternative = "greater")$p.value)
    expect_equal(result$variance_ratio, 1)
})

test_that("corrected t-test runs on a paired design", {
    sce <- simulate_paired()
    set.seed(1)
    sce_single <- predictLOO(sce[, sce$Patient %in% c("P1", "P2", "P3", "P4")],
        target = "Drug", group = "Patient",
        interest_level = "D8", reference_level = "D1",
        types = "Texture", n_threads = 1, progress = FALSE
    )
    expect_message(
        stats <- calculateStats(sce_single,
            target = "Drug", group = "Patient", method = "corrected"
        ),
        "6 leave-two-out fits per feature group"
    )
    t_stats <- calculateStats(sce_single,
        target = "Drug", group = "Patient", method = "t"
    )
    expect_equal(stats$Feature, c("all", "Texture"))
    expect_equal(stats$log2FoldChange, t_stats$log2FoldChange)
    expect_true(all(stats$pvalue >= 0 & stats$pvalue <= 1))
    ratio <- attr(stats, "variance_ratio")
    expect_named(ratio, stats$Feature)
    expect_true(all(ratio >= 1))
})

test_that("corrected t-test needs a paired design", {
    sce <- simulate_paired()
    keep <- (sce$Patient %in% c("P1", "P2", "P3") & sce$Drug == "D1") |
        (sce$Patient %in% c("P4", "P5", "P6") & sce$Drug == "D8")
    set.seed(1)
    sce_single <- predictLOO(sce[, keep],
        target = "Drug", group = "Patient",
        interest_level = "D8", reference_level = "D1",
        n_threads = 1, progress = FALSE
    )
    expect_error(
        calculateStats(sce_single,
            target = "Drug", group = "Patient", method = "corrected"
        ),
        "paired design"
    )
})

test_that("corrected t-test needs the settings of predictLOO()", {
    sce <- simulate_paired()
    set.seed(1)
    sce_single <- predictLOO(sce[, sce$Patient %in% c("P1", "P2", "P3")],
        target = "Drug", group = "Patient",
        interest_level = "D8", reference_level = "D1",
        n_threads = 1, progress = FALSE
    )
    S4Vectors::metadata(sce_single)$predictLOO <- NULL
    expect_error(
        calculateStats(sce_single,
            target = "Drug", group = "Patient", method = "corrected"
        ),
        "refits the models"
    )
})
