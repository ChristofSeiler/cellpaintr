#' Simulate CellProfiler data and write to a temporary file
#'
#' @importFrom readr read_csv write_csv
#' @importFrom stringr str_detect
#' @export
#'
#' @return path to csv file
#'
#' @examples
#' set.seed(23)
#' generate_data()
#'
generate_data <- function() {
    # use column names from CellProfiler and simulate data
    header_file <- system.file(
        "extdata", "header.csv",
        package = "cellpaintr", mustWork = TRUE
    )
    df <- read_csv(header_file, show_col_types = FALSE)

    # create a 384-well plate
    row <- LETTERS[seq(16)]
    col <- sprintf("%02d", seq(24))
    grid <- expand.grid(col = col, row = row)
    wells <- paste0(grid$row, grid$col)

    # set sample sizes
    n_patients <- 6
    n_images <- 100
    n_cells <- 1000
    n_drugs <- 8

    # only keep p features for each type
    p <- 10

    # parameters of distributions
    pois_lambda <- 10
    gamma_shape <- 1
    norm_mean <- 1

    # effect size of drug perturbation
    d7_shift <- 0.4
    d8_shift <- 0.8

    # simulate meta data
    ObjectNumber <- seq(n_cells)
    ImageNumber <- sample(n_images, n_cells, replace = TRUE)
    Metadata_Well <- sample(wells, n_cells, replace = TRUE)
    Metadata_Patient <- sample(paste0("P", seq(n_patients)), n_cells,
        replace = TRUE
    )
    Metadata_Drug <- sample(paste0("D", seq(n_drugs)), n_cells, replace = TRUE)

    # simulate non-negative shape features
    shape <- names(df)[str_detect(names(df), "AreaShape_")]
    shape <- sample(shape, p)
    shape_mat <-
        rpois(n = n_cells * length(shape), lambda = pois_lambda) |>
        matrix(nrow = n_cells, ncol = length(shape)) |>
        as.data.frame()
    colnames(shape_mat) <- shape

    # simulate non-negative intensity features
    intensity <- names(df)[str_detect(names(df), "Intensity_")]
    intensity <- sample(intensity, p)
    intensity_mat <-
        rgamma(n = n_cells * length(intensity), shape = gamma_shape) |>
        matrix(nrow = n_cells, ncol = length(intensity)) |>
        as.data.frame()
    colnames(intensity_mat) <- intensity

    # simulate real-valued texture features
    texture <- names(df)[str_detect(names(df), "Texture_")]
    texture <- sample(texture, p)
    texture_mat <-
        rnorm(n = n_cells * length(texture), mean = norm_mean) |>
        matrix(nrow = n_cells, ncol = length(texture)) |>
        as.data.frame()
    colnames(texture_mat) <- texture

    # combine into one data frame
    df <- data.frame(
        ObjectNumber, ImageNumber,
        Metadata_Well, Metadata_Patient, Metadata_Drug,
        shape_mat, intensity_mat, texture_mat
    )

    # shift texture features of D7
    cells_d07 <- which(df$Metadata_Drug == "D7")
    df[cells_d07, texture] <- df[cells_d07, texture] + d7_shift

    # shift texture features of D8
    cells_d08 <- which(df$Metadata_Drug == "D8")
    df[cells_d08, texture] <- df[cells_d08, texture] + d8_shift

    # write the simulated data frame to a temporary file
    cell_file <- tempfile(fileext = ".csv")
    write_csv(df, cell_file)
    cell_file
}

#' Load cell painting data from file and convert to a SingleCellExperiment
#'
#' @importFrom SingleCellExperiment SingleCellExperiment
#' @importFrom readr read_csv
#' @importFrom stringr str_detect str_remove
#' @export
#'
#' @param cell_file path to csv file from CellProfiler, e.g., MyExpt_Cells.csv
#' @return \code{\link[SingleCellExperiment]{SingleCellExperiment}} object
#'
#' @examples
#' set.seed(23)
#' cell_file <- generate_data()
#' loadData(cell_file)
#'
loadData <- function(cell_file) {
    # read file with cell data
    data <- read_csv(cell_file, show_col_types = FALSE)

    # remove features if they exist
    rmv_vars <- c("Number_Object_Number", "Parent_FilteredCell")
    keep_vars <- names(data)[!names(data) %in% rmv_vars]
    data <- data[, keep_vars]

    # split into meta data and features
    meta_cols <- str_detect(names(data), "Metadata_")
    meta_vars <- c(
        "ImageNumber", "ObjectNumber", names(data)[meta_cols],
        "Location_Center_X", "Location_Center_Y", "Location_Center_Z"
    )
    feature_vars <- names(data)[!names(data) %in% meta_vars]
    meta_cell <- data[, intersect(names(data), meta_vars)]
    features <- data[, feature_vars] |> t() # columns represent cells

    # add NAs to location if it doesn't exist
    if (sum(names(meta_cell) == "Location_Center_X") == 0) {
        meta_cell$Location_Center_X <- NA
    }
    if (sum(names(meta_cell) == "Location_Center_Y") == 0) {
        meta_cell$Location_Center_Y <- NA
    }
    if (sum(names(meta_cell) == "Location_Center_Z") == 0) {
        meta_cell$Location_Center_Z <- NA
    }

    # remove meta label from names
    names(meta_cell) <- str_remove(names(meta_cell), "Metadata_")

    # store meta and features in bioconductor object
    sce <- SingleCellExperiment(
        list(features = features),
        colData = meta_cell
    )
    sce
}

#' Plot number of cells per image
#'
#' @importFrom SingleCellExperiment colData
#' @importFrom ggplot2 ggplot geom_histogram scale_x_log10 ggtitle
#' @importFrom dplyr group_by tally
#' @importFrom tibble as_tibble
#' @export
#'
#' @param sce \code{\link[SingleCellExperiment]{SingleCellExperiment}} object
#' @param bins Number of histogram bins
#' @return \code{\link[ggplot2]{ggplot2}} object
#'
#' @examples
#' set.seed(23)
#' cell_file <- generate_data()
#' sce <- loadData(cell_file)
#' plotCellsPerImage(sce)
#'
plotCellsPerImage <- function(sce, bins = 100) {
    colData(sce) |>
        as_tibble() |>
        group_by(ImageNumber) |>
        tally() |>
        ggplot(aes(n)) +
        geom_histogram(bins = bins) +
        scale_x_log10() +
        ggtitle("Number of Cells Per Image")
}

#' Plot number of cells per image
#'
#' @importFrom SummarizedExperiment assay
#' @importFrom SingleCellExperiment reducedDim
#' @importFrom ggplot2 ggplot geom_tile aes scale_fill_gradient2 theme_minimal
#' @importFrom ggplot2 coord_fixed ylab
#' @importFrom dplyr mutate top_n
#' @importFrom tibble as_tibble
#' @importFrom tidyr pivot_longer
#' @export
#'
#' @param sce \code{\link[SingleCellExperiment]{SingleCellExperiment}} object
#' @param filter_by PC to use for feature selection
#' @param top Number of top features to select
#' @param pcs Number of PCs to plot
#' @return \code{\link[ggplot2]{ggplot2}} object
#'
#' @examples
#' set.seed(23)
#' cell_file <- generate_data()
#' sce <- loadData(cell_file)
#' sce <- transformLogScale(sce, robust = TRUE)
#' sce <- scater::runPCA(sce, exprs_values = "tfmfeatures", ncomponents = 10)
#' plotPCACor(sce, filter_by = 1)
#'
plotPCACor <- function(sce, filter_by = 1, top = 20, pcs = seq(5)) {
    features <- assay(sce, "tfmfeatures")
    scores <- reducedDim(sce, "PCA")[, pcs]
    features_pcs <- cor(t(as.matrix(features)), scores)

    keep <- rownames(features_pcs)
    importance <- abs(features_pcs[, filter_by])

    features_pcs <- features_pcs |>
        as_tibble() |>
        mutate(name = keep, importance = importance) |>
        top_n(top, importance)

    features_pcs <- features_pcs |>
        pivot_longer(!c(name, importance),
            names_to = "PC", values_to = "cor"
        ) |>
        mutate(PC = str_remove(PC, "PC"))

    ggplot(features_pcs, aes(PC, name, fill = cor)) +
        geom_tile(color = "white") +
        scale_fill_gradient2(
            low = "blue", high = "red", mid = "white",
            midpoint = 0, limit = c(-1, 1), space = "Lab",
            name = "Pearson\nCorrelation"
        ) +
        theme_minimal() +
        coord_fixed() +
        ylab("feature name")
}

#' Remove cells if not enough or too many in one image
#'
#' @importFrom SingleCellExperiment colData
#' @importFrom dplyr group_by tally filter pull
#' @importFrom tibble as_tibble
#' @export
#'
#' @param sce \code{\link[SingleCellExperiment]{SingleCellExperiment}} object
#' @param min Remove cells below that number
#' @param max Remove cells above that number
#' @return \code{\link[SingleCellExperiment]{SingleCellExperiment}} object
#'
#' @examples
#' set.seed(23)
#' cell_file <- generate_data()
#' sce <- loadData(cell_file)
#' sce <- removeOutliers(sce, min = 0, max = 300)
#'
removeOutliers <- function(sce, min, max) {
    cell_ids <-
        colData(sce) |>
        as_tibble() |>
        group_by(ImageNumber) |>
        tally() |>
        filter(n >= min & n <= max) |>
        pull(ImageNumber)
    sce[, sce$ImageNumber %in% cell_ids]
}

#' Remove cells with missing features
#'
#' @importFrom SummarizedExperiment assay
#' @export
#'
#' @param sce \code{\link[SingleCellExperiment]{SingleCellExperiment}} object
#' @return \code{\link[SingleCellExperiment]{SingleCellExperiment}} object
#'
#' @examples
#' set.seed(23)
#' cell_file <- generate_data()
#' sce <- loadData(cell_file)
#' sce <- removeNAs(sce)
#'
removeNAs <- function(sce) {
    mat <- assay(sce, "features")

    miss_cells <- apply(mat, 2, function(x) sum(is.na(x)))
    cell_ids <- which(miss_cells > 0)
    if (length(cell_ids) > 0) {
        sce <- sce[, -cell_ids]
    }
    sce
}

#' Filter low variance features
#'
#' @importFrom SummarizedExperiment assay
#' @export
#'
#' @param sce \code{\link[SingleCellExperiment]{SingleCellExperiment}} object
#' @param threshold Keep features that have larger variance than this threshold
#' @param robust If true use median absolute deviation, otherwise use variance
#' @return \code{\link[SingleCellExperiment]{SingleCellExperiment}} object
#'
#' @examples
#' set.seed(23)
#' cell_file <- generate_data()
#' sce <- loadData(cell_file)
#' sce <- removeLowVariance(sce)
#'
removeLowVariance <- function(sce, threshold = 0, robust = FALSE) {
    mat <- assay(sce, "features")

    if (!robust) {
        spread <- apply(mat, 1, var)
    } else {
        spread <- apply(mat, 1, mad)
    }
    feature_ids <- which(spread > threshold)
    sce[feature_ids, ]
}

#' Remove zero-inflated features
#'
#' @importFrom SummarizedExperiment assay
#' @export
#'
#' @param sce \code{\link[SingleCellExperiment]{SingleCellExperiment}} object
#' @param proportion Remove features exceeding this zero-inflation proportion
#' @return \code{\link[SingleCellExperiment]{SingleCellExperiment}} object
#'
#' @examples
#' set.seed(23)
#' cell_file <- generate_data()
#' sce <- loadData(cell_file)
#' sce <- removeZeroInflation(sce)
#'
removeZeroInflation <- function(sce, proportion = 0.2) {
    mat <- assay(sce, "features")
    prop_zeros <- apply(mat, 1, function(x) mean(x == 0))
    sce[!(prop_zeros > proportion), ]
}

#' Filter low variance features
#'
#' @importFrom SummarizedExperiment assay assay<-
#' @export
#'
#' @param sce \code{\link[SingleCellExperiment]{SingleCellExperiment}} object
#' @param robust If true robust z-score, otherwise standard z-score
#' @return \code{\link[SingleCellExperiment]{SingleCellExperiment}} object
#'
#' @examples
#' set.seed(23)
#' cell_file <- generate_data()
#' sce <- loadData(cell_file)
#' sce <- transformLogScale(sce)
#'
transformLogScale <- function(sce, robust = FALSE) {
    mat <- assay(sce, "features")

    # log(x+1) transform on non-negative valued features
    non_neg_features <- apply(mat, 1, function(x) sum(x >= 0) == length(x))
    mat[non_neg_features, ] <- log1p(mat[non_neg_features, ])

    if (!robust) {
        # option: standard z-score
        mat <- mat |>
            t() |>
            scale() |>
            t()
    } else {
        # option: robust z-score
        median_features <- apply(mat, 1, median)
        mad_features <- apply(mat, 1, mad)
        mat <- (mat - median_features) / mad_features
    }

    # add transformed features
    assay(sce, "tfmfeatures") <- mat
    sce
}

#' Internal function to compute y_hat on a subset of the features
#'
#' @importFrom SummarizedExperiment assay
#' @importFrom SingleCellExperiment colData reducedDim reducedDim<-
#' @importFrom dplyr bind_cols bind_rows arrange select
#' @importFrom parsnip rand_forest set_mode set_engine
#' @importFrom rsample group_initial_split training testing
#' @importFrom purrr map
#' @importFrom ranger ranger
#' @importFrom stringr str_starts
#'
#' @param feature_name String with feature name
#' @param sce_feature \code{\link[SingleCellExperiment]{SingleCellExperiment}}
#' @param starts Starting string
#' @param assay_type A string specifying the assay
#' @param target Name of target variable for prediction
#' @param interest_level Factor interest level in `target` variable
#' @param reference_level Factor reference level in `target` variable
#' @param group Grouping variable for cross-validation, e.g., patient
#' @param weights Weights variable when features are aggregated
#' @param n_threads Number of parallel threads for fitting of models
#' @return \code{\link[tibble]{tibble}} data frame
#'
compute_y_hat <- function(feature_name,
                          sce_feature,
                          starts,
                          assay_type,
                          target,
                          interest_level,
                          reference_level,
                          group,
                          weights,
                          n_threads) {
    # convert to regular expression
    if (feature_name == "all") {
        pattern <- ".*"
    } else {
        pattern <- feature_name
    }

    # subset features
    features <- names(sce_feature)
    if (starts) {
        features <- features[str_starts(features, pattern)]
    } else {
        features <- features[str_detect(features, pattern)]
    }
    sce_feature <- sce_feature[features, ]

    # prepare data frame with target variable and predictor matrix
    X <- assay(sce_feature, assay_type) |> t()
    cells <- colData(sce_feature) |> as_tibble()
    ml_data <- data.frame(Target = cells |> pull(target), X)
    ml_data$Target <- factor(ml_data$Target,
        levels = c(interest_level, reference_level)
    )

    # compute y hats
    patients <- unique(sce_feature[[group]])
    result_list <- lapply(patients, function(patient) {
        train_ids <- which(sce_feature[[group]] != patient)
        test_ids <- which(sce_feature[[group]] == patient)

        # train model
        if (!is.null(weights)) {
            case_weights <- sce_feature[[weights]][train_ids]
        } else {
            case_weights <- NULL
        }
        rf_spec <-
            parsnip::rand_forest() |>
            parsnip::set_mode("classification") |>
            parsnip::set_engine("ranger",
                num.threads = n_threads,
                case.weights = case_weights
            )
        rf_fit <- rf_spec |>
            parsnip::fit(Target ~ ., data = ml_data[train_ids, ])

        # evaluate model
        dplyr::bind_cols(
            cell = test_ids,
            predict(rf_fit, ml_data[test_ids, ], type = "prob")
        )
    })

    # combine folds
    result <- result_list |>
        bind_rows() |>
        arrange(cell)

    # rename to feature type / channel
    names(result)[2] <- feature_name
    result |> select(all_of(feature_name))
}

#' Predict target from features
#'
#' @importFrom SummarizedExperiment assay
#' @importFrom SingleCellExperiment colData reducedDim reducedDim<-
#' @importFrom dplyr bind_cols bind_rows arrange select
#' @importFrom parsnip rand_forest set_mode set_engine
#' @importFrom rsample group_initial_split training testing
#' @importFrom purrr map
#' @importFrom ranger ranger
#' @importFrom stringr str_starts
#' @export
#'
#' @param sce \code{\link[SingleCellExperiment]{SingleCellExperiment}} object
#' @param assay_type A string specifying the assay
#' @param target Name of target variable for prediction
#' @param interest_level Factor interest level in `target` variable
#' @param reference_level Factor reference level in `target` variable
#' @param types Vector of strings of feature types
#' @param channels Vector of strings of feature channels
#' @param group Grouping variable for cross-validation, e.g., patient
#' @param weights Weights variable when features are aggregated
#' @param n_threads Number of parallel threads for fitting of models
#' @return \code{\link[tibble]{tibble}} data frame
#'
#' @examples
#' set.seed(23)
#' cell_file <- generate_data()
#' sce <- loadData(cell_file)
#' sce <- transformLogScale(sce)
#'
#' sce$Drug <- as.factor(sce$Drug)
#' sce$Drug <- relevel(sce$Drug, ref = "D1")
#' types <- c("AreaShape", "Intensity", "Texture")
#'
#' sce_single <- predictLOO(
#'     sce,
#'     target = "Drug", group = "Patient",
#'     interest_level = "D7", reference_level = "D1",
#'     types = types,
#'     n_threads = 1
#' )
#'
predictLOO <- function(sce,
                       assay_type = "tfmfeatures",
                       target = "Treatment",
                       interest_level,
                       reference_level,
                       types = NULL,
                       channels = NULL,
                       group = "Patient",
                       weights = NULL,
                       n_threads = 1) {
    # subset for binary classification
    sce_subset <- sce[, sce[[target]] %in% c(reference_level, interest_level)]

    # remove unused levels
    sce_subset[[target]] <- droplevels(sce_subset[[target]])

    # combine and add to reducedDim slot
    y_hat <- compute_y_hat("all", sce_subset,
        starts = TRUE, assay_type,
        target, interest_level, reference_level, group,
        weights, n_threads
    )
    if (length(types) > 0) {
        y_hat <- bind_cols(
            y_hat,
            purrr::map(types, compute_y_hat, sce_subset,
                starts = TRUE, assay_type, target, interest_level,
                reference_level, group, weights, n_threads,
                .progress = TRUE
            ) |>
                bind_cols()
        )
    }
    if (length(channels) > 0) {
        y_hat <- bind_cols(
            y_hat,
            purrr::map(channels, compute_y_hat, sce_subset,
                starts = FALSE, assay_type, target, interest_level,
                reference_level, group, weights, n_threads,
                .progress = TRUE
            ) |>
                bind_cols()
        )
    }
    y_hat <- y_hat |> as.matrix()

    # add y_hat to sce object
    cell_id <- seq(ncol(sce_subset))
    rownames(y_hat) <- cell_id
    colnames(sce_subset) <- cell_id
    reducedDim(sce_subset, "prevalidated") <- y_hat
    sce_subset
}

#' Aggregate predicted leave-one-out probabilities over meta variables
#'
#' @importFrom scater aggregateAcrossCells
#'
#' @param sce \code{\link[SingleCellExperiment]{SingleCellExperiment}} object
#' @param assay_type A string specifying the assay
#' @param meta_vars a vector of variables from `colData`
#' @return \code{data.frame}
#'
aggregateYhat <- function(sce,
                          assay_type = "tfmfeatures",
                          meta_vars = c("Patient", "Treatment")) {
    summed <- aggregateAcrossCells(
        sce,
        id = colData(sce)[, meta_vars],
        use.assay.type = assay_type, statistics = "mean"
    )

    cbind(
        colData(summed)[, meta_vars] |> as.data.frame(),
        reducedDim(summed, type = "prevalidated")
    )
}

#' Plot predicted leave-one-out probabilities
#'
#' @importFrom ggplot2 ggplot aes geom_boxplot geom_jitter ylab facet_wrap
#' @importFrom ggplot2 geom_hline
#' @importFrom tidyr pivot_longer
#' @export
#'
#' @param sce \code{\link[SingleCellExperiment]{SingleCellExperiment}} object
#' @param assay_type A string specifying the assay
#' @param meta_vars a vector of variables from `colData`
#' @param target Name of target variable for prediction
#' @return \code{\link[ggplot2]{ggplot2}} object
#'
#' @examples
#' set.seed(23)
#' cell_file <- generate_data()
#' sce <- loadData(cell_file)
#' sce <- transformLogScale(sce)
#'
#' sce$Drug <- as.factor(sce$Drug)
#' sce$Drug <- relevel(sce$Drug, ref = "D1")
#' types <- c("AreaShape", "Intensity", "Texture")
#'
#' sce_single <- predictLOO(
#'     sce,
#'     target = "Drug", group = "Patient",
#'     interest_level = "D7", reference_level = "D1",
#'     types = types,
#'     n_threads = 1
#' )
#'
#' plotLOO(sce_single, meta_vars = c("Patient", "Drug"), target = "Drug")
#'
plotLOO <- function(sce,
                    assay_type = "tfmfeatures",
                    meta_vars = c("Patient", "Treatment"),
                    target = "Treatment") {
    y_hat <- aggregateYhat(sce, assay_type, meta_vars)

    y_hat |>
        pivot_longer(
            cols = -c(all_of(meta_vars)),
            names_to = "features", values_to = "value"
        ) |>
        ggplot(aes(.data[[target]], value, color = .data[[target]])) +
        geom_boxplot(width = 0.2, outliers = FALSE) +
        geom_jitter(width = 0.1) +
        ylab("predicted leave-one-out probability") +
        facet_wrap(~features)
}

#' Aggregate predicted leave-one-out probabilities over meta variables
#' over a list of SingleCellExperiment objects
#'
#' @importFrom dplyr pull bind_rows all_of
#' @importFrom scater aggregateAcrossCells
#' @importFrom tidyr pivot_wider
#' @importFrom stats t.test
#' @export
#'
#' @param sce A \code{\link[SingleCellExperiment]{SingleCellExperiment}} object
#' @param assay_type A string specifying the assay
#' @param meta_vars a vector of variables from `colData`
#' @param target Name of target variable for prediction
#' @return \code{data.frame}
#'
#' @examples
#' set.seed(23)
#' cell_file <- generate_data()
#' sce <- loadData(cell_file)
#' sce <- transformLogScale(sce)
#'
#' sce$Drug <- as.factor(sce$Drug)
#' sce$Drug <- relevel(sce$Drug, ref = "D1")
#' types <- c("AreaShape", "Intensity", "Texture")
#'
#' sce_single <- predictLOO(
#'     sce,
#'     target = "Drug", group = "Patient",
#'     interest_level = "D7", reference_level = "D1",
#'     types = types,
#'     n_threads = 1
#' )
#'
#' calculateStats(sce_single, meta_vars = c("Patient", "Drug"), target = "Drug")
#'
calculateStats <- function(sce,
                           assay_type = "tfmfeatures",
                           meta_vars = c("Patient", "Treatment"),
                           target = "Treatment") {
    # standard convention in base R to store reference level as first
    reference_level <- levels(sce[[target]])[1]
    interest_level <- levels(sce[[target]])[2]

    y_hat <- aggregateYhat(sce, assay_type, meta_vars)
    feature_vars <- setdiff(names(y_hat), meta_vars)
    y_hat$Target <- interest_level

    # test one feature at a time
    lapply(feature_vars, function(feature) {
        wide <- y_hat |>
            select(all_of(c(meta_vars, feature))) |>
            pivot_wider(
                names_from = all_of(target),
                values_from = all_of(feature)
            )

        x <- pull(wide, all_of(reference_level))
        y <- pull(wide, all_of(interest_level))
        paired <- sum(is.na(c(x, y))) == 0
        result <- t.test(x, y,
            paired = paired, var.equal = TRUE,
            alternative = "less"
        )

        data.frame(
            Target = interest_level,
            Feature = feature,
            pvalue = result$p.value,
            log2FoldChange = log2(mean(y, na.rm = TRUE) / mean(x, na.rm = TRUE))
        )
    }) |> bind_rows()
}

#' Plot predicted leave-one-out probabilities
#'
#' @importFrom ggplot2 ggplot aes geom_vline geom_hline xlab ylab geom_point
#' @importFrom dplyr mutate
#' @importFrom ggrepel geom_text_repel
#' @export
#'
#' @param sce A \code{\link[SingleCellExperiment]{SingleCellExperiment}} object
#' @param assay_type A string specifying the assay
#' @param meta_vars a vector of variables from `colData`
#' @param target Name of target variable for prediction
#' @param p_cutoff Cut-off for statistical significance. A horizontal line
#'                 will be drawn at -log10(p_cutoff).
#' @param fc_cutoff Cut-off for absolute log2 fold-change. A vertical lines
#'                  will be drawn at fc_cutoff.
#' @return \code{\link[ggplot2]{ggplot2}} object
#'
#' @examples
#' set.seed(23)
#' cell_file <- generate_data()
#' sce <- loadData(cell_file)
#' sce <- transformLogScale(sce)
#'
#' sce$Drug <- as.factor(sce$Drug)
#' sce$Drug <- relevel(sce$Drug, ref = "D1")
#' types <- c("AreaShape", "Intensity", "Texture")
#'
#' sce_single <- predictLOO(
#'     sce,
#'     target = "Drug", group = "Patient",
#'     interest_level = "D7", reference_level = "D1",
#'     types = types,
#'     n_threads = 1
#' )
#'
#' volcanoPlot(sce_single,
#'     meta_vars = c("Patient", "Drug"), target = "Drug",
#'     p_cutoff = 0.05, fc_cutoff = 0.5
#' )
#'
volcanoPlot <- function(sce,
                        assay_type = "tfmfeatures",
                        meta_vars = c("Patient", "Treatment"),
                        target = "Treatment",
                        p_cutoff = NULL,
                        fc_cutoff = 1.0) {
    stats <- calculateStats(sce, assay_type, meta_vars, target)

    if (is.null(p_cutoff)) {
        p_cutoff <- 0.01 / nrow(stats)
    }

    stats |>
        mutate(Feature = ifelse(pvalue < p_cutoff & log2FoldChange > fc_cutoff,
            Feature, ""
        )) |>
        ggplot(aes(log2FoldChange, -log10(pvalue), label = Feature)) +
        geom_vline(xintercept = c(0, fc_cutoff), alpha = 0.5, linetype = "dashed") +
        geom_hline(
            yintercept = c(0, -log10(p_cutoff)), alpha = 0.5,
            linetype = "dashed"
        ) +
        xlab("log2 fold change") +
        ylab("-log10 p-value") +
        geom_point() +
        geom_text_repel()
}

#' Plot ROC curves
#'
#' @importFrom ggplot2 ggplot geom_path geom_abline coord_equal
#' @importFrom ggplot2 scale_x_continuous scale_y_continuous xlab ylab ggtitle
#' @importFrom ggplot2 facet_wrap
#' @importFrom dplyr group_by mutate pull
#' @importFrom yardstick roc_auc roc_curve
#' @importFrom SingleCellExperiment reducedDim
#' @importFrom tidyr pivot_longer
#' @export
#'
#' @param sce \code{\link[SingleCellExperiment]{SingleCellExperiment}} object
#' @param assay_type A string specifying the assay
#' @param meta_vars a vector of variables from `colData`
#' @param target Name of target variable for prediction
#' @return \code{\link[ggplot2]{ggplot2}} object
#'
#' @examples
#' set.seed(23)
#' cell_file <- generate_data()
#' sce <- loadData(cell_file)
#' sce <- transformLogScale(sce)
#'
#' sce$Drug <- as.factor(sce$Drug)
#' sce$Drug <- relevel(sce$Drug, ref = "D1")
#' types <- c("AreaShape", "Intensity", "Texture")
#'
#' sce_single <- predictLOO(
#'     sce,
#'     target = "Drug", group = "Patient",
#'     interest_level = "D7", reference_level = "D1",
#'     types = types,
#'     n_threads = 1
#' )
#'
#' plotROC(sce_single, meta_vars = c("Patient", "Drug"), target = "Drug")
#'
plotROC <- function(sce, assay_type = "tfmfeatures",
                    meta_vars = c("Patient", "Treatment"),
                    target = "Treatment") {
    summed <- aggregateAcrossCells(
        sce,
        id = colData(sce)[, meta_vars],
        use.assay.type = assay_type, statistics = "mean"
    )

    result <- cbind(
        colData(summed)[, meta_vars] |> as.data.frame(),
        reducedDim(summed, type = "prevalidated")
    ) |>
        pivot_longer(
            cols = -c(all_of(meta_vars)),
            names_to = "features", values_to = "pred"
        ) |>
        droplevels()

    # title
    reference_level <- levels(result[[target]])[1]
    interest_level <- levels(result[[target]])[2]
    mean_auc <- result |>
        group_by(features) |>
        yardstick::roc_auc(all_of(target), "pred", event_level = "second") |>
        mutate(.estimate = round(.estimate, digits = 3)) |>
        pull(.estimate) |>
        mean() |>
        round(digits = 3)
    title_str <- paste0(
        "Predict ", interest_level, " vs ", reference_level,
        " (mean AUC = ", mean_auc, ")"
    )

    result |>
        group_by(features) |>
        yardstick::roc_curve(all_of(target), "pred", event_level = "second") |>
        ggplot(aes(x = 1 - specificity, y = sensitivity)) +
        geom_path() +
        geom_abline(lty = 3) +
        coord_equal() +
        scale_x_continuous(breaks = c(0, 0.5, 1)) +
        scale_y_continuous(breaks = c(0, 0.5, 1)) +
        xlab("false positive rate (1 - specificity)") +
        ylab("true positive rate (sensitivity)") +
        ggtitle(title_str) +
        facet_wrap(~features)
}

#' Plot AUC comparison
#'
#' @importFrom SingleCellExperiment colData
#' @importFrom ggplot2 ggplot xlab ggtitle geom_jitter aes
#' @importFrom dplyr group_by bind_rows summarize arrange pull
#' @importFrom yardstick roc_auc
#' @importFrom stringr str_remove
#' @importFrom tidyr pivot_longer
#' @importFrom scater aggregateAcrossCells
#' @export
#'
#' @param sce A \code{\link[SingleCellExperiment]{SingleCellExperiment}} object
#' @param assay_type A string specifying the assay
#' @param meta_vars a vector of variables from `colData`
#' @param target Name of target variable for prediction
#' @return \code{\link[ggplot2]{ggplot2}} object
#'
#' @examples
#' set.seed(23)
#' cell_file <- generate_data()
#' sce <- loadData(cell_file)
#' sce <- transformLogScale(sce)
#'
#' sce$Drug <- as.factor(sce$Drug)
#' sce$Drug <- relevel(sce$Drug, ref = "D1")
#' types <- c("AreaShape", "Intensity", "Texture")
#'
#' sce_single <- predictLOO(
#'     sce,
#'     target = "Drug", group = "Patient",
#'     interest_level = "D7", reference_level = "D1",
#'     types = types,
#'     n_threads = 1
#' )
#'
#' plotAUC(sce_single, meta_vars = c("Patient", "Drug"), target = "Drug")
#'
plotAUC <- function(sce,
                    assay_type = "tfmfeatures",
                    meta_vars = c("Patient", "Treatment"),
                    target = "Treatment") {
    summed <- aggregateAcrossCells(
        sce,
        id = colData(sce)[, meta_vars],
        use.assay.type = assay_type, statistics = "mean"
    )

    result <- cbind(
        colData(summed)[, meta_vars] |> as.data.frame(),
        reducedDim(summed, type = "prevalidated")
    ) |>
        pivot_longer(
            cols = -c(all_of(meta_vars)),
            names_to = "features", values_to = "pred"
        ) |>
        droplevels()

    reference_level <- levels(result[[target]])[1]
    interest_level <- levels(result[[target]])[2]
    title_str <- paste0("Predict ", interest_level, " vs ", reference_level)

    aucs <- result |>
        group_by(features) |>
        yardstick::roc_auc(all_of(target), "pred", event_level = "second")

    fct_order <- aucs |>
        group_by(features) |>
        summarize(median = median(.estimate)) |>
        arrange(desc(median)) |>
        pull(features)

    aucs |>
        mutate(features = factor(features, levels = fct_order)) |>
        ggplot(aes(.estimate, features)) +
        geom_point() +
        xlab("AUC") +
        ggtitle(title_str)
}
