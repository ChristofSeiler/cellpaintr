#' Simulate CellProfiler data and write to a temporary file
#'
#' @importFrom readr read_csv write_csv
#' @importFrom stringr str_detect
#' @importFrom stats rgamma rnorm rpois
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

    # effect size of drug perturbation
    d7_shift <- 0.2
    d8_shift <- 0.4

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
        rgamma(n = n_cells * length(texture), shape = gamma_shape) |>
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
#' @param assay_type A string specifying the assay
#' @return \code{\link[ggplot2]{ggplot2}} object
#'
#' @examples
#' set.seed(23)
#' cell_file <- generate_data()
#' sce <- loadData(cell_file)
#' sce <- transformScale(sce, robust = TRUE)
#' sce <- scater::runPCA(sce, exprs_values = "tfmfeatures", ncomponents = 10)
#' plotPCACor(sce, filter_by = 1)
#'
plotPCACor <- function(sce, filter_by = 1, top = 20, pcs = seq(5),
                       assay_type = "tfmfeatures") {
    features <- assay(sce, assay_type)
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

#' Transform, center, and scale features
#'
#' @importFrom SummarizedExperiment assay assay<-
#' @importFrom stats qlogis qnorm
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
#' sce <- transformScale(sce)
#'
transformScale <- function(sce, robust = FALSE) {
    mat <- assay(sce, "features")

    # variance features
    var_ids <- str_detect(
        rownames(mat),
        paste(
            c(
                "Correlation_K_",
                "Variance_",
                "Contrast",
                "MeanIntensity",
                "StdIntensity",
                "IntegratedIntensity",
                "SumAverage",
                "MedianIntensity",
                "MADIntensity",
                "UpperQuartileIntensity",
                "LowerQuartileIntensity",
                "MaxIntensity",
                "MinIntensity"
            ),
            collapse = "|"
        )
    )
    submat <- mat[var_ids, ]
    submat <- pmax(submat, 0)
    submat <- log(submat)
    mat[var_ids, ] <- submat

    # correlation features
    corr_ids <- str_detect(rownames(mat), "_Correlation_")
    submat <- mat[corr_ids, ]
    submat <- pmax(pmin(submat, 1), -1)
    submat <- atanh(submat)
    mat[corr_ids, ] <- submat

    # proportion features
    prop_ids <- str_detect(rownames(mat), "Texture_AngularSecondMoment_")
    submat <- mat[prop_ids, ]
    submat <- pmax(pmin(submat, 1), 0)
    submat <- qlogis(submat)
    mat[prop_ids, ] <- submat

    # extreme value features
    prop_ids <- str_detect(
        rownames(mat),
        "Correlation_Costes_|Correlation_Manders_"
    )
    submat <- mat[prop_ids, ]
    quantile_norm <- function(x) {
        n <- sum(!is.na(x))
        r <- rank(x, na.last = "keep", ties.method = "random")
        # Blom adjustment to avoid -Inf and Inf
        qnorm((r - 0.375) / (n + 0.25))
    }
    submat <- t(apply(submat, 1, quantile_norm))
    mat[prop_ids, ] <- submat

    # fix infinities
    cap_inf <- function(x) {
        finite_x <- x[is.finite(x)]
        max_val <- max(finite_x, na.rm = TRUE)
        min_val <- min(finite_x, na.rm = TRUE)
        x[x == Inf] <- max_val
        x[x == -Inf] <- min_val
        x
    }
    mat <- t(apply(mat, 1, cap_inf))

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

#' Remove row and column effect using the median polish procedure
#'
#' @importFrom SummarizedExperiment assay assay<-
#' @importFrom stats qlogis qnorm medpolish
#' @importFrom scrapper aggregateAcrossCells
#' @importFrom tidyr pivot_wider
#' @importFrom cli cli_abort
#' @importFrom dplyr left_join select join_by
#' @export
#'
#' @param sce \code{\link[SingleCellExperiment]{SingleCellExperiment}} object
#' @param assay_type A string specifying the assay
#' @return \code{\link[SingleCellExperiment]{SingleCellExperiment}} object
#'
#' @examples
#' set.seed(23)
#' cell_file <- generate_data()
#' sce <- loadData(cell_file)
#' sce <- transformScale(sce)
#' row <- substr(sce$Well, 1, 1)
#' row <- factor(row, levels = LETTERS[1:length(unique(row))])
#' col <- substr(sce$Well, 2, 3)
#' col <- factor(as.numeric(col))
#' sce$Row <- row
#' sce$Col <- col
#' sce <- removeRowColEffect(sce)
#'
removeRowColEffect <- function(sce,
                               assay_type = "tfmfeatures") {
    # check input
    if (!is.factor(sce$Col)) {
        cli::cli_abort(c(
            "'Col' in {.code colData(sce)$Col} must be a factor",
            "x" = "You've supplied a {.cls {class(sce$Col)}}."
        ))
    }
    if (!is.factor(sce$Row)) {
        cli::cli_abort(c(
            "'Row' in {.code colData(sce)$Row} must be a factor",
            "x" = "You've supplied a {.cls {class(sce$Row)}}."
        ))
    }

    # retrieve count matrix
    mat <- assay(sce, assay_type)

    # get group index
    aggr <- scrapper::aggregateAcrossCells(
        assay(sce, assay_type),
        factors = colData(sce)[, c("Col", "Row")]
    )

    # calculate median per well
    idx <- aggr$index
    mat_well <- assay(sce, assay_type)
    mat_well <- t(apply(mat_well, 1, function(x) tapply(x, idx, median)))

    # loop over all features
    corrected <- mat
    for (feature_ind in seq_len(nrow(sce))) {
        # convert to matrix
        df_plate_layout_well <- data.frame(
            feature = mat_well[feature_ind, ],
            aggr$combinations
        )
        row_col_well_mat <- df_plate_layout_well |>
            pivot_wider(names_from = Col, values_from = feature) |>
            select(-Row) |>
            as.matrix()
        # NAs occur if no cells are in a well: setting NAs to zero
        row_col_well_mat[is.na(row_col_well_mat)] <- 0
        med <- medpolish(row_col_well_mat, trace.iter = FALSE)

        # subtract medians from full sce object
        df_plate_layout_full <- data.frame(
            feature = mat[feature_ind, ],
            colData(sce)[, c("Col", "Row")]
        )
        df_row <- data.frame(Row = levels(sce$Row), row_median = med$row)
        df_col <- data.frame(Col = levels(sce$Col), col_median = med$col)
        df_sce <-
            df_plate_layout_full |>
            left_join(df_row, by = join_by(Row)) |>
            left_join(df_col, by = join_by(Col)) |>
            mutate(overall_median = med$overall) |>
            mutate(
                corrected = feature - overall_median - row_median - col_median
            )

        # subtract the overall, row, and column medians
        corrected[feature_ind, ] <- df_sce$corrected
    }

    # save in sce object
    assay(sce, "corrected") <- corrected
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
#' @param target Name of target variable for prediction
#' @param group Grouping variable for cross-validation, e.g., patient
#' @param interest_level Factor interest level in `target` variable
#' @param reference_level Factor reference level in `target` variable
#' @param weights Weights variable when features are aggregated
#' @param n_threads Number of parallel threads for fitting of models
#' @param assay_type A string specifying the assay
#' @return \code{\link[tibble]{tibble}} data frame
#'
compute_y_hat <- function(feature_name,
                          sce_feature,
                          starts,
                          target,
                          group,
                          interest_level,
                          reference_level,
                          weights,
                          n_threads,
                          assay_type) {
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
#' @importFrom S4Vectors metadata<-
#' @importFrom dplyr bind_cols bind_rows arrange select
#' @importFrom parsnip rand_forest set_mode set_engine
#' @importFrom rsample group_initial_split training testing
#' @importFrom purrr map
#' @importFrom ranger ranger
#' @importFrom stringr str_starts
#' @export
#'
#' @param sce \code{\link[SingleCellExperiment]{SingleCellExperiment}} object
#' @param target Name of target variable for prediction
#' @param group Grouping variable for cross-validation, e.g., patient
#' @param interest_level Factor interest level in `target` variable
#' @param reference_level Factor reference level in `target` variable
#' @param types Vector of strings of feature types
#' @param channels Vector of strings of feature channels
#' @param weights Weights variable when features are aggregated
#' @param n_threads Number of parallel threads for fitting of models
#' @param assay_type A string specifying the assay
#' @param progress Show a progress bar while fitting the models of feature
#'   types and channels
#' @return \code{\link[SingleCellExperiment]{SingleCellExperiment}} object
#'   with the cells of the two levels, the leave-one-out prediction scores in
#'   \code{reducedDim(sce, "prevalidated")}, and the settings of this call in
#'   \code{metadata(sce)$predictLOO}, which \code{\link{calculateStats}} uses
#'   to refit the models in its permutation test
#'
#' @examples
#' set.seed(23)
#' cell_file <- generate_data()
#' sce <- loadData(cell_file)
#' sce <- transformScale(sce)
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
                       target,
                       group,
                       interest_level,
                       reference_level,
                       types = NULL,
                       channels = NULL,
                       weights = NULL,
                       n_threads = 1,
                       assay_type = "tfmfeatures",
                       progress = TRUE) {
    # subset for binary classification
    sce_subset <- sce[, sce[[target]] %in% c(reference_level, interest_level)]

    # remove unused levels
    sce_subset[[target]] <- droplevels(sce_subset[[target]])

    # combine and add to reducedDim slot
    y_hat <- compute_y_hat("all", sce_subset,
        starts = TRUE,
        target, group, interest_level, reference_level,
        weights, n_threads, assay_type
    )
    if (length(types) > 0) {
        y_hat <- bind_cols(
            y_hat,
            purrr::map(types, compute_y_hat, sce_subset,
                starts = TRUE,
                target, group, interest_level, reference_level,
                weights, n_threads, assay_type,
                .progress = progress
            ) |>
                bind_cols()
        )
    }
    if (length(channels) > 0) {
        y_hat <- bind_cols(
            y_hat,
            purrr::map(channels, compute_y_hat, sce_subset,
                starts = FALSE,
                target, group, interest_level, reference_level,
                weights, n_threads, assay_type,
                .progress = progress
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

    # record the settings, which calculateStats() needs to refit the models
    metadata(sce_subset)$predictLOO <- list(
        target = target,
        group = group,
        interest_level = interest_level,
        reference_level = reference_level,
        types = types,
        channels = channels,
        weights = weights,
        n_threads = n_threads,
        assay_type = assay_type
    )
    sce_subset
}

#' Aggregate predicted leave-one-out probabilities over meta variables
#'
#' @importFrom scrapper aggregateAcrossCells
#'
#' @param sce \code{\link[SingleCellExperiment]{SingleCellExperiment}} object
#' @param target Name of target variable for prediction
#' @param group Grouping variable for cross-validation, e.g., patient
#' @param assay_type A string specifying the assay
#' @return \code{data.frame}
#'
aggregateYhat <- function(sce,
                          target,
                          group,
                          assay_type = "tfmfeatures") {
    meta_vars <- c(target, group)

    # get group index
    aggr <- scrapper::aggregateAcrossCells(
        assay(sce, assay_type),
        factors = colData(sce)[, meta_vars]
    )

    # calculate mean score
    idx <- aggr$index
    pre <- reducedDim(sce, "prevalidated")
    agg_pre <- rowsum(pre, idx) / tabulate(idx)

    # combine it with meta vars
    cbind(
        aggr$combinations,
        agg_pre
    ) |> data.frame()
}

#' Plot predicted leave-one-out probabilities
#'
#' @importFrom ggplot2 ggplot aes geom_boxplot geom_jitter ylab facet_wrap
#' @importFrom ggplot2 geom_hline
#' @importFrom tidyr pivot_longer
#' @export
#'
#' @param sce \code{\link[SingleCellExperiment]{SingleCellExperiment}} object
#' @param target Name of target variable for prediction
#' @param group Grouping variable for cross-validation, e.g., patient
#' @param assay_type A string specifying the assay
#' @return \code{\link[ggplot2]{ggplot2}} object
#'
#' @examples
#' set.seed(23)
#' cell_file <- generate_data()
#' sce <- loadData(cell_file)
#' sce <- transformScale(sce)
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
#' plotLOO(sce_single, target = "Drug", group = "Patient")
#'
plotLOO <- function(sce,
                    target,
                    group,
                    assay_type = "tfmfeatures") {
    meta_vars <- c(target, group)
    y_hat <- aggregateYhat(sce, target, group, assay_type)

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

#' Internal function to compute the test statistic of each feature group
#'
#' One-sided t-statistic for higher mean scores in the interest level,
#' computed on the scores averaged per patient and level: paired if every
#' patient has both levels, with pooled variance otherwise.
#'
#' @importFrom dplyr pull bind_rows all_of select
#' @importFrom tidyr pivot_wider
#' @importFrom stats t.test
#'
#' @param sce Output of \code{\link{predictLOO}}
#' @param target Name of target variable for prediction
#' @param group Grouping variable for cross-validation, e.g., patient
#' @param assay_type A string specifying the assay
#' @param interest_level Factor interest level in `target` variable
#' @param reference_level Factor reference level in `target` variable
#' @return \code{data.frame} with one row per feature group
#' @noRd
#'
score_statistics <- function(sce,
                             target,
                             group,
                             assay_type,
                             interest_level,
                             reference_level) {
    meta_vars <- c(target, group)
    y_hat <- aggregateYhat(sce, target, group, assay_type)
    feature_vars <- setdiff(names(y_hat), meta_vars)

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
        result <- tryCatch(
            t.test(x, y,
                paired = paired, var.equal = TRUE,
                alternative = "less"
            ),
            error = function(e) list(statistic = NA_real_, p.value = NA_real_)
        )

        data.frame(
            Feature = feature,
            statistic = -unname(result$statistic),
            pvalue = result$p.value,
            log2FoldChange = log2(mean(y, na.rm = TRUE) / mean(x, na.rm = TRUE))
        )
    }) |> bind_rows()
}

#' Internal function to list relabelings of whole patients
#'
#' @importFrom utils combn
#'
#' @param design Either "paired", "unpaired", or "mixed"
#' @param n Number of patients
#' @param observed Logical vector of patients with the interest level
#' @param exact Enumerate all relabelings instead of sampling them
#' @param n_permutations Number of random relabelings if not exact
#' @param both Logical vector of patients with both levels (mixed design)
#' @return Logical matrix with one row per relabeling, the observed labels
#'   excluded, and one column per patient: patients whose two levels are
#'   swapped (paired design, and patients with both levels in a mixed
#'   design) or patients with the interest level (unpaired design, and
#'   patients with one level in a mixed design)
#' @noRd
#'
patient_relabelings <- function(design, n, observed, exact, n_permutations,
                                both = NULL) {
    if (design == "mixed") {
        return(mixed_relabelings(observed, both, exact, n_permutations))
    }
    if (design == "paired") {
        if (exact) {
            grid <- as.matrix(expand.grid(rep(list(c(FALSE, TRUE)), n)))
            return(unname(grid[-1, , drop = FALSE]))
        }
        return(matrix(
            sample(c(FALSE, TRUE), n_permutations * n, replace = TRUE),
            nrow = n_permutations
        ))
    }
    n_interest <- sum(observed)
    if (exact) {
        subsets <- combn(n, n_interest)
        grid <- t(apply(subsets, 2, function(subset) seq_len(n) %in% subset))
        is_observed <- apply(grid, 1, function(row) all(row == observed))
        return(grid[!is_observed, , drop = FALSE])
    }
    t(replicate(n_permutations, seq_len(n) %in% sample(n, n_interest)))
}

#' Internal function to list all subsets of k out of m patients
#'
#' @importFrom utils combn
#' @return Logical matrix with one row per subset and one column per patient
#' @noRd
#'
all_subsets <- function(m, k) {
    if (k == 0 || k == m) {
        return(matrix(rep(k == m, m), nrow = 1))
    }
    t(apply(combn(m, k), 2, function(subset) seq_len(m) %in% subset))
}

#' Internal function to list relabelings in a mixed design
#'
#' Patients with both levels keep or swap their two levels, and the levels
#' of the patients with one level are permuted among them.
#'
#' @param observed Logical vector of patients with the interest level
#' @param both Logical vector of patients with both levels
#' @param exact Enumerate all relabelings instead of sampling them
#' @param n_permutations Number of random relabelings if not exact
#' @return Logical matrix as returned by \code{patient_relabelings}
#' @noRd
#'
mixed_relabelings <- function(observed, both, exact, n_permutations) {
    single <- !both
    n_both <- sum(both)
    n_single <- sum(single)
    n_interest <- sum(observed[single])
    if (!exact) {
        grid <- matrix(FALSE, nrow = n_permutations, ncol = length(both))
        grid[, both] <- sample(c(FALSE, TRUE), n_permutations * n_both,
            replace = TRUE
        )
        grid[, single] <- matrix(
            replicate(
                n_permutations,
                seq_len(n_single) %in% sample(n_single, n_interest)
            ),
            ncol = n_single, byrow = TRUE
        )
        return(grid)
    }
    swaps <- as.matrix(expand.grid(rep(list(c(FALSE, TRUE)), n_both)))
    assignments <- all_subsets(n_single, n_interest)
    index <- expand.grid(
        swap = seq_len(nrow(swaps)),
        assignment = seq_len(nrow(assignments))
    )
    grid <- matrix(FALSE, nrow = nrow(index), ncol = length(both))
    grid[, both] <- swaps[index$swap, , drop = FALSE]
    grid[, single] <- assignments[index$assignment, , drop = FALSE]
    is_observed <- !apply(grid[, both, drop = FALSE], 1, any) &
        apply(grid[, single, drop = FALSE], 1, function(row) {
            all(row == observed[single])
        })
    grid[!is_observed, , drop = FALSE]
}

#' Internal function to compute permutation p-values by relabeling patients
#'
#' @importFrom cli cli_abort cli_progress_bar cli_progress_update
#' @importFrom cli cli_progress_done
#'
#' @param sce Output of \code{\link{predictLOO}}
#' @param observed Output of \code{score_statistics} for the observed labels
#' @param settings Settings of \code{\link{predictLOO}} from its metadata
#' @param n_permutations Maximal number of relabelings besides the observed
#' @return List with p-values and the number of relabelings
#' @noRd
#'
permutation_pvalues <- function(sce, observed, settings, n_permutations) {
    target <- settings$target
    group <- settings$group
    interest_level <- settings$interest_level
    reference_level <- settings$reference_level
    labels <- as.character(sce[[target]])
    patients <- as.character(sce[[group]])
    ids <- sort(unique(patients))
    n <- length(ids)
    in_interest <- ids %in% patients[labels == interest_level]
    in_reference <- ids %in% patients[labels == reference_level]

    # paired: every patient has both levels; unpaired: every patient has one;
    # mixed: some patients have both levels, the others have one
    both <- in_interest & in_reference
    if (all(both)) {
        design <- "paired"
        n_total <- 2^n
    } else if (all(xor(in_interest, in_reference))) {
        design <- "unpaired"
        n_total <- choose(n, sum(in_interest))
    } else {
        design <- "mixed"
        n_total <- 2^sum(both) * choose(sum(!both), sum(in_interest & !both))
    }
    exact <- n_total <= n_permutations + 1
    relabelings <- patient_relabelings(
        design, n, in_interest, exact, n_permutations, both
    )
    n_relabelings <- nrow(relabelings) + 1
    message(
        "Permutation test with ", if (exact) "all " else "", n_relabelings,
        " relabelings of ", n, " levels of '", group, "' (", design,
        " design); smallest attainable p-value 1/", n_relabelings, "."
    )
    if (design == "mixed") {
        message(
            "Levels of '", group, "' with cells of only one level of '",
            target, "': ", paste(ids[!both], collapse = ", "), "."
        )
    }

    # labels of the cells after relabeling whole patients
    relabel <- function(row) {
        if (design == "unpaired") {
            return(ifelse(patients %in% ids[row],
                interest_level, reference_level
            ))
        }
        # swap the two levels of the selected patients with both levels
        swap <- if (design == "paired") row else row & both
        swapped <- labels
        selected <- patients %in% ids[swap]
        swapped[selected] <- ifelse(labels[selected] == interest_level,
            reference_level, interest_level
        )
        if (design == "mixed") {
            # permute the levels among the patients with one level
            single <- patients %in% ids[!both]
            swapped[single] <- ifelse(patients[single] %in% ids[row & !both],
                interest_level, reference_level
            )
        }
        swapped
    }

    # refit all leave-one-out models for each relabeling
    null_statistics <- matrix(NA_real_,
        nrow = nrow(relabelings), ncol = nrow(observed)
    )
    cli_progress_bar("Refitting models", total = nrow(relabelings))
    for (b in seq_len(nrow(relabelings))) {
        relabeled <- sce
        relabeled[[target]] <- factor(relabel(relabelings[b, ]),
            levels = levels(sce[[target]])
        )
        refit <- do.call(
            predictLOO,
            c(list(sce = relabeled), settings, list(progress = FALSE))
        )
        statistics <- score_statistics(refit, target, group,
            settings$assay_type, interest_level, reference_level
        )
        null_statistics[b, ] <- statistics$statistic[
            match(observed$Feature, statistics$Feature)
        ]
        cli_progress_update()
    }
    cli_progress_done()

    # a relabeling counts against the observed labels if its statistic is at
    # least as large; undefined statistics count against them too
    exceeds <- sweep(null_statistics, 2, observed$statistic, ">=") |
        is.na(null_statistics)
    list(
        pvalue = (1 + colSums(exceeds)) / n_relabelings,
        relabelings = n_relabelings
    )
}

#' Test each feature group for a difference between two levels
#'
#' Average the leave-one-out prediction scores of \code{\link{predictLOO}}
#' per patient and level, and test in each feature group whether the
#' interest level has higher scores than the reference level.
#'
#' The test statistic is the one-sided t-statistic for higher mean scores in
#' the interest level. It is paired if every patient (`group`) has cells of
#' both levels, and uses a pooled variance otherwise.
#'
#' With `method = "permutation"` (the default), the p-value comes from
#' relabeling whole patients: in a paired design, the two levels of a patient
#' are swapped; in an unpaired design, the levels are permuted across
#' patients; and in a mixed design, where only some patients have cells of
#' both levels, the levels of those patients are swapped and the levels of
#' the other patients are permuted among them. For each relabeling, all
#' leave-one-out models are refitted with
#' the settings of the call to \code{\link{predictLOO}}, and the p-value is
#' the proportion of relabelings, the observed one included, with a statistic
#' at least as large. If there are at most `n_permutations + 1` relabelings
#' (\eqn{2^n} for \eqn{n} patients in a paired design), all of them are
#' enumerated and the test is exact; otherwise, `n_permutations` relabelings
#' are drawn at random. The smallest attainable p-value is one over the
#' number of relabelings. Because all models are refitted for every
#' relabeling, the test takes that many times longer than
#' \code{\link{predictLOO}}. Use \code{\link[base]{set.seed}} for
#' reproducible results.
#'
#' With `method = "t"`, the p-value comes from the t-distribution. This is
#' fast but anti-conservative, because the models of different held-out
#' patients share most of their training data, which makes the scores of
#' different patients positively correlated.
#'
#' @importFrom S4Vectors metadata
#' @importFrom cli cli_abort
#' @export
#'
#' @param sce Output of \code{\link{predictLOO}}
#' @param target Name of target variable for prediction
#' @param group Grouping variable for cross-validation, e.g., patient
#' @param assay_type A string specifying the assay
#' @param method Either `"permutation"` (default) for a permutation test that
#'   relabels whole patients, or `"t"` for p-values from the t-distribution
#' @param n_permutations Largest number of relabelings besides the observed
#'   labels; if there are more, this many are drawn at random
#' @return \code{data.frame} with one row per feature group and the columns
#'   `Target` (interest level), `Feature` (feature group), `pvalue`, and
#'   `log2FoldChange` (log2 ratio of the mean scores of the interest and the
#'   reference level). With the permutation test, the attribute
#'   `relabelings` holds the number of relabelings, the observed one
#'   included.
#'
#' @examples
#' set.seed(23)
#' cell_file <- generate_data()
#' sce <- loadData(cell_file)
#' sce <- transformScale(sce)
#'
#' sce$Drug <- as.factor(sce$Drug)
#' sce$Drug <- relevel(sce$Drug, ref = "D1")
#'
#' sce_single <- predictLOO(
#'     sce,
#'     target = "Drug", group = "Patient",
#'     interest_level = "D7", reference_level = "D1",
#'     n_threads = 1
#' )
#'
#' # permutation test; by default, it enumerates all 64 relabelings of the
#' # six patients, here 9 random relabelings keep the example fast
#' calculateStats(sce_single,
#'     target = "Drug", group = "Patient",
#'     n_permutations = 9
#' )
#'
#' # p-values from the t-distribution: fast but anti-conservative
#' calculateStats(sce_single,
#'     target = "Drug", group = "Patient",
#'     method = "t"
#' )
#'
calculateStats <- function(sce,
                           target,
                           group,
                           assay_type = "tfmfeatures",
                           method = c("permutation", "t"),
                           n_permutations = 999) {
    method <- match.arg(method)
    settings <- metadata(sce)$predictLOO

    # levels as passed to predictLOO(); objects from earlier versions of
    # cellpaintr store the reference level first
    if (!is.null(settings) && identical(settings$target, target)) {
        reference_level <- settings$reference_level
        interest_level <- settings$interest_level
    } else {
        reference_level <- levels(sce[[target]])[1]
        interest_level <- levels(sce[[target]])[2]
    }

    observed <- score_statistics(sce, target, group, assay_type,
        interest_level, reference_level
    )
    stats <- data.frame(
        Target = interest_level,
        Feature = observed$Feature,
        pvalue = observed$pvalue,
        log2FoldChange = observed$log2FoldChange
    )
    if (method == "t") {
        return(stats)
    }

    if (is.null(settings) || !identical(settings$target, target) ||
        !identical(settings$group, group)) {
        cli_abort(c(
            "The permutation test refits the models of {.fn predictLOO}.",
            "i" = paste(
                "Pass its output with the same {.arg target} and",
                "{.arg group}, or use {.code method = \"t\"}."
            )
        ))
    }
    if (!is.numeric(n_permutations) || length(n_permutations) != 1 ||
        n_permutations < 1) {
        cli_abort("{.arg n_permutations} must be a positive number.")
    }
    permutation <- permutation_pvalues(sce, observed, settings, n_permutations)
    stats$pvalue <- permutation$pvalue
    attr(stats, "relabelings") <- permutation$relabelings
    stats
}

#' Plot predicted leave-one-out probabilities
#'
#' @importFrom ggplot2 ggplot aes geom_vline geom_hline xlab ylab geom_point
#' @importFrom dplyr mutate
#' @importFrom ggrepel geom_text_repel
#' @export
#'
#' @param sce A \code{\link[SingleCellExperiment]{SingleCellExperiment}} object
#' @param target Name of target variable for prediction
#' @param group Grouping variable for cross-validation, e.g., patient
#' @param p_cutoff Cut-off for statistical significance. A horizontal line
#'                 will be drawn at -log10(p_cutoff). Defaults to 0.01
#'                 divided by the number of feature groups. With few
#'                 patients, permutation p-values cannot fall below this
#'                 default, because the smallest attainable p-value is one
#'                 over the number of relabelings.
#' @param fc_cutoff Cut-off for absolute log2 fold-change. A vertical lines
#'                  will be drawn at fc_cutoff.
#' @param assay_type A string specifying the assay
#' @param stats Output of \code{\link{calculateStats}}; if `NULL`, it is
#'              computed from `sce`
#' @param ... Further arguments to \code{\link{calculateStats}}, such as
#'            `method` or `n_permutations`
#' @return \code{\link[ggplot2]{ggplot2}} object
#'
#' @examples
#' set.seed(23)
#' cell_file <- generate_data()
#' sce <- loadData(cell_file)
#' sce <- transformScale(sce)
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
#' # p-values from the t-distribution keep the example fast; the default
#' # permutation test refits all models for every relabeling of the patients
#' volcanoPlot(sce_single,
#'     target = "Drug", group = "Patient",
#'     p_cutoff = 0.05, fc_cutoff = 0.5, method = "t"
#' )
#'
volcanoPlot <- function(sce,
                        target,
                        group,
                        p_cutoff = NULL,
                        fc_cutoff = 1.0,
                        assay_type = "tfmfeatures",
                        stats = NULL,
                        ...) {
    if (is.null(stats)) {
        stats <- calculateStats(sce, target, group, assay_type, ...)
    }

    if (is.null(p_cutoff)) {
        p_cutoff <- 0.01 / nrow(stats)
    }

    stats |>
        mutate(Feature = ifelse(pvalue < p_cutoff & log2FoldChange > fc_cutoff,
            Feature, ""
        )) |>
        ggplot(aes(log2FoldChange, -log10(pvalue), label = Feature)) +
        geom_vline(
            xintercept = c(0, fc_cutoff), alpha = 0.5,
            linetype = "dashed"
        ) +
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
#' @param target Name of target variable for prediction
#' @param group Grouping variable for cross-validation, e.g., patient
#' @param assay_type A string specifying the assay
#' @return \code{\link[ggplot2]{ggplot2}} object
#'
#' @examples
#' set.seed(23)
#' cell_file <- generate_data()
#' sce <- loadData(cell_file)
#' sce <- transformScale(sce)
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
#' plotROC(sce_single, target = "Drug", group = "Patient")
#'
plotROC <- function(sce,
                    target,
                    group,
                    assay_type = "tfmfeatures") {
    meta_vars <- c(target, group)
    y_hat <- aggregateYhat(sce, target, group, assay_type)
    result <- y_hat |>
        pivot_longer(
            cols = -c(all_of(meta_vars)),
            names_to = "features", values_to = "pred"
        ) |>
        droplevels()
    result[[target]] <- factor(result[[target]], levels = levels(sce[[target]]))

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
#' @export
#'
#' @param sce A \code{\link[SingleCellExperiment]{SingleCellExperiment}} object
#' @param target Name of target variable for prediction
#' @param group Grouping variable for cross-validation, e.g., patient
#' @param assay_type A string specifying the assay
#' @return \code{\link[ggplot2]{ggplot2}} object
#'
#' @examples
#' set.seed(23)
#' cell_file <- generate_data()
#' sce <- loadData(cell_file)
#' sce <- transformScale(sce)
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
#' plotAUC(sce_single, target = "Drug", group = "Patient")
#'
plotAUC <- function(sce,
                    target,
                    group,
                    assay_type = "tfmfeatures") {
    meta_vars <- c(target, group)
    y_hat <- aggregateYhat(sce, target, group, assay_type)
    result <- y_hat |>
        pivot_longer(
            cols = -c(all_of(meta_vars)),
            names_to = "features", values_to = "pred"
        ) |>
        droplevels()
    result[[target]] <- factor(result[[target]], levels = levels(sce[[target]]))

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
