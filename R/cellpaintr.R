#' Load cell painting data from file and convert to a SingleCellExperiment
#'
#' @import SingleCellExperiment
#' @import dplyr
#' @importFrom readr read_csv
#' @importFrom stringr str_detect str_remove
#' @export
#'
#' @param path meta_file
#' @param path features_file
#' @return \code{\link[SingleCellExperiment]{SingleCellExperiment}} object
#'
loadData <- function(cell_file) {

  # read file with cell data
  data      <- read_csv(cell_file, show_col_types = FALSE)

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
  meta_cell <- data[,intersect(names(data), meta_vars)]
  features  <- data[,feature_vars] |> t() # columns represent cells

  # add NAs to location if it doesn't exist
  if(sum(names(meta_cell) == "Location_Center_X") == 0)
    meta_cell$Location_Center_X <- NA
  if(sum(names(meta_cell) == "Location_Center_Y") == 0)
    meta_cell$Location_Center_Y <- NA
  if(sum(names(meta_cell) == "Location_Center_Z") == 0)
    meta_cell$Location_Center_Z <- NA

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
#' @import SingleCellExperiment
#' @import ggplot2
#' @import dplyr
#' @importFrom tibble as_tibble
#' @export
#'
#' @param sce \code{\link[SingleCellExperiment]{SingleCellExperiment}} object
#' @param bins Number of histogram bins
#' @return \code{\link[ggplot2]{ggplot2}} object
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
#' @import SingleCellExperiment
#' @import ggplot2
#' @import dplyr
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
plotPCACor <- function(sce, filter_by = 1, top = 20, pcs = 1:5) {

  features <- assay(sce, "tfmfeatures")
  scores <- reducedDim(sce, "PCA")[,pcs]
  features_pcs <- cor(t(features), scores)

  keep <- rownames(features_pcs)
  importance <- abs(features_pcs[,filter_by])

  features_pcs <- features_pcs |>
    as_tibble() |>
    mutate(name = keep, importance = importance) |>
    top_n(top, importance)

  features_pcs <- features_pcs |>
    pivot_longer(!c(name, importance), names_to = "PC", values_to = "cor") |>
    mutate(PC = str_remove(PC, "PC"))

  ggplot(features_pcs, aes(PC, name, fill = cor)) +
    geom_tile(color = "white") +
    scale_fill_gradient2(low = "blue", high = "red", mid = "white",
                         midpoint = 0, limit = c(-1,1), space = "Lab",
                         name = "Pearson\nCorrelation") +
    theme_minimal() +
    coord_fixed() +
    ylab("feature name")

}

#' Remove cells if not enough or too many in one image
#'
#' @import SingleCellExperiment
#' @import dplyr
#' @importFrom tibble as_tibble
#' @export
#'
#' @param sce \code{\link[SingleCellExperiment]{SingleCellExperiment}} object
#' @param min Remove cells below that number
#' @param max Remove cells above that number
#' @return \code{\link[SingleCellExperiment]{SingleCellExperiment}} object
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

#' Remove missing values
#'
#' @import SingleCellExperiment
#' @export
#'
#' @param sce \code{\link[SingleCellExperiment]{SingleCellExperiment}} object
#' @return \code{\link[SingleCellExperiment]{SingleCellExperiment}} object
#'
removeMissingValues <- function(sce) {

  mat <- assay(sce, name = "features")

  # remove missing values in features
  miss_perc_features <- apply(mat, 1, function(x) mean(is.na(x)))
  feature_ids <- which(miss_perc_features > 0)

  # remove missing values in cells
  miss_perc_cells <- apply(mat, 2, function(x) mean(is.na(x)))
  cell_ids <- which(miss_perc_cells > 0)

  sce[-feature_ids, -cell_ids]

}

#' Filter low variance features
#'
#' @import SingleCellExperiment
#' @export
#'
#' @param sce \code{\link[SingleCellExperiment]{SingleCellExperiment}} object
#' @return \code{\link[SingleCellExperiment]{SingleCellExperiment}} object
#'
removeLowVariance <- function(sce) {

  mat <- assay(sce, name = "features")

  var_features <- apply(mat, 1, var)
  feature_ids <- which(var_features == 0)

  sce[-feature_ids, ]

}

#' Filter low variance features
#'
#' @import SingleCellExperiment
#' @export
#'
#' @param sce \code{\link[SingleCellExperiment]{SingleCellExperiment}} object
#' @return \code{\link[SingleCellExperiment]{SingleCellExperiment}} object
#'
transformLogScale <- function(sce) {

  mat <- assay(sce, name = "features")

  # log(x+1) transform on non-negative valued features
  non_neg_features <- apply(mat, 1, function(x) sum(x >= 0) == length(x) )
  mat[non_neg_features, ] <- log(1 + mat[non_neg_features, ])

  # center and standardize
  mat <- mat |> t() |> scale() |> t()

  # add transformed features
  assays(sce)$tfmfeatures <- mat
  sce

}

#' Normalize by median absolute deviation and exclude zero features
#'
#' @import SingleCellExperiment
#' @export
#'
#' @param sce \code{\link[SingleCellExperiment]{SingleCellExperiment}} object
#' @param plate Name of the plate variable
#' @param treatment Name of the treatment variable
#' @return \code{\link[SingleCellExperiment]{SingleCellExperiment}} object
#'
normalizeExclude <- function(sce, plate = "Plate", treatment = "Treatment") {

  plates <- unique(sce[[plate]])

  plate_list <- lapply(plates, function(plate_id) {

    # normalize per plate
    sce_plate <- sce[, sce[[plate]] == plate_id]

    # calculate median and mad on control cells
    ref_level <- levels(sce_plate[[treatment]])[1]
    sce_ref <- sce_plate[, sce_plate[[treatment]] == ref_level]
    mat_ref <- assay(sce_ref, name = "features")
    median_features <- apply(mat_ref, 1, function(x) median(x) )
    mad_features <- apply(mat_ref, 1, function(x) mad(x) )

    # normalize
    mat <- assay(sce_plate, name = "features")
    mat <- (mat-median_features)/mad_features

    # add transformed features
    assays(sce_plate)$tfmfeatures <- mat

    # exclude
    nonzero_features <- which(median_features != 0)
    sce_plate[nonzero_features, ]

  })

  # exclude features that have zero median in at least one plate
  feature_name_list <- lapply(plate_list,
                              function(sce_plate) rownames(sce_plate))
  include_features <- Reduce(intersect, feature_name_list)
  plate_list <- lapply(plate_list,
                       function(sce_plate) sce_plate[include_features, ])
  do.call(cbind, plate_list)

}

#' Predict target from features
#'
#' @import SingleCellExperiment
#' @import dplyr
#' @importFrom parsnip rand_forest set_mode set_engine
#' @importFrom rsample group_initial_split training testing
#' @export
#'
#' @param sce \code{\link[SingleCellExperiment]{SingleCellExperiment}} object
#' @param assay_type A string specifying the assay
#' @param target Name of target variable for prediction
#' @param interest_level Factor interest level in `target` variable
#' @param reference_level Factor reference level in `target` variable
#' @param group Grouping variable for cross-validation, e.g., patient
#' @param strata Stratification variable for cross-validation, e.g., joint
#' @param n_folds Number of cross-validation folds
#' @param n_threads Number of parallel threads for fitting of models
#' @return \code{\link[tibble]{tibble}} data frame
#'
fitModel <- function(sce, assay_type = "tfmfeatures",
                     target, interest_level, reference_level,
                     group, strata,
                     n_folds = 20, n_threads = 1) {

  # subset for binary classification
  sce_subset <- sce[, sce[[target]] %in% c(reference_level, interest_level)]

  # create matching variable
  cells <- colData(sce_subset) |> as_tibble()
  Match <- apply(cells[, strata], 1,
                 function(row) paste(row, collapse = "_"))
  cells$Match <- Match
  cells$.cell <- 1:nrow(cells)

  # prepare data frame with target variable and predictor matrix
  X <- assay(sce_subset, name = assay_type) |> t()
  ml_data <- data.frame(Target = cells |> pull(target), X)
  ml_data$Target <- factor(ml_data$Target,
                           levels = c(interest_level, reference_level))

  # cross-validation
  rf_list <- lapply(1:n_folds, function(i) {

    system(paste("echo fold = ", i))

    # split into train and test sets
    cells_split <- rsample::group_initial_split(
      cells,
      prop = 1/2,
      group = group,
      strata = "Match"
    )

    # prepare train data
    train_ids <- rsample::training(cells_split) |>
      ungroup() |>
      pull(.cell) |>
      as.integer()

    # prepare test data
    test_ids <- rsample::testing(cells_split) |>
      ungroup() |>
      pull(.cell) |>
      as.integer()

    # keep track of samples info
    train_info <- cells[train_ids, ] |>
      count(across(c(target, group, strata)))
    test_info <- cells[test_ids, ] |>
      count(across(c(target, group, strata)))

    # train model
    rf_spec <-
      parsnip::rand_forest() |>
      parsnip::set_mode("classification") |>
      parsnip::set_engine("ranger", num.threads = n_threads)
    rf_fit <- rf_spec |>
      parsnip::fit(Target ~ ., data = ml_data[train_ids, ])

    # evaluate model
    dplyr::bind_cols(
      .target = ml_data[test_ids, ]$Target,
      predict(rf_fit, ml_data[test_ids, ]),
      predict(rf_fit, ml_data[test_ids, ], type = "prob"),
      fold = i
    ) |> mutate(
      train_info = list(train_info),
      test_info = list(test_info)
    )

  })

  dplyr::bind_rows(rf_list)

}

#' Plot ROC curves
#'
#' @import ggplot2
#' @import dplyr
#' @importFrom yardstick roc_auc roc_curve
#' @importFrom stringr str_remove
#' @export
#'
#' @param result \code{\link[tibble]{tibble}} data frame
#' @return \code{\link[ggplot2]{ggplot2}} object
#'
plotROC <- function(result) {

  mean_auc <- result |>
    group_by(fold) |>
    roc_auc(names(result)[1], names(result)[3]) |>
    pull(.estimate) |>
    mean() |>
    round(digits = 3)

  result |>
    group_by(fold) |>
    roc_curve(names(result)[1], names(result)[3]) |>
    ggplot(aes(x = 1 - specificity, y = sensitivity, group = fold)) +
    geom_path(alpha = 0.2) +
    geom_abline(lty = 3) +
    coord_equal() +
    xlab("false positive rate (1 - specificity)") +
    ylab("true positive rate (sensitivity)") +
    ggtitle(
      paste0("Predict ",
             str_remove(names(result)[3], ".pred_"),
             " vs ",
             str_remove(names(result)[4], ".pred_"),
             " (mean AUC = ",
             mean_auc,
             ")")
    )

}

#' Print sample info of the ROC curves
#'
#' @import dplyr
#' @export
#'
#' @param result \code{\link[tibble]{tibble}} data frame
#' @param fold_id Fold id to print
#' @return NULL.
#'
printROC <- function(result, fold_id) {

  info <- result |>
    filter(fold == fold_id) |>
    slice(1)

  train_info <- info[["train_info"]][[1]]
  test_info <- info[["test_info"]][[1]]

  cat("-----------------\n")
  cat("training samples:\n\n")
  print(train_info)

  cat("\n")
  cat("testing samples:\n\n")
  print(test_info)

  cat("\n")
  cat("class balance:\n\n")

  treatment_var <- names(test_info)[1]
  train_n <- sum(train_info$n)
  test_n <- sum(test_info$n)

  print(
    left_join(
      train_info |>
        group_by(across(treatment_var)) |>
        summarize(train = sum(n)/train_n),
      test_info |>
        group_by(across(treatment_var)) |>
        summarize(test = sum(n)/test_n),
      by = treatment_var
    )
  )

  cat("-----------------")

}

#' Plot AUC comparison
#'
#' @import ggplot2
#' @import dplyr
#' @importFrom yardstick roc_auc
#' @importFrom stringr str_remove
#' @export
#'
#' @param result_list \code{\link[list]{list}} of data frames from `fitModel`
#' @return \code{\link[ggplot2]{ggplot2}} object
#'
plotAUC <- function(result_list) {

  aucs <- lapply(result_list, function(result) {

    level <- names(result)[3] |> stringr::str_remove(".pred_")
    result |>
      group_by(fold) |>
      yardstick::roc_auc(names(result)[1], names(result)[3]) |>
      mutate(predict = level)

  }) |> dplyr::bind_rows()

  fct_order <- aucs |>
    group_by(predict) |>
    summarize(median = median(.estimate)) |>
    arrange(desc(median)) |>
    pull(predict)

  aucs |>
    mutate(predict = factor(predict, levels = fct_order)) |>
    ggplot(aes(predict, .estimate)) +
    geom_boxplot(outlier.shape = NA, width = 0.2) +
    geom_jitter(width = 0.1) +
    ylab("AUC") +
    ggtitle("Area under the ROC curves")

}

#' Evaluate treatment effect
#'
#' @import SingleCellExperiment
#' @import dplyr
#' @import ggplot2
#' @importFrom BiocParallel MulticoreParam bplapply
#' @importFrom tidyr pivot_longer
#' @export
#'
#' @param sce \code{\link[SingleCellExperiment]{SingleCellExperiment}} object
#' @param assay_type A string specifying the assay
#' @param group Grouping variable, e.g., patient
#' @param treatment Treatment variable
#' @param n_perm Number of permutations for the null distribution
#' @param n_threads Number of parallel threads for fitting of models
#' @return \code{\link[ggplot2]{ggplot2}} data frame
#'
treatmentEffect <- function(sce, assay_type = "tfmfeatures", group, treatment,
                            n_perm = 100, n_threads = 1) {

  group_ids <- unique(sce[[group]])

  # treatment label shuffle
  calculate_scores <- function(sce) {

    # calculate perturbation score
    sce_aggr <- aggregateAcrossCells(
      sce,
      use.assay.type = assay_type,
      statistics = "mean",
      ids = colData(sce)[ , c(group, treatment)]
    )
    weights <- attr(reducedDim(sce, "PCA"), "varExplained")
    scores_aggr <- reducedDim(sce_aggr, "PCA")
    scores_weighted <- sweep(scores_aggr, 2, weights, FUN = "*")
    dist_mat <- as.matrix(dist(scores_weighted))

    # extract distances per group
    pert_scores <- sapply(group_ids, function(current_group) {
      patient_pair <- which(sce_aggr[[group]] == current_group)
      dist_mat[patient_pair[1], patient_pair[2]]
    })
    pert_scores

  }

  permute <- function() {

    sce[[treatment]] <- sample(sce[[treatment]], replace = FALSE)
    calculate_scores(sce)

  }

  # compute null distribution of perturbation distances
  param <- MulticoreParam(workers = n_threads, progressbar = TRUE)
  pert_null <- bplapply(seq(n_perm), function(i) permute(), BPPARAM = param)
  pert_null <- simplify2array(pert_null)
  if(is.vector(pert_null)) {
    # if group variables has only one level,
    # then pert_null is a vector
    pert_null <- data.frame(seq(n_perm), pert_null)
  } else {
    # if group variables has more than one level,
    # then pert_null is a matrix
    pert_null <- data.frame(seq(n_perm), t(pert_null))
  }
  names(pert_null) <- c("perm_id", group_ids)

  # compute observed perturbations distances
  pert_obsv <- calculate_scores(sce)
  pert_obsv <- data.frame(name = group_ids, value = pert_obsv)

  # relate null distribution to observed distances
  pivot_longer(pert_null, cols = -perm_id) |>
    ggplot(aes(value)) +
    geom_density() +
    geom_vline(aes(xintercept = value), pert_obsv, colour = "red") +
    facet_wrap(~name) +
    xlab("perturbation distance")

}

#' Predict target from features
#'
#' @import SingleCellExperiment
#' @import dplyr
#' @importFrom parsnip rand_forest set_mode set_engine
#' @importFrom rsample group_initial_split training testing
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
#' @param n_threads Number of parallel threads for fitting of models
#' @return \code{\link[tibble]{tibble}} data frame
#'
predictYhat <- function(sce, assay_type = "tfmfeatures", target = "Treatment",
                        interest_level = "RBPJ", reference_level = "cont",
                        types = NULL, channels = NULL,
                        group = "Patient", n_threads = 1, n_sub = NULL) {

  # subset for binary classification
  sce_subset <- sce[, sce[[target]] %in% c(reference_level, interest_level)]
  if(!is.null(n_sub)) {
    sce_subset <- sce_subset[,sample(ncol(sce_subset), n_sub)]
  }

  # function to compute y_hat on a subset of the features
  compute_y_hat <- function(feature_name, sce_feature, starts = TRUE) {

    # convert to regular expression
    if(feature_name == "all") {
      pattern <- ".*"
    } else {
      pattern <- feature_name
    }

    # subset features
    features <- names(sce_feature)
    if(starts) {
      features <- features[str_starts(features, pattern)]
    } else {
      features <- features[str_detect(features, pattern)]
    }
    sce_feature <- sce_feature[features,]

    # prepare data frame with target variable and predictor matrix
    X <- assay(sce_feature, name = assay_type) |> t()
    cells <- colData(sce_feature) |> as_tibble()
    ml_data <- data.frame(Target = cells |> pull(target), X)
    ml_data$Target <- factor(ml_data$Target,
                             levels = c(interest_level, reference_level))

    # compute y hats
    patients <- unique(sce_feature[[group]])
    result_list <- lapply(patients, function(patient) {

      train_ids <- which(sce_feature[[group]] != patient)
      test_ids <- which(sce_feature[[group]] == patient)

      # train model
      rf_spec <-
        parsnip::rand_forest() |>
        parsnip::set_mode("classification") |>
        parsnip::set_engine("ranger", num.threads = n_threads)
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

  # combine and add to reducedDim slot
  y_hat <- compute_y_hat("all", sce_subset)
  if(length(types) > 0) {
    y_hat <- bind_cols(
      y_hat,
      lapply(types, compute_y_hat, sce_subset, starts = TRUE) |>
        bind_cols()
    )
  }
  if(length(channels) > 0) {
    y_hat <- bind_cols(
      y_hat,
      lapply(channels, compute_y_hat, sce_subset, starts = FALSE) |>
        bind_cols()
    )
  }
  y_hat <- y_hat |> as.matrix()

  # add y_hat to sce object
  cell_id <- 1:ncol(sce_subset)
  rownames(y_hat) <- cell_id
  colnames(sce_subset) <- cell_id
  reducedDim(sce_subset, "prevalidated") <- y_hat
  sce_subset

}

#' Aggregate predicted leave-one-out probabilities over meta variables
#'
#' @importFrom scater aggregateAcrossCells
#' @export
#'
#' @param sce \code{\link[SingleCellExperiment]{SingleCellExperiment}} object
#' @param assay_type A string specifying the assay
#' @param meta_vars a vector of variables from `colData`
#' @return \code{\link[data.frame]{data.frame}}
#'
aggregateYhat <- function(sce,
                          assay_type = "tfmfeatures",
                          meta_vars = c("Patient", "Treatment", "Gender")) {

  summed <- aggregateAcrossCells(
    sce, id = colData(sce)[, meta_vars],
    use.assay.type = assay_type, statistics = "mean"
  )

  cbind(
    colData(summed)[, meta_vars] |> as.data.frame(),
    reducedDim(summed, type = "prevalidated")
  )

}

#' Plot predicted leave-one-out probabilities
#'
#' @import ggplot2
#' @import dplyr
#' @importFrom tidyr pivot_longer
#' @export
#'
#' @param merged a \code{\link[data.frame]{data.frame}} from `aggregateYhat`
#' @param target Name of target variable for prediction
#' @param meta_vars a vector of variables from `colData`
#' @return \code{\link[ggplot2]{ggplot2}} object
#'
plotYhat <- function(y_hat,
                     target = "Treatment",
                     meta_vars = c("Patient", "Treatment", "Gender")) {

  y_hat |>
    pivot_longer(cols = -c(all_of(meta_vars)),
                 names_to = "features", values_to = "value") |>
    ggplot(aes(.data[[target]], value, color = .data[[target]])) +
    geom_boxplot(width = 0.2, outliers = FALSE) +
    geom_jitter(width = 0.1) +
    ylab("predicted leave-one-out probability") +
    facet_wrap(~features) +
    geom_hline(yintercept = 0.5, alpha = 0.5, linetype = "dashed")

}
