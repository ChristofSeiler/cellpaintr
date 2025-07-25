% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/cellpaintr.R
\name{predictYhat}
\alias{predictYhat}
\title{Predict target from features}
\usage{
predictYhat(
  sce,
  assay_type = "tfmfeatures",
  target = "Treatment",
  interest_level = "RBPJ",
  reference_level = "cont",
  types = NULL,
  channels = NULL,
  group = "Patient",
  n_threads = 1,
  n_sub = NULL
)
}
\arguments{
\item{sce}{\code{\link[SingleCellExperiment]{SingleCellExperiment}} object}

\item{assay_type}{A string specifying the assay}

\item{target}{Name of target variable for prediction}

\item{interest_level}{Factor interest level in `target` variable}

\item{reference_level}{Factor reference level in `target` variable}

\item{types}{Vector of strings of feature types}

\item{channels}{Vector of strings of feature channels}

\item{group}{Grouping variable for cross-validation, e.g., patient}

\item{n_threads}{Number of parallel threads for fitting of models}
}
\value{
\code{\link[tibble]{tibble}} data frame
}
\description{
Predict target from features
}
