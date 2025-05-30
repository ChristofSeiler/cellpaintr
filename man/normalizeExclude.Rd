% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/cellpaintr.R
\name{normalizeExclude}
\alias{normalizeExclude}
\title{Normalize by median absolute deviation and exclude zero features}
\usage{
normalizeExclude(sce, plate = "Plate", treatment = "Treatment")
}
\arguments{
\item{sce}{\code{\link[SingleCellExperiment]{SingleCellExperiment}} object}

\item{plate}{Name of the plate variable}

\item{treatment}{Name of the treatment variable}
}
\value{
\code{\link[SingleCellExperiment]{SingleCellExperiment}} object
}
\description{
Normalize by median absolute deviation and exclude zero features
}
