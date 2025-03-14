% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/sortObservations.R
\name{sortObservations}
\alias{sortObservations}
\title{Sort observations by sample ID and pseudotime.}
\usage{
sortObservations(expr.mat = NULL, pt.vec = NULL, id.vec = NULL)
}
\arguments{
\item{expr.mat}{Either a \code{SingleCellExperiment}, \code{Seurat}, or \code{cell_data_set} object from which cell-level metadata can be extracted, or a matrix of integer-valued counts with genes as rows & cells as columns. Defaults to NULL.}

\item{pt.vec}{A vector of pseudotime values used to sort the observations. May contain NA values. Defaults to NULL.}

\item{id.vec}{A vector of subject IDs used to sort the observations. Defaults to NULL.}
}
\value{
An object of the same class as the input \code{expr.mat}, but sorted by sample ID & pseudotime.
}
\description{
Since the GEE & GLMM modes require data to be sorted by sample ID and pseudotime, this function provides a simple way to do so for a range of inputs.
}
\details{
\itemize{
\item If the input is a matrix, it is assumed that the columns are cells - and are named as such - and the rows are genes.
\item If the input is a \code{Seurat} object, sorting requires converting to \code{SingleCellExperiment} object first, then ordering, then converting back to a \code{Seurat} object. Some information might be lost, so it is recommended not to overwrite your original \code{Seurat} object.
}
}
\examples{
data(sim_counts)
data(sim_pseudotime)
sorted_counts <- sortObservations(sim_counts,
                                  pt = sim_pseudotime$PT,
                                  id.vec = sim_counts$subject)
}
\author{
Jack R. Leary
}
