#' A \code{\link[SingleCellExperiment]{SingleCellExperiment}} object containing simulated counts.
#'
#' Data simulated using the \code{scaffold} R package for 50 dynamic and 50 static genes across 1200 cells from 3 subjects.
#'
#' @format An object of class \code{\link[SingleCellExperiment]{SingleCellExperiment}}.
#' @source https://www.rhondabacher.com/scaffold-vignette.pdf
#' @usage data(sim_counts)
"sim_counts"

#' A data.frame containing ground-truth pseudotime.
#'
#' The true ordering of the 1200 cells contained in \code{sim_counts}.
#'
#' @format An object of class \code{data.frame} with 1200 rows and one variable:
#' \itemize{
#'  \item PT: the true pseudotime (0.0025--1)
#' }
#' @usage data(sim_pseudotime)
"sim_pseudotime"
