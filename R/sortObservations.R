#' Sort observations by sample ID and pseudotime.
#'
#' @name sortObservations
#' @author Jack Leary
#' @importFrom dplyr arrange
#' @description Since the GEE & GLMM modes require data to be sorted by sample ID and pseudotime, this function provides a simple way to do so for a range of inputs.
#' @param expr.mat Either a \code{SingleCellExperiment} or \code{Seurat} object from which counts can be extracted, or a matrix of integer-valued counts with genes as rows & cells as columns. Defaults to NULL.
#' @param pt.vec A vector of pseudotime values used to sort the observations. May contain NA values. Defaults to NULL.
#' @param id.vec A vector of subject IDs used to sort the observations. Defaults to NULL.
#' @return An object of the same class as the input \code{expr.mat}, but sorted by sample ID & pseudotime.
#' @details
#' \itemize{
#' \item If the input is a matrix, it is assumed that the columns are cells - and are named as such - and the rows are genes.
#' }
#' @export
#' @examples
#' data(sim_counts)
#' data(sim_pseudotime)
#' sorted_counts <- sortObservations(sim_counts,
#'                                   pt = sim_pseudotime$PT,
#'                                   id.vec = sim_counts$subject)

sortObservations <- function(expr.mat = NULL,
                             pt.vec = NULL,
                             id.vec = NULL) {
  # check inputs
  if (is.null(expr.mat) || is.null(pt) || is.null(id.vec)) { stop("You forgot some inputs to sortObservations().") }
  if (any(is.na(id.vec))) { stop("The subject ID vector must not contain any NA values.") }
  # create table with subject ID and pseudotime
  subj_df <- data.frame(ID = id.vec,
                        PT = pt.vec,
                        Cell = colnames(expr.mat))
  # arrange by subject ID then pseudotime
  subj_df <- dplyr::arrange(subj_df,
                            ID,
                            PT)
  # sort object by cells
  expr.mat <- expr.mat[, subj_df$Cell]
  if (!is.unsorted(colnames(expr.mat))) {
    stop("Sorting has failed.")
  }
  return(expr.mat)
}
