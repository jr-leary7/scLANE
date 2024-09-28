#' Generate a smoothed matrix of gene expression using \code{scLANE} models.
#'
#' @name smoothedCountsMatrix
#' @author Jack Leary
#' @description This function takes as input the output from \code{\link{testDynamic}} and returns the fitted values from each model in a wide format, with one column per-gene and one row-per cell. This matrix can be use as input to cell or gene clustering and / or visualizations such as heatmaps.
#' @import magrittr
#' @importFrom future plan multisession sequential
#' @importFrom purrr map discard map2 reduce
#' @importFrom furrr future_map
#' @importFrom stats setNames
#' @param test.dyn.res The list returned by \code{\link{testDynamic}} - no extra processing required. Defaults to NULL.
#' @param pt A data.frame of pseudotime values for each cell. Defaults to NULL.
#' @param size.factor.offset (Optional) An offset to be used to rescale the fitted values. Can be generated easily with \code{\link{createCellOffset}}. No need to provide if the GEE backend was used. Defaults to NULL.
#' @param genes (Optional) A character vector of genes with which to subset the results. Defaults to NULL.
#' @param log1p.norm A boolean specifying whether the smoothed counts should be log1p-transformed after depth normalization. Defaults to FALSE.
#' @param n.cores If parallel execution is desired, how many cores should be utilized? Defaults to 2.
#' @return A list of matrices of smoothed counts, with each element of the list being a single pseudotime lineage.
#' @seealso \code{\link{testDynamic}}
#' @export
#' @examples
#' data(sim_pseudotime)
#' data(scLANE_models)
#' smoothed_dynamics <- smoothedCountsMatrix(scLANE_models,
#'                                           pt = sim_pseudotime,
#'                                           n.cores = 1L)

smoothedCountsMatrix <- function(test.dyn.res = NULL,
                                 size.factor.offset = NULL,
                                 pt = NULL,
                                 genes = NULL,
                                 log1p.norm = FALSE,
                                 n.cores = 2L) {
  # check inputs
  if (is.null(test.dyn.res) || is.null(pt)) { stop("Please provide the scLANE output from testDynamic().") }
  # set up parallel execution
  if (n.cores > 1L) {
    future::plan(future::multisession, workers = n.cores)
  } else {
    future::plan(future::sequential)
  }
  # extract fitted values for all genes of interest across pseudotime lineages
  if (!is.null(genes)) {
    test.dyn.res <- test.dyn.res[genes]
  } else {
    genes <- names(test.dyn.res)
  }
  lineages <- LETTERS[seq(length(test.dyn.res[[1]]))]
  colnames(pt) <- paste0("Lineage_", lineages)
  lineage_mat_list <- purrr::map(lineages, \(x) {
    lineage_name <- paste0("Lineage_", x)
    fitted_vals_list <- furrr::future_map(test.dyn.res, \(y) y[[lineage_name]]$MARGE_Preds) %>%
                        stats::setNames(names(test.dyn.res)) %>%
                        purrr::discard(rlang::is_na) %>%
                        purrr::discard(rlang::is_null) %>%
                        purrr::discard(\(p) rlang::inherits_only(p, "try-error"))
    fitted_vals_mat <- furrr::future_map(fitted_vals_list, \(z) {
                         if (is.null(size.factor.offset)) {
                           exp(z$marge_link_fit)
                         } else {
                           exp(z$marge_link_fit) * unname(size.factor.offset)[!is.na(pt[, lineage_name])]
                         }
                       }) %>%
                       purrr::reduce(cbind)
    colnames(fitted_vals_mat) <- names(fitted_vals_list)
    if (log1p.norm) {
      fitted_vals_mat <- log1p(fitted_vals_mat)
    }
    return(fitted_vals_mat)
  })
  names(lineage_mat_list) <- paste0("Lineage_", lineages)
  return(lineage_mat_list)
}
