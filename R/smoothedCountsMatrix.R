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
#' @param genes (Optional) A character vector of genes with which to subset the results. Defaults to NULL.
#' @param parallel.exec Should \code{furrr} be used to speed up execution? Defaults to TRUE.
#' @param n.cores If parallel execution is desired, how many cores should be utilized? Defaults to 2.
#' @return A list of matrices of smoothed counts, with each element of the list being a single pseudotime lineage.
#' @seealso \code{\link{testDynamic}}
#' @export
#' @examples
#' \dontrun{
#' smoothedCountsMatrix(gene_stats)
#' smoothedCountsMatrix(gene_stats,
#'                      genes = names(gene_stats)[1:100],
#'                      parallel.exec = TRUE,
#'                      n.cores = 4)
#' }

smoothedCountsMatrix <- function(test.dyn.res = NULL,
                                 genes = NULL,
                                 parallel.exec = TRUE,
                                 n.cores = 2) {
  # check inputs
  if (is.null(test.dyn.res)) { stop("Please provide the scLANE output from testDynamic().") }
  # set up parallel execution
  if (parallel.exec) {
    future::plan(future::multisession, workers = n.cores)
  } else {
    future::plan(future::sequential)
  }
  # extract fitted values for all genes of interest
  if (!is.null(genes)) {
    test.dyn.res <- test.dyn.res[genes]
  }
  lineages <- LETTERS[1:length(test.dyn.res[[1]])]
  lineage_mat_list <- purrr::map(lineages, \(x) {
    lineage_name <- paste0("Lineage_", x)
    fitted_vals_list <- furrr::future_map(test.dyn.res, \(y) y[[lineage_name]]$MARGE_Preds) %>%
                        stats::setNames(names(test.dyn.res)) %>%
                        purrr::discard(rlang::is_na) %>%
                        purrr::discard(\(p) rlang::inherits_only(p, "try-error")) %>%
                        purrr::map(\(z) exp(z$marge_link_fit))
    fitted_vals_mat <- purrr::reduce(fitted_vals_list, cbind) %>%
                       stats::setNames(names(fitted_vals_list))
    return(fitted_vals_mat)
  })
  names(lineage_mat_list) <- paste0("Lineage_", lineages)
  return(lineage_mat_list)
}
