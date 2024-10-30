#' Tidy the results of \code{\link{testDynamic}}.
#'
#' @name getResultsDE
#' @author Jack R. Leary
#' @description This function turns the nested list differential expression results of \code{\link{testDynamic}} and turns them into a tidy data.frame.
#' @import magrittr
#' @importFrom future plan multisession sequential 
#' @importFrom furrr future_map_dfr
#' @importFrom purrr map_dfr
#' @importFrom dplyr arrange desc mutate across if_else with_groups relocate
#' @importFrom tidyselect everything
#' @importFrom stats p.adjust p.adjust.methods
#' @param test.dyn.res The nested list returned by \code{\link{testDynamic}}. Defaults to NULL.
#' @param p.adj.method (Optional) The method used to adjust \emph{p}-values for multiple hypothesis testing. Defaults to "fdr".
#' @param fdr.cutoff (Optional) The FDR threshold for determining statistical significance. Defaults to 0.01.
#' @param n.cores (Optional) If running in parallel, how many cores should be used? Defaults to 2L.
#' @return A data.frame containing differential expression results & test statistics for each gene.
#' @export
#' @seealso \code{\link{testDynamic}}
#' @seealso \code{\link[stats]{p.adjust}}
#' @examples
#' data(scLANE_models)
#' scLANE_de_res <- getResultsDE(scLANE_models)

getResultsDE <- function(test.dyn.res = NULL,
                         p.adj.method = "fdr",
                         fdr.cutoff = 0.01,
                         n.cores = 2L) {
  # check inputs
  if (is.null(test.dyn.res)) { stop("Please provide a result list from testDynamic().") }
  if (!p.adj.method %in% stats::p.adjust.methods) { stop("Please choose a valid p-value adjustment method.") }
  # set up parallel processing
  if (n.cores > 1L) {
    future::plan(future::multisession, workers = n.cores)
  } else {
    future::plan(future::sequential)
  }
  # iterates first over genes, then over lineages per-gene & coerces to final data.frame after unlisting everything
  result_df <- furrr::future_map_dfr(test.dyn.res,
                                     function(x) {
                                       purrr::map_dfr(x,
                                                      function(y) {
                                                        as.data.frame(rbind(y[seq(15)])) %>%
                                                          dplyr::mutate(dplyr::across(tidyselect::everything(), \(z) unname(unlist(z))))
                                                      })
                                     }) %>%
               dplyr::arrange(dplyr::desc(Test_Stat)) %>%
               dplyr::mutate(P_Val_Adj = stats::p.adjust(P_Val, method = p.adj.method),
                             Gene_Dynamic_Lineage = dplyr::if_else(P_Val_Adj < fdr.cutoff, 1, 0, missing = 0)) %>%
               dplyr::with_groups(Gene,
                                  dplyr::mutate,
                                  Gene_Dynamic_Overall = max(Gene_Dynamic_Lineage)) %>%
               dplyr::relocate(Gene,
                               Lineage,
                               Test_Stat,
                               P_Val)
  # shut down parallel processing
  if (n.cores > 1L) {
    future::plan(future::sequential)
  }
  return(result_df)
}
