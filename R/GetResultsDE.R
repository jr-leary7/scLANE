#' Tidy the results of \code{\link{testDynamic}}.
#'
#' @name getResultsDE
#' @author Jack Leary
#' @description This function turns the nested list differential expression results of \code{\link{testDynamic}} and turns them into a tidy \code{data.frame}.
#' @import magrittr
#' @importFrom purrr map_dfr
#' @importFrom dplyr arrange mutate across case_when with_groups relocate
#' @importFrom tidyselect everything
#' @importFrom stats p.adjust p.adjust.methods
#' @param test.dyn.results The nested list returned by \code{\link{testDynamic}}. Defaults to NULL.
#' @param p.adj.method (Optional) The method used to adjust \emph{p}-values for multiple hypothesis testing. Defaults to "holm".
#' @param fdr.cutoff (Optional) The FDR threshold for determining statistical significance. Defaults to 0.01.
#' @return A data.frame containing differential expression results & modeling statistics for each gene.
#' @export
#' @seealso \code{\link{testDynamic}}
#' @seealso \code{\link[stats]{p.adjust}}
#' @examples
#' \dontrun{
#' getResultsDE(test.dyn.results = gene_stats)
#' getResultsDE(test.dyn.results = gene_stats,
#'              p.adj.method = "BH",
#'              fdr.cutoff = 5e-3)
#' }

getResultsDE <- function(test.dyn.results = NULL,
                         p.adj.method = "holm",
                         fdr.cutoff = 0.01) {
  # check inputs
  if (is.null(test.dyn.results)) { stop("Please provide a result list.") }
  if (!p.adj.method %in% stats::p.adjust.methods) { stop("Please choose a valid p-value adjustment method.") }
  # iterates first over genes, then over lineages per-gene & coerces to final data.frame after unlisting everything
  result_df <- purrr::map_dfr(test.dyn.results,
                              function(x) {
                                purrr::map_dfr(x,
                                               function(y) {
                                                 as.data.frame(rbind(y[c(1:12)])) %>%
                                                   dplyr::mutate(dplyr::across(tidyselect::everything(), \(z) unname(unlist(z))))
                                               })
                              }) %>%
               dplyr::arrange(P_Val) %>%
               dplyr::mutate(P_Val_Adj = stats::p.adjust(P_Val, method = p.adj.method),
                             Gene_Dynamic_Lineage = dplyr::case_when(P_Val_Adj < fdr.cutoff ~ 1, TRUE ~ 0)) %>%
               dplyr::with_groups(Gene,
                                  dplyr::mutate,
                                  Gene_Dynamic_Overall = max(Gene_Dynamic_Lineage)) %>%
               dplyr::relocate(Gene,
                               Lineage,
                               Test_Stat,
                               P_Val)
  return(result_df)
}
