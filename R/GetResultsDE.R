#' Tidy the results of \code{\link{testDynamic}}.
#'
#' @name getResultsDE
#' @author Jack Leary
#' @description This function turns the nested list differential expression results of \code{\link{testDynamic}} and turns them into a tidy data.frame, minus the \code{ggplot} objects.
#' @param test.dyn.results The nested list returned by \code{\link{testDynamic}}. Defaults to NULL.
#' @param p.adj.method (Optional) The method used to adjust \emph{p}-values for MHT. Defaults to "bonferroni".
#' @param fdr.cutoff The FDR threshold for determining statistical significance. Defaults to 0.01.
#' @import magrittr
#' @importFrom purrr map reduce
#' @importFrom dplyr arrange mutate case_when with_groups relocate
#' @importFrom stats p.adjust
#' @return A data.frame containing differential expression results & modeling statistics for each gene.
#' @export
#' @seealso \code{\link{testDynamic}} \code{\link{p.adjust}}
#' @examples
#' \dontrun{GetResultsDE(test.dyn.results = gene_stats)}
#' #' \dontrun{GetResultsDE(test.dyn.results = gene_stats, p.adj.method = "BH")}

getResultsDE <- function(test.dyn.results = NULL,
                         p.adj.method = "bonferroni",
                         fdr.cutoff = 0.01) {
  # check inputs
  if (is.null(test.dyn.results)) stop("Please provide a result list.")
  # take non-dataframe items of nested list structure from testDynamic() and turn them into a dataframe
  result_df <- purrr::map(gene_stats, function(x) { purrr::map(x, function(y) data.frame(y[1:9])) %>% purrr::reduce(rbind) }) %>%
               purrr::reduce(rbind) %>%
               dplyr::arrange(P_Val) %>%
               dplyr::mutate(P_Val_Adj = stats::p.adjust(P_Val, method = p.adj.method),
                             Gene_Dynamic_Lineage = dplyr::case_when(P_Val_Adj < fdr.cutoff ~ 1, TRUE ~ 0)) %>%
               dplyr::with_groups(Gene, dplyr::mutate, Gene_Dynamic_Overall = max(Gene_Dynamic_Lineage)) %>%
               dplyr::relocate(Gene, Lineage, LRT_Stat, P_Val)
  return(result_df)
}
