#' A function to tidy the results of \code{scLANE()}.
#'
#' @name getResultsDE
#' @description This function turns the nested list differential expression results of \code{scLANE()} and turns them into a tidy data.frame, minus the \code{ggplot} objects.
#' @param scLANE.results The nested list returned by \code{scLANE()}. Defaults to NULL.
#' @param p.adj.method (Optional) The method used to adjust \emph{p}-values for MHT. Defaults to "BH", for Benjamini-Hochberg.
#' @importFrom dplyr arrange mutate `%>%`
#' @return A data.frame containing differential expression results & modeling statistics for each gene.
#' @export
#' @seealso \code{\link{scLANE}} \code{\link{p.adjust}}
#' @examples
#' GetResultsDE(scLANE.results = scLANE_res_list)

getResultsDE <- function(scLANE.results = NULL, p.adj.method = "BH") {
  if (is.null(scLANE.results)) stop("Please provide a result list.")
  result_df <- do.call(rbind, lapply(scLANE.results, function(x) data.frame(x[-c(11:13)])))
  result_df <- result_df %>%
               dplyr::arrange(Boot_P_Val) %>%
               dplyr::mutate(Boot_P_Val_Adj = p.adjust(Boot_P_Val, method = p.adj.method)) %>%
               dplyr::relocate(Gene, LRT_Stat, P_Val, Boot_P_Val, Boot_P_Val_SE, Boot_P_Val_Adj)
  return(result_df)
}
