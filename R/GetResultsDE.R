#' A function to tidy the results of \code{testDynamic()}.
#'
#' @name getResultsDE
#' @description This function turns the nested list differential expression results of \code{testDynamic()} and turns them into a tidy data.frame, minus the \code{ggplot} objects.
#' @param scLANE.results The nested list returned by \code{testDynamic()}. Defaults to NULL.
#' @param p.adj.method (Optional) The method used to adjust \emph{p}-values for MHT. Defaults to "BH", for Benjamini-Hochberg.
#' @importFrom dplyr arrange mutate `%>%`
#' @return A data.frame containing differential expression results & modeling statistics for each gene.
#' @export
#' @seealso \code{\link{testDynamic}} \code{\link{p.adjust}}
#' @examples
#' GetResultsDE(scLANE.results = scLANE_res_list)

getResultsDE <- function(test.results = NULL, p.adj.method = "BH") {
  if (is.null(test.results)) stop("Please provide a result list.")
  result_df <- do.call(rbind, lapply(test.results, function(x) data.frame(x[1:8])))
  result_df <- result_df %>%
               dplyr::arrange(P_Val) %>%
               dplyr::mutate(P_Val_Adj = p.adjust(P_Val, method = p.adj.method)) %>%
               dplyr::relocate(Gene, LRT_Stat, P_Val)
  return(result_df)
}
