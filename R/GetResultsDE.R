#' A function to tidy the results of \code{scLANE()}.
#'
#' @name GetResultsDE
#' @description This function turns the nested list differential expression results of \code{scLANE()} and turns them into a tidy data.frame, minus the \code{ggplot} objects.
#' @param scLANE.results The nested list returned by \code{scLANE()}. Defaults to NULL
#' @return A data.frame containing differential expression results & modeling statistics for each gene.
#' @export
#' @examples
#' GetResultsDE(scLANE_res_list)

GetResultsDE <- function(scLANE.results = NULL) {
  # check inputs
  if (is.null(scLANE.results)) stop("Please provide a result list.")
  # tidy results
  result_df <- do.call(rbind, lapply(scLANE.results, function(x) data.frame(x[-11])))
  return(result_df)
}
