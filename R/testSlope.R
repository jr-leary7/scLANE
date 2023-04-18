#' Test whether a gene is dynamic over a pseudotime interval.
#'
#' @name testSlope
#' @author Jack Leary
#' @description This function tests whether each gene's estimated \eqn{\beta} for pseudotime differs significantly from 0 over each empirically estimated sets of knots / pseudotime interval using a Wald test.
#' @import magrittr
#' @importFrom purrr map_dfr
#' @importFrom dplyr arrange mutate case_when with_groups
#' @importFrom stats p.adjust
#' @param test.dyn.results The list returned by \code{\link{testDynamic}} - no extra processing required. Defaults to NULL.
#' @param p.adj.method The method used to adjust the \emph{p}-values for each slope. Defaults to "holm".
#' @param fdr.cutoff The FDR threshold for determining statistical significance. Defaults to 0.01.
#' @return A dataframe containing the genes, breakpoints, and slope \emph{p}-values from each model.
#' @seealso \code{\link{testDynamic}}
#' @seealso \code{\link[stats]{p.adjust}}
#' @export
#' @examples
#' \dontrun{
#' testSlope(test.dyn.results = gene_stats)
#' testSlope(test.dyn.results = gene_stats, method = "BH", fdr.cutoff = 0.05)
#' }

testSlope <- function(test.dyn.results = NULL,
                      p.adj.method = "holm",
                      fdr.cutoff = 0.01) {
  # check inputs
  if (is.null(test.dyn.results)) { stop("You forgot to provide results from testDynamic() to testSlope().") }
  # create table of results
  slope_df <- purrr::map_dfr(test.dyn.results, \(x) { purrr::map_dfr(x, \(y) data.frame(y["MARGE_Slope_Data"][[1]])) }) %>%
              dplyr::arrange(P_Val) %>%
              dplyr::mutate(P_Val_Adj = stats::p.adjust(P_Val, method = p.adj.method)) %>%
              dplyr::arrange(Gene, Breakpoint) %>%
              dplyr::mutate(Gene_Dynamic_Lineage_Slope = dplyr::case_when(P_Val_Adj < fdr.cutoff ~ 1, TRUE ~ 0)) %>%
              dplyr::with_groups(c(Gene, Lineage),
                                 dplyr::mutate,
                                 Gene_Dynamic_Lineage = dplyr::case_when(any(Gene_Dynamic_Lineage_Slope == 1) ~ 1, TRUE ~ 0)) %>%
              dplyr::with_groups(Gene,
                                 dplyr::mutate,
                                 Gene_Dynamic_Overall = dplyr::case_when(any(Gene_Dynamic_Lineage == 1) ~ 1, TRUE ~ 0))
  return(slope_df)
}
