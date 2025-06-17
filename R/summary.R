#' Summary method for scLANE objects.
#'
#' @name summary.scLANE
#' @author Jack R. Leary
#' @importFrom purrr map reduce
#' @importFrom stats p.adjust
#' @param test.dyn.res The nested list returned by \code{\link{testDynamic}}. Defaults to NULL.
#' @param p.adj.method (Optional) The method used to adjust \emph{p}-values for multiple hypothesis testing. Defaults to "fdr".
#' @param fdr.cutoff (Optional) The FDR threshold for determining statistical significance. Defaults to 0.01.
#' @return A summary list with aggregated statistics concerning the trajectory DE tests from \code{scLANE}.
#' @export
#' @examples
#' data(scLANE_models)
#' summary(scLANE_models)

summary.scLANE <- function(test.dyn.res = NULL,
                           p.adj.method = "fdr",
                           fdr.cutoff = 0.01) {
  if (!inherits(test.dyn.res, "scLANE")) { stop("The input must be an object of class 'scLANE'.") }
  summary_stats <- list()
  summary_stats$n_lineages <- length(test.dyn.res[[1]])
  summary_stats$n_genes <- length(test.dyn.res)
  p_values <- purrr::map(test.dyn.res, \(g) {
    purrr::map(g, \(l) {
      l$P_Val
    })
  })
  p_values <- purrr::reduce(purrr::reduce(p_values, c), c)
  adj_p_values <- stats::p.adjust(sort(p_values), method = p.adj.method)
  summary_stats$n_significant_genes <- sum(adj_p_values < fdr.cutoff, na.rm = TRUE)
  summary_stats$mean_adj_p_value <- mean(adj_p_values, na.rm = TRUE)
  summary_stats$test_type <- ifelse(test.dyn.res[[1]][[1]]$Test_Stat_Type == "LRT",
                                    "Likelihood Ratio Test",
                                    ifelse(test.dyn.res[[1]][[1]]$Test_Stat_Type == "Wald",
                                           "Wald Test",
                                           "Score Test"))
  class(summary_stats) <- "summary.scLANE"
  return(summary_stats)
}

#' Print method for summary.scLANE objects.
#'
#' @name print.summary.scLANE
#' @author Jack R. Leary
#' @param x An object of class summary.scLANE.
#' @export

print.summary.scLANE <- function(x) {
  cat("Summary of scLANE Model Results\n")
  cat("-------------------------------\n")
  cat("Test used: ", x$test_type, "\n")
  cat("Total genes analyzed: ", x$n_genes, "\n")
  cat("Total lineages analyzed: ", x$n_lineages, "\n")
  cat("Mean adjusted p-value: ", format(x$mean_adj_p_value, digits = 5), "\n")
  cat("Number of significant gene-lineage tests: ", x$n_significant_genes, "\n")
  invisible(x)
}
