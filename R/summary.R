#' Summary method for scLANE objects.
#'
#' @name summary.scLANE
#' @author Jack R. Leary
#' @importFrom purrr map reduce
#' @importFrom stats p.adjust
#' @param test.dyn.res The nested list returned by \code{\link{testDynamic}}.
#' Defaults to NULL.
#' @return A summary list with aggregated statistics concerning the trajectory
#' DE tests from \code{scLANE}.
#' @export
#' @examples
#' data(scLANE_models)
#' summary(scLANE_models)
summary.scLANE <- function(object, ...) {
    if (!inherits(object, "scLANE")) {
      stop("The input must be an object of class 'scLANE'.")
    }
    summary_stats <- list()
    summary_stats$n_lineages <- length(object[[1]])
    summary_stats$n_genes <- length(object)
    p_values <- purrr::map(object, \(g) {
        purrr::map(g, \(l) {
            l$P_Val
        })
    })
  p_values <- purrr::reduce(purrr::reduce(p_values, c), c)
  adj_p_values <- stats::p.adjust(sort(p_values), method = "fdr")
  summary_stats$n_significant_genes <- sum(adj_p_values < 0.01, na.rm = TRUE)
  summary_stats$mean_adj_p_value <- mean(adj_p_values, na.rm = TRUE)
  if (object[[1]][[1]]$Test_Stat_Type == "LRT") {
    summary_stats$test_type <- "Likelihood Ratio Test"
  } else if (object[[1]][[1]]$Test_Stat_Type == "Wald") {
    summary_stats$test_type <- "Wald"
  } else {
    summary_stats$test_type <- "Score Test"
  }
  class(summary_stats) <- "summary.scLANE"
  return(summary_stats)
}

#' Print method for summary.scLANE objects.
#'
#' @name print.summary.scLANE
#' @author Jack R. Leary
#' @param x An object of class summary.scLANE.
#' @return A printed summary of overall scLANE results
#' @export

print.summary.scLANE <- function(x, ...) {
    cat("Summary of scLANE Model Results\n")
    cat("-------------------------------\n")
    cat("Test used: ", x$test_type, "\n")
    cat("Total genes analyzed: ", x$n_genes, "\n")
    cat("Total lineages analyzed: ", x$n_lineages, "\n")
    cat(
        "Mean FDR-adjusted p-value: ", format(x$mean_adj_p_value, digits = 5),
        "\n"
    )
    cat(
        "Number of significant gene-lineage tests at the alpha = 0.01 level: ", x$n_significant_genes,
        "\n"
    )
    invisible(x)
}
