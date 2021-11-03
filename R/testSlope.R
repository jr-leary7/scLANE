#' Test whether a gene's expression is increasing or decreasing over a window of pseudotime.
#'
#' @name testSlope
#' @description This function tests whether the slope of a gene's \eqn{\beta} for pseudotime has a significant effect on expression within a pseudotime window. In short, it tells us whether a gene's expression is changing over an interval or not.
#' @importFrom dplyr %>% arrange mutate
#' @param model.list A list of \code{marge} models. Defaults to NULL.
#' @param expr A matrix or data.frame of raw integer counts. Defaults to NULL.
#' @param pt A data.frame of pseudotime values for each cell. Defaults to NULL.
#' @param adj.method The method used to adjust the \emph{p}-values for each slope. Defaults to "BH" for Benjamini-Hochberg.
#' @param fdr.alpha The FDR statistical significance cutoff used after MHT adjustment. Defaults to 0.05.
#' @return A dataframe containing the genes, breakpoints, and slope \emph{p}-values from each model.
#' @seealso \code{\link{p.adjust}}
#' @export
#' @examples
#'  testSlope(model.list = marge_list, fdr.alpha = 0.01)

testSlope <- function(model.list = NULL, adj.method = "BH") {
  # check inputs
  if (is.null(model.list)) { stop("You forgot to provide a model list to testSlope().") }
  if (!all(as.logical(lapply(model.list, function(x) "glm" %in% class(x))))) { stop("All models must be of class glm.") }
  # loops on loops
  genes <- names(model.list)
  rounded_brkpts <- c()
  brkpts <- c()
  p_vals <- c()
  for (m in seq_along(model.list)) {
    marge_model <- model.list[m][[1]]
    model_breakpoints_rounded <- extractBreakpoints(marge_model, sort = FALSE)
    model_breakpoints <- sapply(model_breakpoints_rounded, function(x) model_df$PT[which.min(abs(model_df$PT - x))])
    coef_pvals <- unname(summary(marge_model)$coefficients[, "Pr(>|z|)"][-1])
    rounded_brkpts <- c(rounded_brkpts, model_breakpoints_rounded)
    brkpts <- c(brkpts, model_breakpoints)
    p_vals <- c(p_vals, coef_pvals)
  }
  slope_df <- data.frame(Gene = genes,
                         Breakpoint = brkpts,
                         Rounded_Breakpoint = rounded_brkpts,
                         P_Value = p_vals) %>%
              arrange(P_Value) %>%
              mutate(Adj_P_Value = p.adjust(P_Value, method = adj.method)) %>%
              arrange(Gene, Breakpoint)
  return(slope_df)
}
