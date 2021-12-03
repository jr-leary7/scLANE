#' Test whether a gene's expression is increasing or decreasing over a window of pseudotime.
#'
#' @name testSlope
#' @description This function tests whether the slope of a gene's \eqn{\beta} for pseudotime has a significant effect on expression within a pseudotime window. In short, it tells us whether a gene's expression is changing over an interval or not.
#' @importFrom dplyr %>% arrange mutate case_when with_groups
#' @param model.list A list of \code{marge} models. Defaults to NULL.
#' @param pt A data.frame of pseudotime values for each cell. Defaults to NULL.
#' @param adj.method The method used to adjust the \emph{p}-values for each slope. Defaults to "bonferroni".
#' @param fdr.cutoff The FDR threshhold for determining statistical significance. Defaults to 0.01.
#' @return A dataframe containing the genes, breakpoints, and slope \emph{p}-values from each model.
#' @seealso \code{\link{p.adjust}}
#' @export
#' @examples
#'  testSlope(model.list = marge_list, pt = pt_df)

testSlope <- function(model.list = NULL,
                      pt = NULL,
                      adj.method = "bonferroni",
                      fdr.cutoff = 0.01) {
  # check inputs
  if (is.null(model.list)) { stop("You forgot to provide a model list to testSlope().") }
  if (!all(as.logical(lapply(model.list, function(x) "glm" %in% class(x))))) { stop("All models must be of class glm.") }
  # loops on loops
  genes <- names(model.list)
  rounded_brkpts <- c()
  brkpts <- c()
  brkpt_dirs <- c()
  p_vals <- c()
  for (m in seq_along(model.list)) {
    marge_model <- model.list[m][[1]]
    model_breakpoints_rounded <- extractBreakpoints(marge_model, directions = TRUE)
    model_breakpoints <- sapply(model_breakpoints_rounded$Breakpoint, function(x) pt$PT[which.min(abs(pt$PT - x))])
    coef_pvals <- unname(summary(marge_model)$coefficients[, "Pr(>|z|)"][-1])
    rounded_brkpts <- c(rounded_brkpts, model_breakpoints_rounded$Breakpoint)
    brkpts <- c(brkpts, model_breakpoints)
    brkpt_dirs <- c(brkpt_dirs, model_breakpoints_rounded$Direction)
    p_vals <- c(p_vals, coef_pvals)
  }
  slope_df <- data.frame(Gene = genes,
                         Breakpoint = brkpts,
                         Rounded_Breakpoint = rounded_brkpts,
                         Direction = brkpt_dirs,
                         P_Value = p_vals) %>%
              dplyr::arrange(P_Value) %>%
              dplyr::mutate(Adj_P_Value = p.adjust(P_Value, method = adj.method)) %>%
              dplyr::arrange(Gene, Breakpoint) %>%
              dplyr::mutate(Adj_P_Value_Signif = dplyr::case_when(Adj_P_Value < fdr.cutoff ~ 1, TRUE ~ 0)) %>%
              dplyr::with_groups(Gene, dplyr::mutate, Gene_Dynamic = dplyr::case_when(any(Adj_P_Value_Signif == 1) ~ 1, TRUE ~ 0))
  return(slope_df)
}
