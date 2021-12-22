#' Test whether a gene's expression is increasing or decreasing over a window of pseudotime.
#'
#' @name testSlope
#' @author Jack Leary
#' @description This function tests whether the slope of a gene's \eqn{\beta} for pseudotime has a significant effect on expression within a pseudotime window. In short, it tells us whether a gene's expression is changing over an interval or not.
#' @import magrittr
#' @importFrom dplyr arrange mutate case_when with_groups
#' @param model.list A list of \code{marge} models. Defaults to NULL.
#' @param pt A data.frame of pseudotime values for each cell. Defaults to NULL.
#' @param adj.method The method used to adjust the \emph{p}-values for each slope. Defaults to "bonferroni".
#' @param fdr.cutoff The FDR threshold for determining statistical significance. Defaults to 0.01.
#' @return A dataframe containing the genes, breakpoints, and slope \emph{p}-values from each model.
#' @seealso \code{\link{p.adjust}}
#' @export
#' @examples
#' \dontrun{testSlope(model.list = marge_list, pt = pt_df)}

testSlope <- function(model.list = NULL,
                      pt = NULL,
                      adj.method = "bonferroni",
                      fdr.cutoff = 0.01) {
  # check inputs
  if (is.null(model.list)) { stop("You forgot to provide a model list to testSlope().") }
  # if (!all(as.logical(lapply(model.list, function(x) "glm" %in% class(x))))) { stop("All models must be of class glm.") }
  # loops on loops
  genes <- names(model.list)
  genes_long <- c()
  rounded_brkpts <- c()
  brkpts <- c()
  brkpt_dirs <- c()
  p_vals <- c()
  for (m in seq_along(model.list)) {
    marge_model <- model.list[m][[1]]
    # checks to see if marge_model was set to NA in testDynamic()
    if (any(class(marge_model)) == "logical") {
      next
    }
    # grab the k breakpoints from the MARGE model
    model_breakpoints_rounded <- extractBreakpoints(marge_model, directions = TRUE)
    model_breakpoints <- sapply(model_breakpoints_rounded$Breakpoint, function(x) pt$PT[which.min(abs(pt$PT - x))])
    # extract Wald test p-values from summary.glm()
    coef_pvals <- unname(summary(marge_model)$coefficients[, "Pr(>|z|)"][-1])
    rounded_brkpts <- c(rounded_brkpts, model_breakpoints_rounded$Breakpoint)
    brkpts <- c(brkpts, model_breakpoints)
    brkpt_dirs <- c(brkpt_dirs, model_breakpoints_rounded$Direction)
    p_vals <- c(p_vals, coef_pvals)
    genes_long <- c(genes_long, rep(genes[m], length(coef_pvals)))
  }
  # create table of results
  slope_df <- data.frame(Gene = genes_long,
                         Breakpoint = unlist(brkpts),
                         Rounded_Breakpoint = rounded_brkpts,
                         Direction = brkpt_dirs,
                         P_Val = p_vals) %>%
              dplyr::arrange(P_Val) %>%
              dplyr::mutate(P_Val_Adj = p.adjust(P_Val, method = adj.method)) %>%
              dplyr::arrange(Gene, Breakpoint) %>%
              dplyr::mutate(P_Val_Adj_Signif = dplyr::case_when(P_Val_Adj < fdr.cutoff ~ 1, TRUE ~ 0)) %>%
              dplyr::with_groups(Gene, dplyr::mutate, Gene_Dynamic = dplyr::case_when(any(P_Val_Adj_Signif == 1) ~ 1, TRUE ~ 0))
  return(slope_df)
}
