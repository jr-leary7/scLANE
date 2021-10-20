#' Test whether a gene's expression is increasing or decreasing over a window of pseudotime.
#'
#'  @name testChange
#'  @description This function tests whether the slope of a gene's \eqn{\hat{\bata}} for pseudotime has a significant effect on expression within a pseudotime window. In short, it tells us wether a gene's expression is changing over an interval or not.
#'  @importFrom
#'  @param model.list A list of \code{marge} models. Defaults to NULL.
#'  @param expr A matrix or data.frame of raw integer counts. Defaults to NULL.
#'  @param pt A vector of pseudotime values for each cell. Defaults to NULL.
#'  @param fdr.alpha The FDR statistical significance cutoff used after MHT adjustment. Defaults to 0.05.
#'  @export
#'  @examples
#'  testChange(model.list = marge_list, fdr.alpha = 0.01)

testChange <- function(model.list = NULL, expr = NULL, pt = NULL, fdr.alpha = 0.05) {
  # check inputs
  if (any(sapply(c(model.list, expr, pt)))) { stop("You forgot to provide a model list to testChange().") }
  if (fdr.alpha > 0.1) { warning(sprintf("fdr.alpha is set to %s. Are you sure about this?", fdr.alpha)) }
  if (!all(lapply(model.list), function(x) "glm" %in% class(x))) { stop("All models must be of class glm.") }
  # loops on loops
  for (m in seq_along(model.list)) {
    marge_model <- model.list[m][[1]]
    model_breakpoints <- extractBreakpoints(marge_model)
    n_breakpoints <- length(model_breakpoints)
    model_df <- data.frame(PT = pt, Expr = expr[, m])
    for (k in seq_along(model_breakpoints)) {
      if (k == 1) {
        # min case
        model_df_sub <- model_df[model_df$PT <= model_breakpoints[k], ]
      } else if (k == n_breakpoints) {
        # max case
        model_df_sub <- model_df[model_df$PT > model_breakpoints[k], ]
      } else {
        # all other cases
        model_df_sub <- model_df[model_df$PT <= model_breakpoints[k] & model_df$PT > model_breakpoints[k - 1], ]
      }
      sub_model <- MASS::glm.nb(Expr ~ PT, data = model_df_sub)
      sub_model_pval <- summary(sub_model)$coefficients["PT", "Pr(>|z|)"]
    }
  }
}
