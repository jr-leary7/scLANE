#' Test whether a gene's expression is increasing or decreasing over a window of pseudotime.
#'
#'  @name testChange
#'  @description This function tests whether the slope of a gene's \eqn{\hat{\bata}} for pseudotime has a significant effect on expression within a pseudotime window. In short, it tells us wether a gene's expression is changing over an interval or not.
#'  @importFrom
#'  @param model.list A list of \code{marge} models. Defaults to NULL.
#'  @param fdr.alpha The FDR statistical significance cutoff used after MHT adjustment. Defaults to 0.05.
#'  @export
#'  @examples
#'  testChange(model.list = marge_list, fdr.alpha = 0.01)

testChange <- function(model.list = NULL, expr = NULL, fdr.alpha = 0.05) {
  # check inputs
  if (is.null(model.list)) { stop("You forgot to provide a model list to testChange().") }
  if (fdr.alpha > 0.1) { warning(sprintf("fdr.alpha is set to %s. Are you sure about this?", fdr.alpha)) }
  if (!all(lapply(model.list), function(x) "glm" %in% class(x))) { stop("All models must be of class glm.") }
  # loops on loops
  for (m in model.list) {
    model_breakpoints <- extractBreakpoints(m)
    for (k in model_breakpoints) {

    }
  }
}
