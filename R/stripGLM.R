#' Make GLM objects much smaller.
#'
#' @name stripGLM
#' @author Jack R. Leary
#' @description This function removes a \emph{lot} of components from the default GLM object in order to make it take up less memory. It does however retain enough pieces for \code{predict()} to still work. No promises beyond that.
#' @param glm.obj An object of class GLM from which you'd like to strip out unnecessary components. Defaults to NULL.
#' @return A slimmed-down \code{glm} object.
#' @seealso \code{\link[MASS]{glm.nb}}

stripGLM <- function(glm.obj = NULL) {
  # check inputs
  if (is.null(glm.obj)) { stop("You forgot to supply inputs to stripGLM().") }
  if (inherits(glm.obj, "try-error")) { return(glm.obj) }
  if (!inherits(glm.obj, "glm")) { stop("Input to stripGLM() must be of class glm.") }

  # strip out unnecessary glm pieces
  glm.obj$effects <- NULL
  glm.obj$weights <- NULL
  glm.obj$residuals <- NULL
  glm.obj$prior.weights <- NULL
  glm.obj$fitted.values <- NULL

  # if glm() wasn't run with model = FALSE and y = FALSE (defaults for scLANE)
  if ("model" %in% names(glm.obj)) {
    glm.obj$model <- NULL
  }
  if ("y" %in% names(glm.obj)) {
    glm.obj$y <- NULL
  }

  return(glm.obj)
}
