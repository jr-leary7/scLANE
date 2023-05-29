#' Make GLM objects much smaller.
#'
#' @name stripGLM
#' @author Jack Leary
#' @description This function removes a \emph{lot} of components from the default GLM object in order to make it take up less memory. It does however retain enough pieces for \code{predict()} to still work. No promises beyond that.
#' @param glm.obj An object of class GLM from which you'd like to strip out unnecessary components. Defaults to NULL.
#' @return A slimmed-down \code{glm} object.
#' @export
#' @seealso \code{\link{glm}}
#' @examples
#' \dontrun{
#' stripGLM(marge_model)
#' }

stripGLM <- function(glm.obj = NULL) {
  # check inputs
  if (is.null(glm.obj)) { stop("You forgot to supply inputs to stripGLM().") }
  if (!inherits(glm.obj, "glm")) { stop("Input to stripGLM() must be of class glm.") }

  # strip out unnecessary glm pieces
  glm.obj$effects <- c()
  glm.obj$weights <- c()
  glm.obj$residuals <- c()
  glm.obj$prior.weights <- c()
  glm.obj$fitted.values <- c()

  # if glm() wasn't run with model = FALSE and y = FALSE (defaults for scLANE)
  if ("model" %in% names(glm.obj)) {
    glm.obj$model <- c()
  }
  if ("y" %in% names(glm.obj)) {
    glm.obj$y <- c()
  }

  return(glm.obj)
}
