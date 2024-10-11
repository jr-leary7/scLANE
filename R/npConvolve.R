#' Convolution that matches \code{np.convolve}.
#'
#' @name npConvolve
#' @author Jack R. Leary
#' @importFrom stats convolve
#' @description Convolve a vector with a user-specified kernel. Can be useful for heatmap smoothing, weighted moving means, etc.
#' @param x The vector to be convolved. Defaults to NULL.
#' @param conv.kernel The kernel to be used in the convolution. If unspecified, defaults to a vector of \eqn{1 / n} where \eqn{n} is the length of the input vector. Defaults to NULL.
#' @return A convolution with same length as the input vector.
#' @details
#' \itemize{
#' \item The convolution here uses \code{\link[stats]{convolve}}, but creates the kernel and padding in such a way that it matches the output from \code{np.convolve} in Python's \code{numpy} matrix algebra package.
#' }
#' @seealso \code{\link[stats]{convolve}}
#' @export
#' @examples
#' convolved_vec <- npConvolve(x = rnorm(20), conv.kernel = rep(1/5, 5))
#'

npConvolve <- function(x = NULL, conv.kernel = NULL) {
  # check inputs & generate kernel if necessary
  if (is.null(x)) { stop("Cannot convolve with no input.") }
  n <- length(x)
  if (is.null(conv.kernel)) {
    conv.kernel <- rep(1, n) / n
  }
  m <- length(conv.kernel)
  # convolution
  padded_x <- c(rep(0, m %/% 2), x, rep(0, m %/% 2))
  padded_conv <- stats::convolve(padded_x,
                                 conv.kernel,
                                 type = "open")
  start_idx <- (length(padded_conv) - n) %/% 2 + 1
  end_idx <- start_idx + n - 1
  filtered_conv <- padded_conv[start_idx:end_idx]
  # check convolution output
  if (length(filtered_conv) != n) {
    warning("Filtered convolution does not match original vector length, returning unfiltered convolution.")
    res <- padded_conv
  } else {
    res <- filtered_conv
  }
  return(res)
}
