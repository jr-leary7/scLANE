#' A truncation function applied on the predictor variable for knot selection.
#'
#' @name min_span
#' @param X_red : a vector of reduced predictor variable values. Defaults to NULL.
#' @param q : the number of predictor variables used. Defaults to NULL.
#' @param minspan : the set minimum span value. Defaults to \code{round((-log2(-(1 / (q * N)) * log(1 - alpha)) / 2.5))}.
#' @param alpha : see Friedman (1991) equation (43). The default is 0.05.
#' @details  This function selects a minimum span between the knots to mitigate runs of correlated noise in the input data and hence avoiding estimation issues, this equation comes from Friedman (1991) equation 43.
#' @return \code{min_span} returns a vector of truncated predictor variable values.
#' @author Jakub Stoklosa and David I. Warton.
#' @references Friedman, J. (1991). Multivariate adaptive regression splines. \emph{The Annals of Statistics}, \strong{19}, 1--67.
#' @references Stoklosa, J. and Warton, D.I. (2018). A generalized estimating equation approach to multivariate adaptive regression splines. \emph{Journal of Computational and Graphical Statistics}, \strong{27}, 245--253.

min_span <- function(X_red = NULL,
                     q = NULL,
                     minspan = NULL,
                     alpha = 0.05) {
  # check inputs
  if (is.null(X_red) || is.null(q)) { stop("Some inputs to min_span() are missing.") }
  N <- length(X_red)
  x <- sort(X_red)
  minspan <- ifelse(is.null(minspan),
                    round((-log2(-(1 / (q * N)) * log(1 - alpha)) / 2.5)),
                    minspan)
  # run function
  okA <- TRUE
  x_new <- min(x, na.rm = TRUE)
  cc <- 1
  while (okA) {
    if ((cc + minspan) > length(x)) {
      break
      okA <- FALSE
    }
    x_new1 <- x[cc + (minspan + 1)]
    x_new <- c(x_new, x_new1)
    cc <- cc + (minspan + 1)
  }
  res <- unique(x_new)
  return(res)
}
