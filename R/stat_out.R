#' Fits a linear regression model and calculates RSS/GCV measures (used for MARS linear models).
#'
#' @name stat_out
#' @author Jakub Stoklosa
#' @author David I. Warton
#' @author Jack Leary
#' @importFrom stats lm.fit fitted.values
#' @description Calculate the final RSS / GCV for a fitted model.
#' @param Y The response variable. Defaults to NULL.
#' @param B1 The model matrix of predictor variables. Defaults to NULL.
#' @param TSS Total sum of squares. Defaults to NULL.
#' @param GCV.null GCV value for the intercept model. Defaults to NULL.
#' @param pen The set/fixed penalty used for the GCV. Defaults to 2.
#' @return A list consisting of the RSS, RSSq1, GCV1 and GCVq1 values for the fitted model.
#' @references Friedman, J. (1991). Multivariate adaptive regression splines. \emph{The Annals of Statistics}, \strong{19}, 1--67.
#' @references Milborrow, S. (2017a). Notes on the \code{\link[earth]{earth}} package. Package vignette. Available at: \url{http://www.milbo.org/doc/earth-notes.pdf}.
#' @references Milborrow, S. (2017b). \code{\link[earth]{earth}}: Multivariate Adaptive Regression Splines. R package version 4.4.7. Available at \url{https://cran.r-project.org/package=earth.}
#' @references Stoklosa, J. and Warton, D.I. (2018). A generalized estimating equation approach to multivariate adaptive regression splines. \emph{Journal of Computational and Graphical Statistics}, \strong{27}, 245--253.
#' @seealso \code{\link{stat_out_score_glm_null}}
#' @seealso \code{\link[earth]{earth}}

stat_out <- function(Y = NULL,
                     B1 = NULL,
                     TSS = NULL,
                     GCV.null = NULL,
                     pen = 2) {
  # check inputs
  if (is.null(Y) || is.null(B1) || is.null(TSS) || is.null(GCV.null)) { stop("Some of the arguments to stat_out() are missing.") }
  if (GCV.null == 0) { stop("GCV.null argument to stat_out() cannot be 0.") }
  N <- length(Y)
  reg <- try({ stats::lm.fit(B1, Y) }, silent = TRUE)
  if (inherits(reg, "try-error") || any(is.na(reg$coefficients))) {
    RSS1 <- RSSq1 <- GCV1 <- GCVq1 <- NA_real_
  } else {
    df1a <- ncol(B1) + pen * (ncol(B1) - 1) / 2  # This matches the earth() package, SAS and Friedman (1991) penalty.
    RSS1 <- sum((Y - stats::fitted.values(reg))^2)
    RSSq1 <- 1 - RSS1 / TSS
    GCV1 <- RSS1 / (N * (1 - (df1a) / N)^2)
    GCVq1 <- 1 - GCV1 / GCV.null
  }
  res <- list(RSS1 = RSS1,
              RSSq1 = round(RSSq1, 10),
              GCV1 = GCV1,
              GCVq1 = round(GCVq1, 10))
  return(res)
}
