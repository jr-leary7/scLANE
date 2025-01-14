#' Use a Wald test to compare nested GEE models.
#'
#' @name waldTestGEE
#' @author Jack R. Leary
#' @description Performs a basic Wald test to determine whether an alternate model is significantly better than a nested null model. This is the GEE equivalent (kind of) of \code{\link{modelLRT}}. Be careful with small sample sizes.
#' @importFrom stats pchisq
#' @param mod.1 The model under the alternative hypothesis. Must be of class \code{geem}. Defaults to NULL.
#' @param mod.0 The model under the null hypothesis. Must be of class \code{geem}. Defaults to NULL.
#' @param correction.method A string specifying the correction method to be used. Currently supported options are "df", "kc", and NULL. Defaults to NULL.
#' @param id.vec A vector of subject IDs. Required when \code{correction.method} is "kc". Defaults to NULL.
#' @param verbose (Optional) A Boolean specifying whether or not verbose output should be printed to the console. Occasionally useful for debugging. Defaults to FALSE.
#' @return A list containing the Wald test statistic, a \emph{p}-value, and the degrees of freedom used in the test.
#' @details
#' \itemize{
#' \item Calculating the test statistic involves taking the inverse of the variance-covariance matrix of the coefficients. Ideally this would be done using the "true" inverse with something like \code{\link{solve}}, \code{\link{qr.solve}}, or \code{\link{chol2inv}}, but in practice this can cause issues when the variance-covariance matrix is near-singular. With this in mind, we use the Moore-Penrose pseudoinverse as implemented in \code{\link[MASS]{ginv}} instead, which leads to more stable results.
#' \item The \emph{p}-value is calculated using an asymptotic \eqn{\Chi^2} distribution, with the degrees of freedom equal to the number of non-intercept coefficients in the alternative model.
#' }
#' @seealso \code{\link[multgee]{waldts}}
#' @seealso \code{\link[aod]{wald.test}}
#' @seealso \code{\link[geeM]{geem}}
#' @seealso \code{\link{biasCorrectGEE}}
#' @seealso \code{\link{modelLRT}}

waldTestGEE <- function(mod.1 = NULL,
                        mod.0 = NULL,
                        correction.method = NULL,
                        id.vec = NULL,
                        verbose = FALSE) {
  # check inputs
  if (is.null(mod.1) || is.null(mod.0)) { stop("Please provide both models to waldTestGEE().") }
  if (inherits(mod.1, "try-error") || inherits(mod.0, "try-error")) {
    res <- list(Wald_Stat = NA_real_,
                DF = NA_real_,
                P_Val = NA_real_,
                Notes = "No test performed due to model failure.")
    return(res)
  }
  if (!is.null(correction.method)) {
    correction.method <- tolower(correction.method)
  }
  if (!is.null(correction.method)) {
    if (!(correction.method %in% c("df", "kc"))) { stop("Unsupported bias correction method in waldTestGEE().") }
    if (correction.method == "kc" && is.null(id.vec)) { stop("The Kauermann and Carroll bias correction method requires a vector of subject IDs.") }
  }
  mod.1 <- mod.1$final_mod
  if (!(inherits(mod.1, "geem") && inherits(mod.0, "geem"))) { stop("You must provide two geeM models to waldTestGee().") }
  if (length(coef(mod.0)) != 1) { stop("Null GEE model must be intercept-only.") }
  if (length(coef(mod.1)) <= length(coef(mod.0))) {
    # can't calculate Wald statistic if both models are intercept-only
    res <- list(Wald_Stat = 0,
                DF = 0,
                P_Val = 1,
                Notes = NA_character_)
  } else {
    # compute test statistic, optionally bias-adjust, & estimate asymptotic p-value
    coef_vals <- as.matrix(coef(mod.1))
    if (!is.null(correction.method)) {
      vcov_mat <- as.matrix(mod.1$var)
      vcov_mat <- biasCorrectGEE(mod.1,
                                 correction.method = correction.method,
                                 id.vec = id.vec,
                                 cor.structure = mod.1$corr,
                                 verbose = verbose)
    } else {
      vcov_mat <- as.matrix(mod.1$naiv.var)
    }
    p_alt <- length(coef(mod.1))
    Lpmat <-  diag(1, nrow = p_alt, ncol = p_alt)
    Lpmat <- Lpmat[-1,,drop=FALSE]
    Lmat <- t(Lpmat)
    middle <- Lpmat %*% vcov_mat %*% Lmat
    middle_inv <- try({ eigenMapMatrixInvert(middle) }, silent = TRUE)
    if (inherits(middle_inv, "try-error")) {
      middle_inv <- eigenMapPseudoInverse(middle)
    }
    sides <- Lpmat %*% coef(mod.1)
    wald_test_stat <- t(sides) %*% middle_inv %*% sides
    if (inherits(wald_test_stat, "try-error")) {
      wald_note <- wald_test_stat[1]  # this is the error message
      wald_test_stat <- 0
      p_value <- 1
    } else {
      p_value <- as.numeric(1 - stats::pchisq(wald_test_stat, df = p_alt - 1))
      wald_note <- NA_character_
    }
    # format results
    res <- list(Wald_Stat = wald_test_stat,
                DF = p_alt - 1,
                P_Val = p_value,
                Notes = wald_note)
  }
  return(res)
}
