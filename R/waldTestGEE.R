#' Use a Wald test to compare nested GEE models.
#'
#' @name waldTestGEE
#' @author Jack Leary
#' @description Performs a basic Wald test to determine whether an alternate model is significantly better than a nested null model. This is the GEE equivalent (kind of) of \code{\link{modelLRT}}. Be careful with small sample sizes.
#' @importFrom MASS ginv
#' @importFrom stats pchisq
#' @param mod.1 The model under the alternative hypothesis. Must be of class \code{geem}. Defaults to NULL.
#' @param mod.0 The model under the null hypothesis. Must be of class \code{geem}. Defaults to NULL.
#' @param bias.correct Boolean specifying whether a small-sample bias correction should be applied to the estimated sandwich variance-covariance matrix. Useful when the number of subjects is small. Defaults to FALSE.
#' @param correction.method A string specifying the correction method to be used. Currently supported options are "df" and "kc". Defaults to "kc".
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
                        bias.correct = FALSE,
                        correction.method = "kc",
                        id.vec = NULL,
                        verbose = FALSE) {
  # check inputs
  if (inherits(mod.1, "try-error") || inherits(mod.0, "try-error")) {
    res <- list(Wald_Stat = NA_real_,
                DF = NA_real_,
                P_Val = NA_real_,
                Notes = "No test performed due to model failure.")
    return(res)
  }
  correction.method <- tolower(correction.method)
  if (!correction.method %in% c("df", "kc")) { stop("Unsupported bias correction method in waldTestGEE().") }
  if (bias.correct && correction.method == "kc" && is.null(id.vec)) { stop("The Kauermann and Carroll bias correction method requires a vector of subject IDs.") }
  mod.1 <- mod.1$final_mod
  if (is.null(mod.1) || is.null(mod.0) || !(inherits(mod.1, "geem") && inherits(mod.0, "geem"))) { stop("You must provide two geeM models to waldTestGee().") }
  if (length(coef(mod.0)) != 1) { stop("Null GEE model must be intercept-only.") }
  if (length(coef(mod.1)) <= length(coef(mod.0))) {
    # can't calculate Wald statistic if both models are intercept-only
    res <- list(Wald_Stat = 0,
                DF = 0,
                P_Val = 1,
                Notes = NA_character_)
  } else {
    # compute test statistic & asymptotic p-value
    coef_alt_mod <- names(coef(mod.1))
    coef_null_mod <- c("Intercept")
    coef_diff <- setdiff(coef_alt_mod, coef_null_mod)
    coef_idx <- rep(0, length(coef_diff))
    for (i in seq_len(length(coef_diff))) {
      coef_idx[i] <- which(coef_diff[i] == coef_alt_mod)
    }
    coef_vals <- as.matrix(coef(mod.1)[coef_idx])
    vcov_mat <- as.matrix(mod.1$var)
    if (bias.correct) {
      vcov_mat <- biasCorrectGEE(mod.1,
                                 correction.method = correction.method,
                                 id.vec = id.vec,
                                 cor.structure = mod.1$corr,
                                 verbose = verbose)
    }
    vcov_mat <- vcov_mat[coef_idx, coef_idx]
    wald_test_stat <- try({
      as.numeric(crossprod(coef_vals, MASS::ginv(vcov_mat)) %*% coef_vals)
    }, silent = TRUE)
    if (inherits(wald_test_stat, "try-error")) {
      wald_note <- wald_test_stat[1]  # this is the error message
      wald_test_stat <- 0
      p_value <- 1
    } else {
      p_value <- as.numeric(1 - stats::pchisq(wald_test_stat, df = length(coef_diff)))
      wald_note <- NA_character_
    }
    res <- list(Wald_Stat = wald_test_stat,
                DF = length(coef_diff),
                P_Val = p_value,
                Notes = wald_note)
  }
  return(res)
}
