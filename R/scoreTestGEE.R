#' Use a Lagrange multiplier (score) test to compare nested GEE models.
#'
#' @name scoreTestGEE
#' @author Jack R. Leary
#' @description Performs a basic Lagrange multiplier test to determine whether an alternate model is significantly better than a nested null model. This is the GEE equivalent (kind of) of \code{\link{modelLRT}}. Be careful with small sample sizes.
#' @importFrom stats model.matrix predict pchisq
#' @param mod.1 The model under the alternative hypothesis. Must be of class \code{geem}. Defaults to NULL.
#' @param mod.0 The model under the null hypothesis. Must be of class \code{geem}. Defaults to NULL.
#' @param alt.df The dataframe used to fit the alternative model. Defaults to NULL. 
#' @param null.df The dataframe used to fit the null model. Defaults to NULL. 
#' @param sandwich.var Boolean specifying whether the sandwich variance-covariance matrix should be used. Defaults to FALSE.
#' @return A list containing the Score test statistic, a \emph{p}-value, and the degrees of freedom used in the test.
#' @details
#' \itemize{
#' \item Calculating the test statistic involves taking the inverse of the variance-covariance matrix of the coefficients. Ideally this would be done using the "true" inverse with something like \code{\link{solve}}, \code{\link{qr.solve}}, or \code{\link{chol2inv}}, but in practice this can cause issues when the variance-covariance matrix is near-singular. With this in mind, we use the Moore-Penrose pseudoinverse as implemented in \code{\link[MASS]{ginv}} instead, which leads to more stable results.
#' \item The \emph{p}-value is calculated using an asymptotic \eqn{\Chi^2} distribution, with the degrees of freedom equal to the number of non-intercept coefficients in the alternative model.
#' }
#' @seealso \code{\link[geeM]{geem}}
#' @seealso \code{\link{waldTestGEE}}
#' @seealso \code{\link{modelLRT}}

scoreTestGEE <- function(mod.1 = NULL, 
                         mod.0 = NULL, 
                         alt.df = NULL, 
                         null.df = NULL, 
                         sandwich.var = FALSE) {
  # check inputs 
  if (is.null(mod.1) || is.null(mod.0) || is.null(alt.df) || is.null(null.df)) { stop("Please provide all inputs to scoreTestGEE().") }
  if (inherits(mod.1, "try-error") || inherits(mod.0, "try-error")) {
    res <- list(Score_Stat = NA_real_,
                DF = NA_real_,
                P_Val = NA_real_,
                Notes = "No test performed due to model failure.")
    return(res)
  }
  mod.1 <- mod.1$final_mod
  if (!(inherits(mod.1, "geem") && inherits(mod.0, "geem"))) { stop("You must provide two geeM models to scoreTestGee().") }
  if (length(coef(mod.0)) != 1) { stop("Null GEE model must be intercept-only.") }
  if (length(coef(mod.1)) <= length(coef(mod.0))) {
    # can't calculate Score statistic if both models are intercept-only
    res <- list(Score_Stat = 0,
                DF = 0,
                P_Val = 1,
                Notes = NA_character_)
  } else {
    X_null <- stats::model.matrix(mod.0, data = null.df)
    X_alt <- stats::model.matrix(mod.1, data = alt.df)
    residuals_null <- null.df$Y - exp(stats::predict(mod.0))
    if (sandwich.var) {
      V_null <- as.matrix(mod.0$var)
    } else {
      V_null <- as.matrix(mod.0$naiv.var)
    }
    X_alt_t <- t(X_alt)
    U <- X_alt_t %*% residuals_null
    V_U <- (X_alt_t * as.numeric(V_null)) %*% X_alt
    V_U_inv <- try({ eigenMapMatrixInvert(V_U) }, silent = TRUE)
    if (inherits(V_U_inv, "try-error")) {
      V_U_inv <- eigenMapPseudoInverse(V_U)
    }
    test_stat <- as.numeric(t(U) %*% V_U_inv %*% U)
    df1 <- ncol(X_alt) - ncol(X_null)
    p_value <- 1 - stats::pchisq(test_stat, df = df1)
    res <- list(Score_Stat = test_stat,
                DF = df1,
                P_Val = p_value,
                Notes = NA_character_)
  }
  return(res)
}
