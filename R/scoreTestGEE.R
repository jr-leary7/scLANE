#' Use a Lagrange multiplier (score) test to compare nested GEE models.
#'
#' @name scoreTestGEE
#' @author Jack R. Leary
#' @description Performs a basic Lagrange multiplier test to determine whether an alternate model is significantly better than a nested null model. This is the GEE equivalent (kind of) of \code{\link{modelLRT}}. Be careful with small sample sizes.
#' @importFrom stats model.matrix predict pchisq
#' @importFrom Matrix bdiag
#' @param mod.1 The model under the alternative hypothesis. Must be of class \code{geem}. Defaults to NULL.
#' @param mod.0 The model under the null hypothesis. Must be of class \code{geem}. Defaults to NULL.
#' @param alt.df The dataframe used to fit the alternative model. Defaults to NULL.
#' @param null.df The dataframe used to fit the null model. Defaults to NULL.
#' @param id.vec A vector of subject IDs to use as input to \code{\link{marge2}}. Defaults to NULL.
#' @param cor.structure A string specifying the working correlation structure used to fit each model. Must be one of "ar1", "independence", or "exchangeable". Defaults to "ar1".
#' @return A list containing the Score test statistic, a \emph{p}-value, and the degrees of freedom used in the test.
#' @details
#' \itemize{
#' \item Calculating the test statistic involves taking the inverse of the variance of the score vector. Ideally this would be done using the true inverse, but in practice this can cause issues when the matrix is near-singular. With this in mind, we use the Moore-Penrose pseudoinverse if the original matrix inversion fails.
#' \item The \emph{p}-value is calculated using an asymptotic \eqn{\Chi^2} distribution, with the degrees of freedom equal to the number of non-intercept coefficients in the alternative model.
#' }
#' @seealso \code{\link[geeM]{geem}}
#' @seealso \code{\link[glmtoolbox]{anova2}}
#' @seealso \code{\link{waldTestGEE}}
#' @seealso \code{\link{modelLRT}}

scoreTestGEE <- function(mod.1 = NULL,
                         mod.0 = NULL,
                         alt.df = NULL,
                         null.df = NULL,
                         id.vec = NULL,
                         cor.structure = "ar1") {
  # check inputs
  if (is.null(mod.1) || is.null(mod.0) || is.null(alt.df) || is.null(null.df) || is.null(id.vec)) { stop("Please provide all inputs to scoreTestGEE().") }
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
    rho <- summary(mod.0)$alpha
    phi <- summary(mod.0)$phi
    theta <- as.numeric(gsub("\\)", "", gsub(".*\\(", "", mod.1$FunList$family)))
    X_null <- stats::model.matrix(mod.0, data = null.df)
    X_alt <- stats::model.matrix(mod.1, data = alt.df)
    r_null <- null.df$Y - exp(stats::predict(mod.0))
    p_alt <- ncol(X_alt)
    groups <- unique(id.vec)
    W <- K <- vector("list", length = length(groups))
    for (i in seq(groups)) {
      group_idx <- which(id.vec == groups[i])
      n_i <- length(group_idx)
      # create working correlation matrix R_i
      if (cor.structure == "independence") {
        R_i <- diag(1, nrow = n_i, ncol = n_i)
      } else if (cor.structure == "exchangeable") {
        R_i <- matrix(rho, nrow = n_i, ncol = n_i)
        diag(R_i) <- 1
      } else if (cor.structure == "ar1") {
        R_i <- matrix(rho^abs(outer(seq(n_i), seq(n_i), "-")), nrow = n_i, ncol = n_i)
      }
      # create working covariance matrix V_i
      mu_i <- exp(stats::predict(mod.0)[group_idx])
      V_mu_i <- mu_i * (1 + mu_i / theta)
      A_i <- diag(V_mu_i)
      A_i_sqrt <- sqrt(A_i)  # same as taking the 1/2 power of A_i since A_i is diagonal
      V_i <- A_i_sqrt %*% R_i %*% A_i_sqrt
      # create matrices W_i and K_i
      K_i <- diag(mu_i)
      V_i_inv <- try({ eigenMapMatrixInvert(V_i) }, silent = TRUE)
      if (inherits(V_i_inv, "try-error")) {
        V_i_inv <- eigenMapPseudoInverse(V_i)
      }
      W_i <- K_i %*% V_i_inv %*% K_i
      W[[i]] <- W_i
      K[[i]] <- K_i
    }
    W <- as.matrix(Matrix::bdiag(W))
    K <- as.matrix(Matrix::bdiag(K))
    K_inv <- try({ eigenMapMatrixInvert(K) }, silent = TRUE)
    if (inherits(K_inv, "try-error")) {
      K_inv <- eigenMapPseudoInverse(K)
    }
    # generate score vector
    U <- phi^(-1) * t(X_alt) %*% W %*% K_inv %*% r_null
    # generate variance of score vector and invert it
    V_U <- phi^(-2) * t(X_alt) %*% W %*% X_alt
    V_U_inv <- try({ eigenMapMatrixInvert(V_U) }, silent = TRUE)
    if (inherits(V_U_inv, "try-error")) {
      V_U_inv <- eigenMapPseudoInverse(V_U)
    }
    # estimate test statistic and accompanying p-value
    S <- t(U) %*% V_U_inv %*% U
    S_df <- ncol(X_alt) - ncol(X_null)
    p_value <- 1 - stats::pchisq(S, df = S_df)
    res <- list(Score_Stat = S,
                DF = S_df,
                P_Val = p_value,
                Notes = NA_character_)
  }
  return(res)
}
