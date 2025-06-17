#' Calculate part of the score statistic for a GLM.
#'
#' @name stat_out_score_glm_null
#' @author Jakub Stoklosa
#' @author David I. Warton
#' @author Jack R. Leary
#' @importFrom gamlss gamlss
#' @importFrom stats fitted.values
#' @description A function that calculates parts of the score statistic for 
#' GLMs only (it is used for the full path for forward selection).
#' @param Y : The response variable Defaults to NULL.
#' @param B_null : Design matrix under the null model. Defaults to NULL.
#' @return A list of values (mainly products of matrices) that make up the 
#' final score statistic calculation (required for another function).
#' @references Stoklosa, J., Gibb, H. and Warton, D.I. (2014). Fast forward 
#' selection for generalized estimating equations with a large number of 
#' predictor variables. \emph{Biometrics}, \strong{70}, 110--120.
#' @references Stoklosa, J. and Warton, D.I. (2018). A generalized estimating 
#' equation approach to multivariate adaptive regression splines. 
#' \emph{Journal of Computational and Graphical Statistics}, 
#' \strong{27}, 245--253.
#' @seealso \code{\link{stat_out}}
#' @seealso \code{\link{stat_out_score_gee_null}}

stat_out_score_glm_null <- function(Y = NULL, B_null = NULL) {
    # check inputs
    if (is.null(Y) || is.null(B_null)) {
        stop("Some inputs to stat_out_score_glm_null() are missing.")
    }
    # extract statistics
    ests <- gamlss::gamlss(Y ~ B_null - 1,
        family = "NBI",
        trace = FALSE,
        data = NULL
    )
    sigma_est <- exp(ests$sigma.coef)
    mu_est <- as.matrix(stats::fitted.values(ests))
    # Type I NB variance = mu (1 + mu * sigma); sigma = 1 / theta
    V_est <- mu_est * (1 + mu_est * sigma_est)
    VS_est_list <- (Y - mu_est) / V_est
    mu_V_diag <- diag(c(mu_est^2 / V_est))
    temp_prod <- eigenMapMatMult(A = t(B_null), B = mu_V_diag)
    A_list_inv <- eigenMapMatMult(A = temp_prod, B = B_null)
    if (ncol(A_list_inv) == 1 && nrow(A_list_inv) == 1) {
        A_list <- 1 / A_list_inv
    } else {
        A_list <- try(
            {
                eigenMapMatrixInvert(A_list_inv)
            },
            silent = TRUE
        )
        if (inherits(A_list, "try-error")) {
            A_list <- eigenMapPseudoInverse(A_list_inv)
        }
    }
    B1_list <- eigenMapMatMult(A = t(B_null), B = mu_V_diag)
    res <- list(
        VS.est_list = VS_est_list,
        A_list = A_list,
        B1_list = B1_list,
        mu.est = mu_est,
        V.est = V_est
    )
    return(res)
}
