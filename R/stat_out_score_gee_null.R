#' Calculate part of the score statistic for a GEE.
#'
#' @name stat_out_score_gee_null
#' @author Jakub Stoklosa
#' @author David I. Warton
#' @author Jack R. Leary
#' @importFrom geeM geem
#' @importFrom MASS negative.binomial
#' @importFrom stats fitted.values
#' @description A function that calculates parts of the score statistic for 
#' GEEs only (it is used for the full path for forward selection).
#' @param Y The response variable Defaults to NULL.
#' @param B_null The design matrix matrix under the null model Defaults to NULL.
#' @param id.vec A vector of ID values for the observations. Defaults to NULL.
#' @param cor.structure A string specifying the desired correlation structure 
#' for the NB GEE. Defaults to NULL.
#' @param theta.hat Estimated value to be treated as the "known" theta when 
#' passing the negative binomial family to \code{\link[geeM]{geeM}}. 
#' Defaults to NULL.
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
#' @seealso \code{\link{stat_out_score_glm_null}}

stat_out_score_gee_null <- function(Y = NULL,
    B_null = NULL,
    id.vec = NULL,
    cor.structure = NULL,
    theta.hat = NULL) {
    # check inputs
    if (is.null(Y) || is.null(B_null) || is.null(id.vec) ||
        is.null(cor.structure) || is.null(theta.hat)) {
        stop("All arguments to stat_out_score_gee_null() must be non-NULL.")
    }
    cor.structure <- tolower(cor.structure)
    if (!cor.structure %in% c("independence", "exchangeable", "ar1")) {
        stop("cor.structure in stat_out_score_gee_null() must be a known type.")
    }
    n_vec <- as.numeric(table(id.vec))
    N <- length(unique(id.vec))
    # fit null NB GEE
    ests <- geeM::geem(Y ~ B_null - 1,
        id = id.vec,
        data = NULL,
        corstr = cor.structure,
        family = MASS::negative.binomial(50, link = "log"),
        scale.fix = FALSE,
        sandwich = FALSE
    )
    alpha_est <- ests$alpha
    sigma_est <- ests$phi
    mu_est <- as.matrix(stats::fitted.values(ests))
    # Type I NB variance = mu (1 + mu * sigma); sigma = 1 / phi
    V_est <- mu_est * (1 + mu_est * sigma_est)
    n_vec1 <- c(0, n_vec)
    Sigma2_list <- vector("list", length = N)
    VS_est_list <- AWA_est_list <- J2_list <- Sigma2_list
    J11 <- Sigma11 <- matrix(0, ncol(B_null), ncol(B_null))

    for (i in seq(N)) {
        k <- sum(n_vec[seq(i)])
        # set up working correlation matrix structure
        if (cor.structure == "independence") {
            R_alpha <- diag(1, n_vec[i], n_vec[i])
        } else if (cor.structure == "exchangeable") {
            R_alpha <- matrix(alpha_est, nrow = n_vec[i], ncol = n_vec[i])
            R_alpha <- R_alpha + diag(1 - alpha_est,
                nrow = n_vec[i],
                ncol = n_vec[i]
            )
        } else if (cor.structure == "ar1") {
            R_alpha <- alpha_est^outer(
                seq(n_vec[i]), seq(n_vec[i]),
                \(x, y) abs(x - y)
            )
        } else {
            stop("Currently unsupported correlation structure.")
        }
        # compute iteration values for each statistic
        temp_seq_n <- (sum(n_vec1[seq(i)]) + 1):k
        mu_est_sub <- mu_est[temp_seq_n]
        diag_sqrt_V_est <- diag(sqrt(V_est[temp_seq_n]),
            nrow = n_vec[i],
            ncol = n_vec[i]
        )
        temp_prod <- eigenMapMatMult(A = diag_sqrt_V_est, B = R_alpha)
        V_est_i <- eigenMapMatMult(A = temp_prod, B = diag_sqrt_V_est)
        V_est_i_inv <- try(
            {
                eigenMapMatrixInvert(V_est_i)
            },
            silent = TRUE
        )
        if (inherits(V_est_i_inv, "try-error")) {
            V_est_i_inv <- eigenMapPseudoInverse(V_est_i)
        }
        S_est_i <- t(Y)[temp_seq_n] - mu_est_sub
        temp_prod <- eigenMapMatMult(A = S_est_i, B = t(S_est_i))
        temp_prod <- eigenMapMatMult(A = V_est_i_inv, B = temp_prod)
        AWA_est_i <- eigenMapMatMult(A = temp_prod, B = V_est_i_inv)
        D_est_i <- eigenMapMatMult(A = diag(mu_est_sub,
            nrow = n_vec[i],
            ncol = n_vec[i]
        ), B = B_null[temp_seq_n, ])
        D_est_i_transpose <- t(D_est_i)
        temp_prod <- eigenMapMatMult(A = D_est_i_transpose, B = V_est_i_inv)
        J1_i <- eigenMapMatMult(A = temp_prod, B = D_est_i)
        J11 <- J11 + J1_i
        J2_i <- eigenMapMatMult(A = D_est_i_transpose, B = V_est_i_inv)
        Sigma2_i <- eigenMapMatMult(A = D_est_i_transpose, B = AWA_est_i)
        Sigma1_i <- eigenMapMatMult(A = Sigma2_i, B = D_est_i)
        Sigma11 <- Sigma11 + Sigma1_i
        V_est_list_elem_2 <- eigenMapMatMult(A = V_est_i_inv, B = S_est_i)
        VS_est_list[[i]] <- V_est_list_elem_2
        AWA_est_list[[i]] <- AWA_est_i
        J2_list[[i]] <- J2_i
        Sigma2_list[[i]] <- Sigma2_i
    }
    if (nrow(J11) == 1 && ncol(J11) == 1) {
        J11_inv <- 1 / J11
    } else {
        J11_inv <- try(
            {
                eigenMapMatrixInvert(J11)
            },
            silent = TRUE
        )
        if (inherits(J11_inv, "try-error")) {
            J11_inv <- eigenMapPseudoInverse(J11)
        }
    }
    temp_prod <- eigenMapMatMult(A = J11_inv, B = Sigma11)
    JSigma11 <- eigenMapMatMult(A = temp_prod, B = J11_inv)
    res <- list(
        VS.est_list = VS_est_list,
        AWA.est_list = AWA_est_list,
        J2_list = J2_list,
        J11.inv = J11_inv,
        Sigma2_list = Sigma2_list,
        JSigma11 = JSigma11,
        mu.est = mu_est,
        V.est = V_est
    )
    return(res)
}
