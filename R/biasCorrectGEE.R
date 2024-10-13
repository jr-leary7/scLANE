#' Bias-correct the GEE sandwich variance-covariance matrix.
#'
#' @name biasCorrectGEE
#' @author Jack R. Leary
#' @description This functions implements several bias-correction methods for the GEE sandwich variance-covariance matrix; they are to be used when the number of subjects is small or the numer of timepoints per-subject is very large.
#' @importFrom stats fitted.values cor
#' @importFrom dplyr with_groups summarise
#' @importFrom Matrix bdiag
#' @param fitted.model The fitted model of class \code{geem} returned by \code{\link{marge2}}. Defaults to NULL.
#' @param correction.method A string specifying the correction method to be used. Currently supported options are "df" and "kc". Defaults to "kc".
#' @param id.vec A vector of subject IDs. Defaults to NULL.
#' @param cor.structure A string specifying the correlation structure used in fitting the model. Defaults to "ar1".
#' @param verbose (Optional) A Boolean specifying whether or not verbose output should be printed to the console. Occasionally useful for debugging. Defaults to FALSE.
#' @return An object of class \code{matrix} containing the bias-corrected variance-covariance estimates.
#' @seealso \code{\link{waldTestGEE}}

biasCorrectGEE <- function(fitted.model = NULL,
                           correction.method = "kc",
                           id.vec = NULL,
                           cor.structure = "ar1",
                           verbose = FALSE) {
  # check inputs
  if (is.null(fitted.model)) { stop("Arguments to biasCorrectGEE() are missing.") }
  if (!inherits(fitted.model, "geem")) { stop("The fitted model must be of class 'geem'.") }
  correction.method <- tolower(correction.method)
  if (!correction.method %in% c("df", "kc")) { stop("Unsupported bias correction method in waldTestGEE().") }
  if (correction.method == "kc" && is.null(id.vec)) { stop("Kauermann and Carroll requires the provision of a vector of subject IDs.") }
  cor.structure <- tolower(cor.structure)
  if (!cor.structure %in% c("exchangeable", "ar1")) { stop("Unrecognized correlation structure in biasCorrectGEE().") }
  # extract sandwich variance-covariance matrix
  V_sandwich <- as.matrix(fitted.model$var)
  n_s <- length(fitted.model$clusz)
  if (correction.method == "df") {
    # compute degrees-of-freedom bias-corrected variance-covariance matrix
    p <- ncol(fitted.model$X)
    if (p >= n_s) {
      warning("Cannot perform DF bias correction on sandwich variance-covariance matrix as the number of subjects is less than or equal to the number of covariates.")
      res <- V_sandwich
    } else {
      df_factor <- n_s / (n_s - p)
      if (verbose) {
        message(paste0("DF factor estimated at: ", round(df_factor, 5)))
      }
      res <- df_factor * V_sandwich
    }
  } else if (correction.method == "kc") {
    # compute Kauermann and Carroll bias-corrected variance-covariance matrix
    resid_df <- data.frame(subject = id.vec, r_hat = fitted.model$y - stats::fitted.values(fitted.model))
    if (cor.structure == "exchangeable") {
      rhoExchangeable <- function(residuals = NULL) {
        n_i <- length(residuals)
        if (n_i < 2) {
          res <- NA_real_
        } else {
          sum_r <- sum(residuals)
          sum_r2 <- sum(residuals^2)
          denominator <- n_i * (sum_r2 - (sum_r^2 / n_i))
          if (denominator == 0) {
            res <- NA_real_
          } else {
            res <- (sum_r^2 - sum_r2) / denominator
          }
        }
        return(res)
      }
      rho_estimates <- dplyr::with_groups(resid_df,
                                          subject,
                                          dplyr::summarise,
                                          rho = rhoExchangeable(r_hat))
      rho_avg <- mean(rho_estimates$rho, na.rm = TRUE)
      if (verbose) {
        message(paste0("Exchangeable correlation parameter estimated at: ", round(rho_avg, 5)))
      }
    } else if (cor.structure == "ar1") {
      rhoAR1 <- function(residuals = NULL) {
        n_i <- length(residuals)
        if (n_i < 2) {
          res <- NA_real_
        } else {
          res <- stats::cor(residuals[-n_i],
                            residuals[-1],
                            use = "complete.obs")
        }
        return(res)
      }
      rho_estimates <- dplyr::with_groups(resid_df,
                                          subject,
                                          dplyr::summarise,
                                          rho = rhoAR1(r_hat))
      rho_avg <- mean(rho_estimates$rho, na.rm = TRUE)
      if (verbose) {
        message(paste0("AR-1 correlation parameter estimated at: ", round(rho_avg, 5)))
      }
    }
    createExchCor <- function(n = NULL, rho = NULL) {
      if (n == 1) {
        res <- matrix(1, nrow = 1, ncol = 1)
      } else {
        res <- matrix(rho, nrow = n, ncol = n) + diag(1 - rho, n, n)
      }
      return(res)
    }
    createAR1Cor <- function(n = NULL, rho = NULL) {
      if (n == 1) {
        res <- matrix(1, nrow = 1, ncol = 1)
      } else {
        res <- rho^abs(outer(seq(n), seq(n), "-"))
      }
      return(res)
    }
    cov_matrices <- vector("list", length = n_s)
    subjects <- unique(resid_df$subject)
    for (s in seq(subjects)) {
      subject_resid <- resid_df$r_hat[resid_df$subject == subjects[s]]
      n_i <- length(subject_resid)
      if (cor.structure == "exchangeable") {
        R_i <- createExchCor(n_i, rho = rho_avg)
      } else if (cor.structure == "ar1") {
        R_i <- createAR1Cor(n_i, rho = rho_avg)
      }
      sigma2 <- fitted.model$phi
      Var_Yi <- sigma2 * R_i
      Var_Yi_inv <- try({ eigenMapMatrixInvert(Var_Yi, n_cores = 1L) }, silent = TRUE)
      if (inherits(Var_Yi_inv, "try-error")) {
        Var_Yi_inv <- try({ eigenMapPseudoInverse(Var_Yi, n_cores = 1L) }, silent = TRUE)
      }
      if (inherits(Var_Yi_inv, "try-error")) {
        if (verbose) {
          warning(paste0("Covariance matrix inversion failed for subject ",
                         subjects[s],
                         ", applying diagonal matrix as fallback in biasCorrectGEE()."))
        }
        Var_Yi_inv <- diag(1, n_i, n_i)
      }
      cov_matrices[[s]] <- Var_Yi_inv
    }
    W <- as.matrix(Matrix::bdiag(cov_matrices))
    X <- fitted.model$X
    XWX <- t(X) %*% W %*% X
    XWX_inv <- try({ eigenMapMatrixInvert(XWX, n_cores = 1L) }, silent = TRUE)
    if (inherits(XWX_inv, "try-error")) {
      XWX_inv <- eigenMapPseudoInverse(XWX, n_cores = 1L)
    }
    H <- X %*% XWX_inv %*% t(X) %*% W
    tr_H <- sum(diag(H))
    if (verbose) {
      message(paste0("Trace of projection matrix H estimated at: ", round(tr_H, 5)))
    }
    n <- length(fitted.model$y)
    kc_factor <- (n_s / (n_s - 1)) / (1 - (tr_H / n))
    if (verbose) {
      message(paste0("KC factor estimated at: ", round(kc_factor, 5)))
    }
    res <- kc_factor * V_sandwich
  }
  return(res)
}
