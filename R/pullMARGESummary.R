#' Generate a summary of the MARGE model.
#'
#' @name pullMARGESummary
#' @author Jack R. Leary
#' @author Rhonda Bacher
#' @description This function takes in the \code{MARGE} model fitted during the 
#' running of \code{\link{marge2}} and summarizes it.
#' @import magrittr
#' @importFrom stats predict logLik deviance
#' @importFrom broom.mixed tidy
#' @importFrom dplyr rename
#' @param marge.model The \code{MARGE} model from \code{\link{marge2}}.
#' Defaults to NULL.
#' @param is.gee Boolean specifying whether GEE mode was used in fitting the 
#' null model. Defaults to FALSE.
#' @param sandwich.var Boolean specifying whether the robust sandwich
#' variance-covariance matrix should be used. Defaults to FALSE.
#' @param is.glmm Boolean specifying whether the GLMM mode was used in fitting 
#' the model. Defaults to FALSE.
#' @return A list containing a coefficient summary, fitted values and their 
#' standard errors, and the log-likelihood and deviance of the model.
#' @seealso \code{\link{marge2}}

pullMARGESummary <- function(
        marge.model = NULL,
        is.gee = FALSE,
        sandwich.var = FALSE,
        is.glmm = FALSE) {
    # check inputs
    if (is.null(marge.model)) {
        stop("A null model must be provided to pullMARGESummary().")
    }
    # handle the degenerate case
    if (inherits(marge.model, "try-error")) {
        res <- list(
            marge_pred_df = NULL,
            marge_sumy_df = NULL,
            ll_marge = NA_real_,
            marge_fit_notes = marge.model[1]
        )
    } else {
        # pull model summary
        if (is.gee) {
            if (sandwich.var) {
                vcov_mat <- as.matrix(marge.model$final_mod$var)
            } else {
                vcov_mat <- as.matrix(marge.model$final_mod$naiv.var)
            }
            marge_pred_df <- try(
                {
                    data.frame(
                        marge_link_fit = predict(marge.model$final_mod),
                        marge_link_se = sqrt(apply(tcrossprod(
                            marge.model$final_mod$X,
                            vcov_mat
                        ) * marge.model$final_mod$X, 1, sum))
                    )
                },
                silent = TRUE
            )
            marge_sumy_df <- try(
                {
                    marge_gee_summary <- summary(marge.model$final_mod)
                    data.frame(
                        term = marge_gee_summary$coefnames,
                        estimate = unname(marge_gee_summary$beta),
                        `std.error` = ifelse(sandwich.var,
                            unname(marge_gee_summary$se.robust),
                            unname(marge_gee_summary$se.model)
                        ),
                        statistic = unname(marge_gee_summary$wald.test),
                        `p.value` = unname(marge_gee_summary$p)
                    )
                },
                silent = TRUE
            )
        } else if (is.glmm) {
            marge_pred_df <- try(
                {
                    data.frame(predict(marge.model$final_mod,
                        type = "link",
                        se.fit = TRUE
                    )[seq(2)]) %>%
                        dplyr::rename(
                            marge_link_fit = fit,
                            marge_link_se = se.fit
                        )
                },
                silent = TRUE
            )
            marge_sumy_df <- try({
                marge_glmm_summary <- broom.mixed::tidy(marge.model$final_mod,
                        effects = "fixed"
                    )
                    data.frame(
                        term = marge.model$marge_coef_names,
                        estimate = marge_glmm_summary$estimate,
                        `std.error` = marge_glmm_summary$std.error,
                        statistic = marge_glmm_summary$statistic,
                        `p.value` = marge_glmm_summary$p.value
                    )
                },
                silent = TRUE
            )
        } else {
            marge_pred_df <- try(
                {
                    data.frame(stats::predict(marge.model$final_mod,
                        type = "link",
                        se.fit = TRUE
                    )[seq(2)]) %>%
                        dplyr::rename(
                            marge_link_fit = fit,
                            marge_link_se = se.fit
                        )
                },
                silent = TRUE
            )
            marge_sumy_df <- try(
                {
                    as.data.frame(broom.mixed::tidy(marge.model$final_mod)) %>%
                        lapply(unname) %>%
                        as.data.frame()
                },
                silent = TRUE
            )
        }

        # get log-likelihood for GLMM or GLM cases
        if (is.glmm) {
            ll_marge <- -marge.model$final_mod$fit$objective
        } else if (is.gee) {
            ll_marge <- NA_real_
        } else {
            ll_marge <- as.numeric(stats::logLik(marge.model$final_mod))
        }
        # check positive-definiteness of Hessian for GLMM -- might have an
        # effect on LRT stat / accompanying p-value
        if (is.glmm) {
            if (!marge.model$final_mod$sdr$pdHess) {
                marge_fit_notes <- "Non-positive definite Hessian in GLMM,
                probably due to shallow log-likelihood. Be careful!"
            } else {
                marge_fit_notes <- NA_character_
            }
        } else {
            marge_fit_notes <- NA_character_
        }
        res <- list(
            marge_pred_df = marge_pred_df,
            marge_sumy_df = marge_sumy_df,
            ll_marge = ll_marge,
            marge_dev = ifelse(is.gee, NA_real_,
                stats::deviance(marge.model$final_mod)
            ),
            marge_fit_notes = marge_fit_notes
        )
    }
    return(res)
}
