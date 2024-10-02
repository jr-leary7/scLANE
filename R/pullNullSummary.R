#' Generate a summary of the null model.
#'
#' @name pullNullSummary
#' @author Jack R. Leary
#' @author Rhonda Bacher
#' @description This function takes in the null model fitted during the running of \code{\link{marge2}} and summarizes it.
#' @import magrittr
#' @importFrom stats predict logLik deviance
#' @importFrom broom.mixed tidy
#' @importFrom dplyr rename
#' @param null.model The null model from \code{\link{marge2}}. Defaults to NULL.
#' @param is.gee Boolean specifying whether GEE mode was used in fitting the null model. Defaults to FALSE.
#' @param sandwich.var Boolean specifying whether the robust sandwich variance-covariance matrix should be used. Defaults to FALSE.
#' @param is.glmm Boolean specifying whether the GLMM mode was used in fitting the model. Defaults to FALSE.
#' @return A list containing a coefficient summary, fitted values and their standard errors, and the log-likelihood and deviance of the model.
#' @seealso \code{\link{marge2}}

pullNullSummary <- function(null.model = NULL,
                            is.gee = FALSE,
                            sandwich.var = FALSE,
                            is.glmm = FALSE) {
  # check inputs
  if (is.null(null.model)) { stop("A null model must be provided to pullNullSummary") }
  # handle the degenerate case
  if (inherits(null.model, "try-error")) {
    res <- list(null_sumy_df = NULL,
                null_pred_df = NULL,
                null_ll = NA_real_,
                null_dev = NA_real_,
                null_fit_notes = null.model[1])
  } else {
    # pull model summary
    if (is.gee) {
      null_sumy_df <- try({
        null_gee_summary <- summary(null.model)
        data.frame(term = null_gee_summary$coefnames,
                   estimate = unname(null_gee_summary$beta),
                   `std.error` = ifelse(sandwich.var, unname(null_gee_summary$se.robust), unname(null_gee_summary$se.model)),
                   statistic = unname(null_gee_summary$wald.test),
                   `p.value` = unname(null_gee_summary$p))
      }, silent = TRUE)
      null_pred_df <- try({
        if (sandwich.var) {
          vcov_mat <- as.matrix(null.model$var)
        } else {
          vcov_mat <- as.matrix(null.model$naiv.var)
        }
        data.frame(null_link_fit = stats::predict(null.model),
                   null_link_se = as.numeric(sqrt(tcrossprod(null.model$X, vcov_mat))))
      }, silent = TRUE)
    } else if (is.glmm) {
      null_sumy_df <- try({
        null_glmm_summary <- as.data.frame(broom.mixed::tidy(null.model, effects = "fixed"))
        data.frame(term = null_glmm_summary$term,
                   estimate = null_glmm_summary$estimate[1],
                   `std.error` = null_glmm_summary$std.error[1],
                   statistic = null_glmm_summary$statistic[1],
                   `p.value` = null_glmm_summary$p.value[1])
      }, silent = TRUE)
      null_pred_df <- try({
        data.frame(stats::predict(null.model, type = "link", se.fit = TRUE)[seq(2)]) %>%
          dplyr::rename(null_link_fit = fit, null_link_se = se.fit)
      }, silent = TRUE)
    }  else {
      null_sumy_df <- try({
        as.data.frame(broom.mixed::tidy(null.model)) # saves a few bytes by converting from tibble
      }, silent = TRUE)
      null_pred_df <- try({
        data.frame(stats::predict(null.model, type = "link", se.fit = TRUE)[seq(2)]) %>%
          dplyr::rename(null_link_fit = fit, null_link_se = se.fit)
      }, silent = TRUE)
    }
    res <- list(null_sumy_df = null_sumy_df,
                null_pred_df = null_pred_df,
                null_ll = ifelse(is.gee, NA_real_, as.numeric(stats::logLik(null.model))),
                null_dev = ifelse(is.gee, NA_real_, stats::deviance(null.model)))
  }
  return(res)
}
