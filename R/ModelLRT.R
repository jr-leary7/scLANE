#' Perform a likelihood ratio test for one model against another.
#'
#' @name modelLRT
#' @author Jack Leary
#' @description This function compares two models using a likelihood ratio test (LRT) under the assumption that the test statistic is asymptotically Chi-squared with degrees freedom equal to the difference in the number of parameters between the larger and smaller model.
#' @importFrom stats logLik pchisq
#' @param mod.1 The model corresponding to the alternative hypothesis. Defaults to NULL.
#' @param mod.0 The model corresponding to the null hypothesis. Defaults to NULL.
#' @param is.glmm Are the models being compared GLMMs? Defaults to FALSE.
#' @return A list containing the LRT test statistic, degrees freedom, and the \emph{p}-value computed using the Chi-squared assumption.
#' @export
#' @examples
#' \dontrun{
#' modelLRT(mod.1 = marge_mod, mod.0 = null_mod)
#' modelLRT(mod.1 = marge_mod,
#'          mod.0 = null_mod,
#'          is.glmm = TRUE)
#' }

modelLRT <- function(mod.1 = NULL,
                     mod.0 = NULL,
                     is.glmm = FALSE) {
  # check inputs
  if (is.null(mod.1) || is.null(mod.0)) { stop("You must provide two models to compare.") }
  LRT_notes <- NULL
  # compute LRT
  if (is.glmm) {
    if (!(inherits(mod.1, "glmmTMB") && inherits(mod.0, "glmmTMB"))) { stop("Models must be of class glmmTMB.") }
    mod.1_ll <- -mod.1$fit$objective
    mod.0_ll <- -mod.0$fit$objective
    lrt_stat <- 2 * (mod.1_ll - mod.0_ll)
    if (is.na(logLik(mod.1))) {
      if (!mod.1$sdr$pdHess) {
        LRT_notes <- "Non-positive definite Hessian in alternate GLMM, probably because of shallow log-likelihood."
      } else {
        LRT_notes <- "Log-likelihood for alternate GLMM is very shallow / flat, exercise caution."
      }
    }
  } else {
    if (!(inherits(mod.1, "glm") || inherits(mod.1, "lm"))) { stop("Models must be of class lm or glm.") }
    mod.1_ll <- stats::logLik(mod.1)
    mod.0_ll <- stats::logLik(mod.0)
    lrt_stat <- as.numeric(2 * (mod.1_ll - mod.0_ll))
  }
  dgr_free <- attributes(stats::logLik(mod.1))$df - attributes(stats::logLik(mod.0))$df
  if (dgr_free == 0) {
    dgr_free <- 1
  }
  p_val <- stats::pchisq(lrt_stat,
                         dgr_free,
                         lower = FALSE)
  res <- list(Alt_Mod_LL = as.numeric(mod.1_ll),
              Null_Mod_LL = as.numeric(mod.0_ll),
              LRT_Stat = lrt_stat,
              DF = dgr_free,
              P_Val = p_val,
              Model_Type = ifelse(is.glmm, "GLMM", "GLM"),
              Notes = LRT_notes)
  return(res)
}
