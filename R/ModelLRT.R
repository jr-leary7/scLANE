#' Perform a likelihood ratio test for one model against another.
#'
#' @name modelLRT
#' @description This function compares two models using a likelihood ratio test (LRT) under the assumption that the test statistic is asymptotically Chi-squared with degrees freedom equal to the difference in the number of parameters between the larger and smaller model.
#' @importFrom stats logLik
#' @param mod.1 The model corresponding to the alternative hypothesis. Defaults to NULL.
#' @param mod.0 The model corresponding to the null hypothesis. Defaults to NULL.
#' @return A list containing the LRT test statistic, degrees freedom, and the \emph{p}-value computed using the Chi-squared assumption.
#' @export
#' @examples
#' ModelLRT(mod1 = saturated_glm, mod0 = null_glm)

modelLRT <- function(mod.1 = NULL, mod.0 = NULL) {
  # check inputs
  if (is.null(mod.1) | is.null(mod.0)) stop("You must provide two models to compare.")
  if (!(any(class(mod.1) == "glm") | any(class(mod.1) == "lm"))) stop("Models must be of class lm or glm.")
  # compute LRT
  mod.1_ll <- logLik(mod.1)
  mod.0_ll <- logLik(mod.0)
  lrt_stat <- as.numeric(2 * (mod.1_ll - mod.0_ll))
  dgr_free <- length(coef(mod.1)) - length(coef(mod.0))
  if (dgr_free == 0) { dgr_free <- 1 }
  p_val <- pchisq(lrt_stat, dgr_free, lower = FALSE)
  res <- list(Alt_Mod_LL = mod.1_ll, Null_Mod_LL = mod.0_ll, LRT_Stat = lrt_stat, DF = dgr_free, P_Val = p_val)
  return(res)
}
