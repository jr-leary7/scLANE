#' Perform a likelihood ratio test for one model against another.
#'
#' @name modelLRT
#' @author Jack Leary
#' @description This function compares two models using a likelihood ratio test (LRT) under the assumption that the test statistic is asymptotically Chi-squared with degrees freedom equal to the difference in the number of parameters between the larger and smaller model.
#' @importFrom stats logLik pchisq
#' @param mod.1 The model corresponding to the alternative hypothesis. Defaults to NULL.
#' @param mod.0 The model corresponding to the null hypothesis. Defaults to NULL.
#' @return A list containing the LRT test statistic, degrees freedom, and the \emph{p}-value computed using the Chi-squared assumption.
#' @export
#' @examples
#' \dontrun{ModelLRT(mod1 = marge_mod, mod0 = null_glm)}

modelLRT <- function(mod.1 = NULL, mod.0 = NULL) {
  # check inputs
  if (is.null(mod.1) | is.null(mod.0)) { stop("You must provide two models to compare.") }
  if (!(any(class(mod.1) == "glm") | any(class(mod.1) == "lm"))) { stop("Models must be of class lm or glm.") }
  # compute LRT
  mod.1_ll <- stats::logLik(mod.1)
  mod.0_ll <- stats::logLik(mod.0)
  lrt_stat <- as.numeric(2 * (mod.1_ll - mod.0_ll))
  dgr_free <- attributes(mod.1_ll)$df - attributes(mod.0_ll)$df
  if (dgr_free == 0) { dgr_free <- 1 }
  p_val <- stats::pchisq(lrt_stat, dgr_free, lower = FALSE)
  res <- list(Alt_Mod_LL = as.numeric(mod.1_ll),
              Null_Mod_LL = as.numeric(mod.0_ll),
              LRT_Stat = lrt_stat,
              DF = dgr_free,
              P_Val = p_val)
  return(res)
}
