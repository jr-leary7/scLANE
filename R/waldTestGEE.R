#' Use a Wald test to compare nested GEE models.
#'
#' @name waldTestGEE
#' @author Jack Leary
#' @description Performs a basic Wald test to determine whether an alternate model is significantly better than a nested null model. This is the GEE equivalent (kind of) of \code{\link{modelLRT}}. Be careful with small sample sizes.
#' @importFrom MASS ginv
#' @importFrom stats pchisq
#' @param mod.1 The alternate model. Must be of class \code{geem}. Defaults to NULL.
#' @param mod.0 The null model. Must be of class \code{geem}. Defaults to NULL.
#' @return A list containing the Wald test statistic, a \emph{p}-value, and the degrees of freedom used in the test.
#' @seealso \code{\link[multgee]{waldts}}
#' @seealso \code{\link[geeM]{geem}}
#' @seealso \code{\link{modelLRT}}
#' @export
#' @examples
#' \dontrun{
#' waldTestGEE(mod.1 = full_model, mod.0 = null_model)
#' }

waldTestGEE <- function(mod.1 = NULL, mod.0 = NULL) {
  # check inputs
  if (is.null(mod.1) | is.null(mod.0) | (class(mod.1) != "geem" | class(mod.0) != "geem")) { stop("You must provide two geeM models to wald_test_gee().") }
  if (length(coef(mod.1)) <= length(coef(mod.0))) {
    # can't calculate Wald statistic if both models are intercept-only
    res <- list(Wald_Stat = 0,
                DF = 0,
                P_Val = 1,
                Notes = NA_character_)
  } else {
    # generate wald statistic & p-value -- see multgee::waldts() at https://mran.microsoft.com/snapshot/2016-01-07/web/packages/multgee/vignettes/multgee.pdf
    coef_alt_mod <- names(coef(mod.1))
    coef_alt_mod[coef_alt_mod == "B_finalIntercept"] <- "(Intercept)"  # correct name of intercept to correctly compare coefficient sets
    coef_null_mod <- names(coef(mod.0))
    coef_diff <- setdiff(coef_alt_mod, coef_null_mod)
    coef_idx <- rep(0, length(coef_diff))
    for (i in seq_len(length(coef_diff))) {
      coef_idx[i] <- which(coef_diff[i] == coef_alt_mod)
    }
    coef_vals <- as.matrix(coef(mod.1)[coef_idx])
    vcov_mat <- as.matrix(mod.1$var)[coef_idx, coef_idx]
    wald_test <- try({
      as.numeric(t(coef_vals) %*% solve(vcov_mat) %*% coef_vals)
    }, silent = TRUE)
    if (all(class(wald_test) == "try-error")) {
      wald_test <- try({
        as.numeric(t(coef_vals) %*% MASS::ginv(vcov_mat) %*% coef_vals)  # use Moore-Penrose pseudoinverse if vcov matrix is singular
      }, silent = TRUE)
      if (all(class(wald_test) == "try-error")) {
        wald_res <- 0
        p_value <- 1
        wald_note <- wald_test[1]  # this is the error message
      } else {
        wald_res <- wald_test
        p_value <- as.numeric(1 - stats::pchisq(wald_test, df = length(coef_diff)))
        wald_note <- "vcov matrix was singular, used Moore-Penrose pseudoinverse instead"
      }
    } else {
      wald_res <- wald_test
      p_value <- as.numeric(1 - stats::pchisq(wald_test, df = length(coef_diff)))
      wald_note <- NA_character_
    }
    res <- list(Wald_Stat = wald_res,
                DF = length(coef_diff),
                P_Val = p_value,
                Notes = wald_note)
  }
  return(res)
}
