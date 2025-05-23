#' A helper function to create a dataframe of breakpoints and associated \emph{p}-values from a \code{marge} model.
#'
#' @name createSlopeTestData
#' @author Jack R. Leary
#' @description Creates a data.frame of \code{marge} model breakpoints, \emph{p}-values, and other info.
#' @importFrom purrr map_dbl
#' @param marge.model A \code{marge} model object, like those returned from \code{\link{marge2}}. Defaults to NULL.
#' @param pt A data.frame containing pseudotime or latent time values. Defaults to NULL.
#' @param is.gee Was the GEE mode used? Defaults to FALSE.
#' @param is.glmm Was the GLMM mode used? Defaults to FALSE.
#' @return A data.frame containing model data.
#' @seealso \code{\link{marge2}}
#' @seealso \code{\link{testSlope}}

createSlopeTestData <- function(marge.model = NULL,
                                pt = NULL,
                                is.gee = FALSE,
                                is.glmm = FALSE) {
  # check inputs
  if (is.null(marge.model) || is.null(pt)) { stop("Input to createSlopeTestData() should be non-null.") }
  # run function
  # checks to see if marge model was set to "try-error" in testDynamic()
  if (inherits(marge.model, "try-error")) {
    rounded_brkpts <- NA_real_
    brkpts <- NA_real_
    brkpt_dirs <- NA_character_
    est_betas <- NA_real_
    p_vals <- NA_real_
    test_stats <- NA_real_
    mod_notes <- "MARGE model error"
  } else {
    # check to see if final marge model is intercept-only
    if ((!is.gee && !is.glmm && length(marge.model$final_mod$coefficients) == 1) || (is.gee && length(marge.model$final_mod$beta) == 1) || (is.glmm && nrow(summary(marge.model$final_mod)$coefficients$cond) == 1)) {
      rounded_brkpts <- NA_real_
      brkpts <- NA_real_
      brkpt_dirs <- NA_character_
      est_betas <- NA_real_
      test_stats <- NA_real_
      p_vals <- NA_real_
      mod_notes <- "No non-intercept coefficients"
    } else {
      # grab the k breakpoints from the MARGE model
      model_breakpoints_rounded <- extractBreakpoints(marge.model, directions = TRUE)
      rounded_brkpts <- model_breakpoints_rounded$Breakpoint
      brkpt_dirs <- model_breakpoints_rounded$Direction
      brkpts <- purrr::map_dbl(rounded_brkpts, \(x) pt[, 1][which.min(abs(pt[, 1] - x))])
      # extract p-values for coefficients other than intercept
      model_sumy <- summary(marge.model$final_mod)
      if (is.gee) {
        est_betas <- model_sumy$beta[-1]
        test_stats <- model_sumy$wald.test[-1]
        p_vals <- model_sumy$p[-1]
      } else if (is.glmm) {
        est_betas <- unname(model_sumy$coefficients$cond[, "Estimate"][-1])
        test_stats <- unname(model_sumy$coefficients$cond[, "z value"][-1])
        p_vals <- unname(model_sumy$coefficients$cond[, "Pr(>|z|)"][-1])
      } else {
        est_betas <- unname(model_sumy$coefficients[, "Estimate"][-1])
        test_stats <- unname(model_sumy$coefficients[, "z value"][-1])
        p_vals <- unname(model_sumy$coefficients[, "Pr(>|z|)"][-1])
      }
      mod_notes <- rep(NA_character_, length(p_vals))
    }
  }
  # create table of results
  slope_df <- data.frame(Breakpoint = brkpts,
                         Rounded_Breakpoint = rounded_brkpts,
                         Direction = brkpt_dirs,
                         Beta = est_betas,
                         Test_Stat = test_stats,
                         P_Val = p_vals,
                         Notes = mod_notes)
  return(slope_df)
}
