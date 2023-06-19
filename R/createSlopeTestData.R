#' A helper function to create a dataframe of breakpoints and associated \emph{p}-values from a \code{marge} model.
#'
#' @name createSlopeTestData
#' @author Jack Leary
#' @description Creates a data.frame of \code{marge} model breakpoints, \emph{p}-values, and other info.
#' @importFrom purrr map_dbl
#' @param marge.model A \code{marge} model object, like those returned from \code{\link{marge2}}. Defaults to NULL.
#' @param pt A data.frame containing pseudotime or latent time values. Defaults to NULL.
#' @param is.gee Was the GEE framework used? Defaults to FALSE.
#' @param is.glmm Was the GLMM framework used? Defaults to FALSE.
#' @return A data.frame containing model data.
#' @seealso \code{\link{marge2}}
#' @seealso \code{\link{testSlope}}
#' @examples
#' \dontrun{
#' createSlopeTestData(marge_mod, pt_df)
#' createSlopeTestData(marge_mod,
#'                     pt = pt_df,
#'                     is.glmm = TRUE)
#' }

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
    p_vals <- NA_real_
    mod_notes <- "MARGE model error"
  } else {
    if (!is.glmm && length(marge.model$marge_coef_names) == 1) {
      rounded_brkpts <- NA_real_
      brkpts <- NA_real_
      brkpt_dirs <- NA_character_
      p_vals <- NA_real_
      mod_notes <- "No non-intercept coefficients"
    } else {
      # grab the k breakpoints from the MARGE model
      model_breakpoints_rounded <- extractBreakpoints(marge.model, directions = TRUE)
      rounded_brkpts <- model_breakpoints_rounded$Breakpoint
      brkpt_dirs <- model_breakpoints_rounded$Direction
      brkpts <- purrr::map_dbl(rounded_brkpts, \(x) pt[, 1][which.min(abs(pt[, 1] - x))])
      # extract p-values for coefficients other than intercept
      if (is.gee) {
        p_vals <- summary(marge.model$final_mod)$p[-1]
      } else if (is.glmm) {
        p_vals <- unname(summary(marge.model$final_mod)$coefficients$cond[, "Pr(>|z|)"][-1])
      } else {
        p_vals <- unname(summary(marge.model$final_mod)$coefficients[, "Pr(>|z|)"][-1])
      }
      mod_notes <- rep(NA_character_, length(p_vals))
    }
  }
  # create table of results
  slope_df <- data.frame(Breakpoint = brkpts,
                         Rounded_Breakpoint = rounded_brkpts,
                         Direction = brkpt_dirs,
                         P_Val = p_vals,
                         Notes = mod_notes)
  return(slope_df)
}
