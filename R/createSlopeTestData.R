#' A helper function to create a dataframe of breakpoints and associated \emph{p}-values from a \code{marge} model.
#'
#' @name createSlopeTestData
#' @author Jack Leary
#' @description Creates a data.frame of \code{marge} model breakpoints, \emph{p}-values, and other info.
#' @param marge.model A \code{marge} model object, like those returned from \code{\link{marge2}}. Defaults to NULL.
#' @param pt A data.frame containing pseudotime or latent time values. Defaults to NULL.
#' @return A data.frame containing model data.
#' @seealso \code{\link{marge2}}
#' @examples
#' \dontrun{createSlopeTestData(marge_mod, pt_df)}

createSlopeTestData <- function(marge.model = NULL, pt = NULL) {
  # check inputs
  if (is.null(marge.model) | is.null(pt)) {
    stop("Input to createSlopeTestData() should be non-null.")
  }
  # run function
  # checks to see if marge model was set to "NA"try-error" in testDynamic()
  if (all(class(marge.model) == "try-error")) {
    rounded_brkpts <- NA_real_
    brkpts <- NA_real_
    brkpt_dirs <- NA_character_
    p_vals <- NA_real_
    mod_notes <- "MARGE model error"
  } else {
    if (length(coef(marge.model[[1]])) == 1) {
      rounded_brkpts <- NA_real_
      brkpts <- NA_real_
      brkpt_dirs <- NA_character_
      p_vals <- NA_real_
      mod_notes <- "No non-intercept coefficients"
    } else {
      # grab the k breakpoints from the MARGE model
      model_breakpoints_rounded <- extractBreakpoints(marge.model[[1]], directions = TRUE)
      rounded_brkpts <- model_breakpoints_rounded$Breakpoint
      brkpt_dirs <- model_breakpoints_rounded$Direction
      brkpts <- sapply(rounded_brkpts, function(x) pt$PT[which.min(abs(pt$PT - x))])
      # extract Wald test p-values from summary.glm()
      p_vals <- unname(summary(marge.model[[1]])$coefficients[, "Pr(>|z|)"][-1])
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
