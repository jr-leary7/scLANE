#' Fit a negative-binomial GAM.
#'
#' @name nbGAM
#' @author Jack Leary
#' @description Fits a negative-binomial family GAM using a B-spline of degree 3 on pseudotime. If data are multi-subject in nature, a random intercept is included for each subject.
#' @importFrom gamlss gamlss random gamlss.control
#' @importFrom splines bs
#' @param expr A vector of integer counts. Defaults to NULL.
#' @param pt A dataframe or vector of pesudotime values. Defaults to NULL.
#' @param id.vec (Optional) A vector of subject IDs to be used in creating random intercepts in the GAM. Useful for comparing GAMs to GLMMs. Defaults to NULL.
#' @return An object of class \code{gamlss}
#' @seealso \code{\link[gamlss]{gamlss}}
#' @seealso \code{\link[gamlss]{random}}
#' @seealso \code{\link[gamlss.dist]{NBI}}
#' @export
#' @examples
#' \dontrun{
#' nbGAM(expr_vec, pt = pt_df)
#' nbGAM(expr_vec,
#'       pt = pt_df,
#'       id.vec = subject_ids)
#' }

nbGAM <- function(expr = NULL,
                  pt = NULL,
                  id.vec = NULL) {
  # check inputs
  if (is.null(expr) || is.null(pt)) { stop("Required inputs to nbGAM() are missing.") }
  if (inherits(pt, "data.frame")) {
    pt_vec <- pt[, 1]
  } else if (inherits(pt, "numeric")) {
    pt_vec <- pt
  } else {
    stop("Unsupported class for argument pt to nbGAM().")
  }
  # fit model
  if (!is.null(id.vec)) {
    mod_df <- data.frame(COUNT = expr,
                         PT = pt_vec,
                         ID = as.factor(id.vec))
    nb_gam_mod <- gamlss::gamlss(COUNT ~ splines::bs(PT, degree = 3) + gamlss::random(ID, df = length(unique(ID))),
                                 data = mod_df,
                                 family = "NBI",
                                 control = gamlss::gamlss.control(trace = FALSE))
  } else {
    mod_df <- data.frame(COUNT = expr, PT = pt_vec)
    nb_gam_mod <- gamlss::gamlss(COUNT ~ splines::bs(PT, degree = 3),
                                 data = mod_df,
                                 family = "NBI",
                                 control = gamlss::gamlss.control(trace = FALSE))
  }
  return(nb_gam_mod)
}
