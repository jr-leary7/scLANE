#' Fit a negative-binomial GAM.
#'
#' @name nbGAM
#' @author Jack Leary
#' @description Fits a negative-binomial family GAM using a cubic basis spline on pseudotime. If data are multi-subject in nature, a random intercept is included for each subject.
#' @importFrom dplyr mutate
#' @importFrom stats offset
#' @importFrom gamlss gamlss random gamlss.control pb
#' @importFrom splines bs
#' @param expr A vector of integer counts. Defaults to NULL.
#' @param pt A dataframe or vector of pseudotime values. Defaults to NULL.
#' @param Y.offset (Optional) An offset to be included in the final model fit. Defaults to NULL.
#' @param id.vec (Optional) A vector of subject IDs to be used in creating random intercepts in the GAM. Useful for comparing GAMs to GLMMs. Defaults to NULL.
#' @param penalize.spline (Optional) Should a P-spline be used to fit the GAM? Otherwise the default cubic basis spline is used instead. Defaults to FALSE.
#' @param spline.df (Optional) Degrees of freedom of the cubic basis spline. Unused if a P-spline is being fit, since it's estimated internally. Defaults to 5.
#' @return An object of class \code{gamlss}
#' @seealso \code{\link[gamlss]{gamlss}}
#' @seealso \code{\link[gamlss]{random}}
#' @seealso \code{\link[splines]{bs}}
#' @seealso \code{\link[gamlss]{pb}}
#' @seealso \code{\link[gamlss.dist]{NBI}}
#' @export
#' @examples
#' \dontrun{
#' nbGAM(expr_vec, pt_df)
#' nbGAM(expr_vec,
#'       pt = pt_df,
#'       id.vec = subject_ids,
#'       random.slopes = TRUE)
#' nbGAM(expr_vec,
#'       pt = pt_df,
#'       Y.offset = size_factor_vec,
#'       penalize.spline = TRUE,
#'       spline.df = 10)
#' }

nbGAM <- function(expr = NULL,
                  pt = NULL,
                  Y.offset = NULL,
                  id.vec = NULL,
                  penalize.spline = FALSE,
                  spline.df = 5) {
  # check inputs
  if (is.null(expr) || is.null(pt)) { stop("Required inputs to nbGAM() are missing.") }
  if (inherits(pt, "data.frame")) {
    pt_vec <- pt[, 1]
  } else if (inherits(pt, "numeric")) {
    pt_vec <- pt
  } else {
    stop("Unsupported class for argument pt to nbGAM().")
  }
  if (!is.null(id.vec)) {
    mod_df <- data.frame(COUNT = expr,
                         PT = pt_vec,
                         ID = as.factor(id.vec))
  } else {
    mod_df <- data.frame(COUNT = expr, PT = pt_vec)
  }
  if (!is.null(Y.offset)) {
    mod_df <- dplyr::mutate(mod_df,
                            cell_offset = Y.offset)
  }
  # set up model formula w/ desired parameters
  model_formula <- "COUNT ~ "
  if (penalize.spline) {
    model_formula <- paste0(model_formula, "pb(PT)")
  } else {
    model_formula <- paste0(model_formula, "bs(PT, df = ", spline.df, ")")
  }
  if (!is.null(id.vec)) {
    model_formula <- paste0(model_formula, " + random(ID, df = ", length(unique(id.vec)), ")")
  }
  if (!is.null(Y.offset)) {
    model_formula <- paste0(model_formula, " + offset(log(1 / cell_offset))")
  }
  model_formula <- stats::as.formula(model_formula)
  # fit additive model
  nb_gam_mod <- gamlss::gamlss(model_formula,
                               data = mod_df,
                               family = "NBI",
                               control = gamlss::gamlss.control(trace = FALSE))
  return(nb_gam_mod)
}
