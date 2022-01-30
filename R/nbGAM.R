#' Fit a negative-binomial GAM.
#'
#' @name nbGAM
#' @author Jack Leary
#' @description Fits a negative-binomial family GAM (glorified wrapper function).
#' @importFrom mgcv gam nb negbin
#' @param expr A vector of integer counts. Defaults to NULL.
#' @param pt A dataframe or vector of pesudotime values. Defaults to NULL.
#' @param theta.init (Optional) If you have a numeric value for \eqn{\theta} you'd like to use provide it here, else it can be estimated internally. Defaults to NULL.
#' @return An object of class \code{gam}
#' @export
#' @examples
#' \dontrun{nbGAM(expr_vec, pt_df)}
#' \dontrun{nbGAM(expr_vec, pt_vec, theta.init = 2.5)}

nbGAM <- function(expr = NULL, pt = NULL, theta.init = NULL) {
  # check inputs
  if (is.null(expr) || is.null(pt)) stop("Inputs to nbGAM() are missing.")
  if (class(pt) == "data.frame") {
    pt_vec <- pt[, 1]
  } else if (class(pt) == "numeric") {
    pt_vec <- pt
  } else {
    stop("Unsupported class for argument pt to nbGAM().")
  }
  # fit model
  if (is.null(theta.init)) {
    nb_gam_mod <- mgcv::gam(expr ~ s(pt_vec), family = mgcv::nb(), )
  } else {
    nb_gam_mod <- mgcv::gam(expr ~ s(pt_vec), family = mgcv::negbin(theta = theta.init))
  }
  return(nb_gam_mod)
}
