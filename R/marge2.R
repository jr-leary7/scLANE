#' Fit \code{MARGE} models of single cell counts.
#'
#' @name marge2
#' @author Jakub Stoklosa
#' @author David I. Warton.
#' @author Jack Leary
#' @description MARS fitting function for negative binomial generalized linear models (GLMs) & generalized estimating equations (GEEs).
#' @import glm2
#' @import magrittr
#' @importFrom dplyr mutate ntile group_by sample_frac
#' @importFrom geeM geem
#' @importFrom MASS glm.nb negative.binomial theta.mm
#' @importFrom utils tail
#' @importFrom purrr map_chr
#' @importFrom stats fitted coef offset as.formula
#' @param X_pred A matrix of the predictor variables. Defaults to NULL.
#' @param Y The response variable. Defaults to NULL.
#' @param Y.offset (Optional) An vector of per-cell size factors to be included in the final model fit as an offset. Defaults to NULL.
#' @param M A set threshold for the number of basis functions to be used. Defaults to 5.
#' @param is.gee Should the \code{geeM} package be used to fit a negative binomial GEE? Defaults to FALSE.
#' @param id.vec If \code{is.gee = TRUE}, must be a vector of ID values for the observations. Data must be sorted such that the subjects are in order! Defaults to NULL.
#' @param cor.structure If \code{is.gee = TRUE}, must be a string specifying the desired correlation structure for the NB GEE. Defaults to "ar1".
#' @param approx.knot (Optional) Should the set of candidate knots be reduce in order to speed up computation? This has little effect on the final fit, but can improve computation time significantly. Defaults to TRUE.
#' @param n.knot.max (Optional) The maximum number of candidate knots to consider. Uses uniform sampling to set this number of unique values from the reduced set of all candidate knots. Defaults to 20.
#' @param tols_score (Optional) The set tolerance for monitoring the convergence for the difference in score statistics between the parent and candidate model (this is the lack-of-fit criterion used for MARGE). Defaults to 0.00001.
#' @param minspan (Optional) A set minimum span value. Defaults to NULL.
#' @param return.basis (Optional) Whether the basis model matrix should be returned as part of the \code{marge} model object. Defaults to FALSE.
#' @param return.WIC (Optional) Whether the WIC matrix should be returned as part of the \code{marge} model object. Defaults to FALSE.
#' @param return.GCV (Optional) Whether the final GCV value should be returned as part of the \code{marge} model object. Defaults to FALSE.
#' @details
#' \itemize{
#' \item If models are being fit using an offset (as is recommended), it is assumed that the offset represents a library size factor (or similar quantity) generated using e.g., \code{\link{createCellOffset}} or \code{\link[scuttle]{computeLibraryFactors}}. Since this quantity represents a scaling factor divide by sequencing depth, the offset is formulated as \code{offset(log(1 / cell_offset))}. The inversion is necessary because the rate term, i.e. the sequencing depth, is the denominator of the estimated size factors.
#' }
#' @return An object of class \code{marge} containing the fitted model & other optional quantities of interest (basis function matrix, GCV, etc.).
#' @references Friedman, J. (1991). Multivariate adaptive regression splines. \emph{The Annals of Statistics}, \strong{19}, 1--67.
#' @references Stoklosa, J., Gibb, H. and Warton, D.I. (2014). Fast forward selection for generalized estimating equations with a large number of predictor variables. \emph{Biometrics}, \strong{70}, 110--120.
#' @references Stoklosa, J. and Warton, D.I. (2018). A generalized estimating equation approach to multivariate adaptive regression splines. \emph{Journal of Computational and Graphical Statistics}, \strong{27}, 245--253.
#' @seealso \code{\link{backward_sel_WIC}}
#' @seealso \code{\link{testDynamic}}
#' @seealso \code{\link{createCellOffset}}
#' @seealso \code{\link[MASS]{glm.nb}}
#' @seealso \code{\link[geeM]{geem}}
#' @export
#' @examples
#' \dontrun{
#'   marge2(pseudotime_df,
#'          Y = expr_vec,
#'          M = 3)
#'   marge2(pseudotime_df,
#'          Y = expr_vec,
#'          Y.offset = size_factor_vec,
#'          is.gee = TRUE,
#'          id.vec = subject_vec,
#'          cor.structure = "exchangeable")
#'   marge2(pseudotime_df,
#'          Y = expr_vec,
#'          is.gee = TRUE,
#'          id.vec = subject_vec,
#'          cor.structure = "ar1",
#'          n.knot.max = 10,
#'          return.basis = TRUE)
#' }

marge2 <- function(X_pred = NULL,
                   Y = NULL,
                   Y.offset = NULL,
                   M = 5,
                   is.gee = FALSE,
                   id.vec = NULL,
                   cor.structure = "ar1",
                   approx.knot = TRUE,
                   n.knot.max = 20,
                   tols_score = 1e-5,
                   minspan = NULL,
                   return.basis = FALSE,
                   return.WIC = FALSE,
                   return.GCV = FALSE) {
  # check inputs
  if (is.null(X_pred) || is.null(Y)) { stop("Some required inputs to marge2() are missing.") }
  if (is.gee & is.null(id.vec)) { stop("id.vec in marge2() must be non-null if is.gee = TRUE.") }
  if (is.gee & (!cor.structure %in% c("independence", "exchangeable", "ar1", "unstructured"))) { stop("cor.structure in marge2() must be a known type if is.gee = TRUE.") }
  if (is.gee & is.unsorted(id.vec)) { stop("Your data must be ordered by subject, please do so before running marge2() with is.gee = TRUE.") }

  # Algorithm 2 (forward pass) as in Friedman (1991). Uses score statistics instead of RSS, etc.
  NN <- length(Y)  # Total sample size
  if (is.gee) {
    N <- length(unique(id.vec))
    n_vec <- as.numeric(table(id.vec))
  }
  q <- ncol(X_pred)  # Number of predictor variables
  B <- as.matrix(rep(1, NN))  # Start with the intercept model.
  # estimate "known" theta for GEE models
  if (is.gee) {
    theta_hat <- MASS::theta.mm(y = Y,
                                mu = mean(Y),
                                dfr = NN - 1)
  }

  pen <- 2  # penalty for GCV criterion
  colnames(B) <- "Intercept"
  var_name_vec <- "Intercept"
  var_name_list <- list("Intercept")
  B_names_vec <- "Intercept"
  min_knot_vec <- "Intercept"
  pred.name_vec <- "Intercept"
  cut_vec <- "Intercept"
  trunc.type_vec <- 1
  mod_struct <- 1  # Univariate (1) or interaction (2).
  score_term <- 0
  TSS <- sum((Y - mean(Y))^2)
  GCV.null <- TSS / (NN * (1 - (1 / NN))^2)

  # Null model setup.
  m <- k <- 1
  breakFlag <- FALSE
  ok <- TRUE
  int.count <- 0

  while (ok) {  # this is such egregiously bad code lol
    if (breakFlag) {
      break
    }

    var.mod_temp <- NULL
    score_term_temp <- NULL
    min_knot_vec_temp <- NULL
    int.count1_temp <- NULL
    is.int_temp <- NULL
    trunc.type_temp <- NULL
    B_new_list_temp <- vector("list")
    var_name_list1_temp <- vector("list")
    B_names_temp <- vector("list")
    X_red_temp <- vector("list")
    B_temp_list <- vector("list")

    # Obtain/calculate the null stats here (speeds things up).
    if (is.gee) {
      B_null_stats <- stat_out_score_gee_null(Y = Y,
                                              B_null = B,
                                              id.vec = id.vec,
                                              cor.structure = cor.structure,
                                              theta.hat = theta_hat)
      VS.est_list <- B_null_stats$VS.est_list
      AWA.est_list <- B_null_stats$AWA.est_list
      J2_list <- B_null_stats$J2_list
      J11.inv <- B_null_stats$J11.inv
      Sigma2_list <- B_null_stats$Sigma2_list
      JSigma11 <- B_null_stats$JSigma11
      mu.est <- B_null_stats$mu.est
      V.est <- B_null_stats$V.est
    } else {
      B_null_stats <- stat_out_score_glm_null(Y = Y, B_null = B)
      VS.est_list <- B_null_stats$VS.est_list
      A_list <- B_null_stats$A_list
      B1_list <- B_null_stats$B1_list
      mu.est <- B_null_stats$mu.est
      V.est <- B_null_stats$V.est
    }

    for (v in seq(q)) {
      var_name <- colnames(X_pred)[v]
      if (approx.knot) {
        # enhanced candidate knot selection
        X <- round(X_pred[, v], 4)
        X_red1 <- min_span(X_red = X, q = q, minspan = minspan)
        X_red2 <- max_span(X_red = X, q = q)
        X_red <- intersect(X_red1, X_red2)
        if (length(X_red) > n.knot.max) {
          X_red <- seq(min(X_red), max(X_red), length.out = n.knot.max)
        }
      } else {
        # original candidate knot selection from 2017 Stoklosa & Warton paper
        X <- round(X_pred[, v], 4)
        X_red1 <- min_span(X_red = X, q = q, minspan = minspan)
        X_red2 <- max_span(X_red = X, q = q)
        X_red <- intersect(X_red1, X_red2)
      }

      score_knot_both_int_mat <- NULL
      score_knot_both_add_mat <- NULL
      score_knot_one_int_mat <- NULL
      score_knot_one_add_mat <- NULL

      int.count1 <- 0

      in.set <- ifelse(ncol(B) > 1, sum(!var_name_vec %in% var_name), 0)

      for (t in seq(length(X_red))) {
        b1_new <- matrix(tp1(x = X,  t = X_red[t]), ncol = 1)  # Pairs of truncated functions.
        b2_new <- matrix(tp2(x = X, t = X_red[t]), ncol = 1)

        score_knot_both_int <- NULL
        score_knot_both_add <- NULL
        score_knot_one_int <- NULL
        score_knot_one_add <- NULL

        if (in.set == 0) {
          B_new_both_add <- cbind(B, b1_new, b2_new)  # Additive model with both truncated functions.
          B_new_one_add <- cbind(B, b1_new)  # Additive model with one truncated function (positive part).

          if (is.gee) {
            meas_model_both_add <- score_fun_gee(Y = Y,
                                                 N = N,
                                                 n_vec = n_vec,
                                                 VS.est_list = VS.est_list,
                                                 AWA.est_list = AWA.est_list,
                                                 J2_list = J2_list,
                                                 Sigma2_list = Sigma2_list,
                                                 J11.inv = J11.inv,
                                                 JSigma11 = JSigma11,
                                                 mu.est = mu.est,
                                                 V.est = V.est,
                                                 B1 = B_new_both_add,
                                                 XA = cbind(b1_new, b2_new))
            meas_model_one_add <- score_fun_gee(Y = Y,
                                                N = N,
                                                n_vec = n_vec,
                                                VS.est_list = VS.est_list,
                                                AWA.est_list = AWA.est_list,
                                                J2_list = J2_list,
                                                Sigma2_list = Sigma2_list,
                                                J11.inv = J11.inv,
                                                JSigma11 = JSigma11,
                                                mu.est = mu.est,
                                                V.est = V.est,
                                                B1 = B_new_one_add,
                                                XA = b1_new)
          } else {
            meas_model_both_add <- score_fun_glm(Y = Y,
                                                 VS.est_list = VS.est_list,
                                                 A_list = A_list,
                                                 B1_list = B1_list,
                                                 mu.est = mu.est,
                                                 V.est = V.est,
                                                 B1 = B_new_both_add,
                                                 XA = cbind(b1_new, b2_new))
            meas_model_one_add <- score_fun_glm(Y = Y,
                                                VS.est_list = VS.est_list,
                                                A_list = A_list,
                                                B1_list = B1_list,
                                                mu.est = mu.est,
                                                V.est = V.est,
                                                B1 = B_new_one_add,
                                                XA = b1_new)
          }

          score_knot_both_add <- c(score_knot_both_add, meas_model_both_add$score)
          score_knot_one_add <- c(score_knot_one_add, meas_model_one_add$score)
          score_knot_both_int <- score_knot_one_int <- -1e5  # Interaction set is impossible since there is nothing to interact with, so let the LOF measure be a huge negative number.

          } else {
          var_name_struct <- which(((var_name != var_name_vec) * mod_struct) == 1)
          colnames(B)[1] <- ""
          B2 <- as.matrix(B[, var_name_struct])
          if (k != 1 & any(!var_name_vec[-1] %in% var_name)) {
            B2 <- as.matrix(B2[, -1, drop = FALSE])
          }

          for (nn in seq(ncol(B2))) {
            B2a <- matrix(rep(B2[, nn], 2), ncol = 2)
            B2b <- matrix(B2[, nn], ncol = 1)
            B_new_both_int <- cbind(B, B2a * cbind(b1_new, b2_new))
            B_new_one_int <- cbind(B, B2b * b1_new)  # Interaction model with one truncated function (i.e., the positive part).

            if (is.gee) {
              meas_model_both_int <- score_fun_gee(Y = Y,
                                                   N = N,
                                                   n_vec = n_vec,
                                                   VS.est_list = VS.est_list,
                                                   AWA.est_list = AWA.est_list,
                                                   J2_list = J2_list,
                                                   Sigma2_list = Sigma2_list,
                                                   J11.inv = J11.inv,
                                                   JSigma11 = JSigma11,
                                                   mu.est = mu.est,
                                                   V.est = V.est,
                                                   B1 = B_new_both_int,
                                                   XA = B2a * cbind(b1_new, b2_new))
              meas_model_one_int <- score_fun_gee(Y = Y,
                                                  N = N,
                                                  n_vec = n_vec,
                                                  VS.est_list = VS.est_list,
                                                  AWA.est_list = AWA.est_list,
                                                  J2_list = J2_list,
                                                  Sigma2_list = Sigma2_list,
                                                  J11.inv = J11.inv,
                                                  JSigma11 = JSigma11,
                                                  mu.est = mu.est,
                                                  V.est = V.est,
                                                  B1 = B_new_one_int,
                                                  XA = B2b * b1_new)
            } else {
              meas_model_both_int <- score_fun_glm(Y = Y,
                                                   VS.est_list = VS.est_list,
                                                   A_list = A_list,
                                                   B1_list = B1_list,
                                                   mu.est = mu.est,
                                                   V.est = V.est,
                                                   B1 = B_new_both_int,
                                                   XA = B2a * cbind(b1_new, b2_new))
              meas_model_one_int <- score_fun_glm(Y = Y,
                                                  VS.est_list = VS.est_list,
                                                  A_list = A_list,
                                                  B1_list = B1_list,
                                                  mu.est = mu.est,
                                                  V.est = V.est,
                                                  B1 = B_new_one_int,
                                                  XA = B2b * b1_new)
            }
            score_knot_both_int <- c(score_knot_both_int, meas_model_both_int$score)
            score_knot_one_int <- c(score_knot_one_int, meas_model_one_int$score)
          }

          B_new_both_add <- cbind(B, b1_new, b2_new)
          B_new_one_add <- cbind(B, b1_new)

          if (is.gee) {
            meas_model_both_add <- score_fun_gee(Y = Y,
                                                 N = N,
                                                 n_vec = n_vec,
                                                 VS.est_list = VS.est_list,
                                                 AWA.est_list = AWA.est_list,
                                                 J2_list = J2_list,
                                                 Sigma2_list = Sigma2_list,
                                                 J11.inv = J11.inv,
                                                 JSigma11 = JSigma11,
                                                 mu.est = mu.est,
                                                 V.est = V.est,
                                                 B1 = B_new_both_add,
                                                 XA = cbind(b1_new, b2_new))
            meas_model_one_add <- score_fun_gee(Y = Y,
                                                N = N,
                                                n_vec = n_vec,
                                                VS.est_list = VS.est_list,
                                                AWA.est_list = AWA.est_list,
                                                J2_list = J2_list,
                                                Sigma2_list = Sigma2_list,
                                                J11.inv = J11.inv,
                                                JSigma11 = JSigma11,
                                                mu.est = mu.est,
                                                V.est = V.est,
                                                B1 = B_new_one_add,
                                                XA = b1_new)
          } else {
            meas_model_both_add <- score_fun_glm(Y = Y,
                                                 VS.est_list = VS.est_list,
                                                 A_list = A_list,
                                                 B1_list = B1_list,
                                                 mu.est = mu.est,
                                                 V.est = V.est,
                                                 B1 = B_new_both_add,
                                                 XA = cbind(b1_new, b2_new))
            meas_model_one_add <- score_fun_glm(Y = Y,
                                                VS.est_list = VS.est_list,
                                                A_list = A_list,
                                                B1_list = B1_list,
                                                mu.est = mu.est,
                                                V.est = V.est,
                                                B1 = B_new_one_add,
                                                XA = b1_new)
          }
          score_knot_both_add <- c(score_knot_both_add, meas_model_both_add$score)
          score_knot_one_add <- c(score_knot_one_add, meas_model_one_add$score)
        }
        score_knot_both_int_mat <- rbind(score_knot_both_int_mat, score_knot_both_int)
        score_knot_both_add_mat <- rbind(score_knot_both_add_mat, score_knot_both_add)
        score_knot_one_int_mat <- rbind(score_knot_one_int_mat, score_knot_one_int)
        score_knot_one_add_mat <- rbind(score_knot_one_add_mat, score_knot_one_add)
      }

      # See the LM code above in regards to what the conditions below actually do.

      if (all((apply(score_knot_both_int_mat, 1, is.na))) & all((apply(score_knot_one_int_mat, 1, is.na)))) {
        int <- FALSE
        if (any(!is.na(score_knot_both_add_mat)) & any(!is.na(score_knot_one_add_mat))) {
          if (utils::tail(max(score_knot_both_add_mat, na.rm = TRUE), n = 1) > utils::tail(max(score_knot_one_add_mat, na.rm = TRUE), n = 1)) {
            trunc.type <- 2
            score_knot <- score_knot_both_add_mat
            min_knot1 <- utils::tail(which(utils::tail(max(round(score_knot, 6), na.rm = TRUE), n = 1) == round(score_knot, 6), arr.ind = TRUE), n = 1)[1]
          } else {
            trunc.type <- 1
            score_knot <- score_knot_one_add_mat
            min_knot1 <- utils::tail(which(utils::tail(max(round(score_knot, 6), na.rm = TRUE), n = 1) == round(score_knot, 6), arr.ind = TRUE), n = 1)[1]
          }
        } else if (all(is.na(score_knot_both_add_mat)) & any(!is.na(score_knot_one_add_mat))) {
          trunc.type <- 1
          score_knot <- score_knot_one_add_mat
          min_knot1 <- utils::tail(which.max(round(score_knot, 6)), n = 1)
        } else if (any(!is.na(score_knot_both_add_mat)) & all(is.na(score_knot_one_add_mat))) {
          trunc.type <- 2
          score_knot <- score_knot_one_add_mat
          min_knot1 <- utils::tail(which.max(round(score_knot, 6)), n = 1)
        } else {
          breakFlag <- TRUE
          break
        }
      } else if (all((apply(score_knot_both_int_mat, 1, is.na))) & any(!(apply(score_knot_one_int_mat, 1, is.na)))) {
        if (all(is.na(score_knot_both_add_mat))) {
          trunc.type <- 1
          if (all(is.na(score_knot_one_add_mat))) {
            int <- TRUE
            score_knot <- score_knot_one_int_mat
            temp <- utils::tail(which(utils::tail(max(round(score_knot, 6), na.rm = TRUE), n = 1) == round(score_knot, 6), arr.ind = TRUE), n = 1)
            min_knot1 <- temp[1]
            best.var <- temp[2]
          } else {
            if (utils::tail(max(score_knot_one_int_mat, na.rm = TRUE), n = 1) > utils::tail(max(score_knot_one_add_mat, na.rm = TRUE), n = 1)) {
              int <- TRUE
              score_knot <- score_knot_one_int_mat
              temp <- utils::tail(which(utils::tail(max(round(score_knot, 6), na.rm = T), n = 1) == round(score_knot, 6), arr.ind = TRUE), n = 1)
              min_knot1 <- temp[1]
              best.var <- temp[2]
            } else {
              int <- FALSE
              score_knot <- score_knot_one_add_mat
              min_knot1 <- utils::tail(which.max(round(score_knot, 6)), n = 1)
            }
          }
        } else {
          if (utils::tail(max(score_knot_one_int_mat, na.rm = TRUE), n = 1) > utils::tail(max(score_knot_both_add_mat, na.rm = TRUE), n = 1)) {
            trunc.type <- 1
            if (utils::tail(max(score_knot_one_int_mat, na.rm = TRUE), n = 1) > utils::tail(max(score_knot_one_add_mat, na.rm = TRUE), n = 1)) {
              int <- TRUE
              score_knot <- score_knot_one_int_mat
              temp <- utils::tail(which(utils::tail(max(round(score_knot, 6), na.rm = TRUE), n = 1) == round(score_knot, 6), arr.ind = TRUE), n = 1)
              min_knot1 <- temp[1]
              best.var <- temp[2]
            } else {
              int <- FALSE
              score_knot <- score_knot_one_add_mat
              temp <- utils::tail(which(utils::tail(max(round(score_knot, 6), na.rm = TRUE), n = 1) == round(score_knot, 6), arr.ind = TRUE), n = 1)
              min_knot1 <- temp[1]
              best.var <- temp[2]
            }
          } else {
            int <- FALSE
            if (utils::tail(max(score_knot_both_add_mat, na.rm = TRUE), n = 1) <= utils::tail(max(score_knot_one_add_mat, na.rm = TRUE), n = 1)) {
              trunc.type <- 1
              score_knot <- score_knot_one_add_mat
              min_knot1 <- utils::tail(which.max(round(score_knot, 6)), n = 1)
            } else {
              trunc.type <- 2
              score_knot <- score_knot_both_add_mat
              min_knot1 <- utils::tail(which.max(round(score_knot, 6)), n = 1)
            }
          }
        }
      } else if (any(!(apply(score_knot_both_int_mat, 1, is.na))) & all((apply(score_knot_one_int_mat, 1, is.na)))) {
        if (all(is.na(score_knot_both_add_mat))) {
          if (all(is.na(score_knot_one_add_mat))) {
            int <- TRUE
            trunc.type <- 2
            score_knot <- score_knot_both_int_mat
            temp <- utils::tail(which(utils::tail(max(round(score_knot, 6), na.rm = TRUE), n = 1) == round(score_knot, 6), arr.ind = TRUE), n = 1)
            min_knot1 <- temp[1]
            best.var <- temp[2]
          }
          if (utils::tail(max(score_knot_both_int_mat, na.rm = TRUE), n = 1) > utils::tail(max(score_knot_one_add_mat, na.rm = TRUE), n = 1)) {
            int <- TRUE
            trunc.type <- 2
            score_knot <- score_knot_both_int_mat
            temp <- utils::tail(which(utils::tail(max(round(score_knot, 6), na.rm = TRUE), n = 1) == round(score_knot, 6), arr.ind = TRUE), n = 1)
            min_knot1 <- temp[1]
            best.var <- temp[2]
          } else {
            int <- FALSE
            trunc.type <- 1
            score_knot <- score_knot_one_add_mat
            min_knot1 <- utils::tail(which.max(round(score_knot, 6)), n = 1)
          }
        } else {
          if (utils::tail(max(score_knot_both_int_mat, na.rm = TRUE), n = 1) > utils::tail(max(score_knot_both_add_mat, na.rm = TRUE), n = 1)) {
            trunc.type <- 2
            if (utils::tail(max(score_knot_both_int_mat, na.rm = TRUE), n = 1) > utils::tail(max(score_knot_both_add_mat, na.rm = TRUE), n = 1)) {
              int <- TRUE
              score_knot <- score_knot_both_int_mat
              temp <- utils::tail(which(utils::tail(max(round(score_knot, 6), na.rm = TRUE), n = 1) == round(score_knot, 6), arr.ind = TRUE), n = 1)
              min_knot1 <- temp[1]
              best.var <- temp[2]
            } else {
              int <- FALSE
              score_knot <- score_knot_both_add_mat
              min_knot1 <- utils::tail(which(utils::tail(max(round(score_knot, 6), na.rm = TRUE), n = 1) == round(score_knot, 6), arr.ind = TRUE), n = 1)[1]
            }
          } else {
            int <- FALSE
            if (utils::tail(max(score_knot_both_add_mat, na.rm = TRUE), n = 1) <= utils::tail(max(score_knot_one_add_mat, na.rm = TRUE), n = 1)) {
              trunc.type <- 1
              score_knot <- score_knot_one_add_mat
              min_knot1 <- utils::tail(which.max(round(score_knot, 6)), n = 1)
            } else {
              trunc.type <- 2
              score_knot <- score_knot_both_add_mat
              min_knot1 <- utils::tail(which.max(round(score_knot, 6)), n = 1)
            }
          }
        }
      } else {
        if (utils::tail(max(score_knot_both_int_mat, na.rm = TRUE), n = 1) > utils::tail(max(score_knot_one_int_mat, na.rm = TRUE), n = 1)) {
          if (any(!is.na(score_knot_both_add_mat))) {
            if (utils::tail(max(score_knot_both_add_mat, na.rm = TRUE), n = 1) >= utils::tail(max(score_knot_both_int_mat, na.rm = TRUE), n = 1)) {
              int <- FALSE
              if (utils::tail(max(score_knot_both_add_mat, na.rm = TRUE), n = 1) > utils::tail(max(score_knot_one_add_mat, na.rm = TRUE), n = 1)) {
                trunc.type <- 2
                score_knot <- score_knot_both_add_mat
                min_knot1 <- utils::tail(which.max(round(score_knot, 6)), n = 1)
              } else {
                trunc.type <- 1
                score_knot <- score_knot_one_add_mat
                min_knot1 <- utils::tail(which.max(round(score_knot, 6)))
              }
            } else {
              if (utils::tail(max(score_knot_both_int_mat, na.rm = TRUE), n = 1) > utils::tail(max(score_knot_one_add_mat, na.rm = TRUE), n = 1)) {
                int <- TRUE
                trunc.type <- 2
                score_knot <- score_knot_both_int_mat
                temp <- utils::tail(which(utils::tail(max(round(score_knot, 6), na.rm = TRUE), n = 1) == round(score_knot, 6), arr.ind = TRUE), n = 1)
                min_knot1 <- temp[1]
                best.var <- temp[2]
              } else {
                int <- FALSE
                trunc.type <- 1
                score_knot <- score_knot_one_add_mat
                min_knot1 <- utils::tail(which.max(round(score_knot, 6)), n = 1)
              }
            }
          } else {
            if (utils::tail(max(score_knot_one_add_mat, na.rm = TRUE), n = 1) >= utils::tail(max(score_knot_both_int_mat, na.rm = TRUE), n = 1)) {
              int <- FALSE
              trunc.type <- 1
              score_knot <- score_knot_one_add_mat
              min_knot1 <- utils::tail(which.max(round(score_knot, 6)), n = 1)
            } else {
              int <- TRUE
              trunc.type <- 2
              score_knot <- score_knot_both_int_mat
              temp <- utils::tail(which(utils::tail(max(round(score_knot, 6), na.rm = TRUE), n = 1) == round(score_knot, 6), arr.ind = TRUE), n = 1)
              min_knot1 <- temp[1]
              best.var <- temp[2]
            }
          }
        } else {
          if (any(!is.na(score_knot_both_add_mat))) {
            if (utils::tail(max(score_knot_both_add_mat, na.rm = TRUE), n = 1) >= utils::tail(max(score_knot_one_int_mat, na.rm = TRUE), n = 1)) {
              int <- FALSE
              if (utils::tail(max(score_knot_both_add_mat, na.rm = TRUE), n = 1) <= utils::tail(max(score_knot_one_add_mat, na.rm = TRUE), n = 1)) {
                trunc.type <- 1
                score_knot <- score_knot_one_add_mat
                min_knot1 <- utils::tail(which.max(round(score_knot, 6)), n = 1)
              } else {
                trunc.type <- 2
                score_knot <- score_knot_both_add_mat
                min_knot1 <- utils::tail(which.max(round(score_knot, 6)), n = 1)
              }
            } else {
              trunc.type <- 1
              if (utils::tail(max(score_knot_one_int_mat, na.rm = TRUE), n = 1) > utils::tail(max(score_knot_one_add_mat, na.rm = TRUE), n = 1)) {
                int <- TRUE
                score_knot <- score_knot_one_int_mat
                temp <- utils::tail(which(utils::tail(max(round(score_knot, 6), na.rm = TRUE), n = 1) == round(score_knot, 6), arr.ind = TRUE), n = 1)
                min_knot1 <- temp[1]
                best.var <- temp[2]
              } else {
                int <- FALSE
                score_knot <- score_knot_one_add_mat
                min_knot1 <- utils::tail(which.max(round(score_knot, 6)), n = 1)
              }
            }
          } else {
            trunc.type <- 1
            if (all(is.na(score_knot_one_add_mat))) {
              int <- TRUE
              score_knot <- score_knot_one_int_mat
              temp <- utils::tail(which(utils::tail(max(round(score_knot, 6), na.rm = TRUE), n = 1) == round(score_knot, 6), arr.ind = TRUE), n = 1)
              min_knot1 <- temp[1]
              best.var <- temp[2]
            } else {
              if (utils::tail(max(score_knot_one_add_mat, na.rm = TRUE), n = 1) >= utils::tail(max(score_knot_one_int_mat, na.rm = TRUE), n = 1)) {
                int <- FALSE
                score_knot <- score_knot_one_add_mat
                min_knot1 <- utils::tail(which.max(round(score_knot, 6)))
              } else {
                int <- TRUE
                score_knot <- score_knot_one_int_mat
                temp <- utils::tail(which(utils::tail(max(round(score_knot, 6), na.rm = TRUE), n = 1) == round(score_knot, 6), arr.ind = TRUE), n = 1)
                min_knot1 <- temp[1]
                best.var <- temp[2]
              }
            }
          }
        }
      }

      b1_new <- matrix(tp1(x = X, t = X_red[min_knot1]), ncol = 1)
      b2_new <- matrix(tp2(x = X, t = X_red[min_knot1]), ncol = 1)
      colnames(b1_new) <- colnames(b2_new) <- var_name

      B_name1 <- paste("(", var_name, "-", signif(X_red[min_knot1], 4), ")", sep = "")
      B_name2 <- paste("(", signif(X_red[min_knot1], 4), "-", var_name, ")", sep = "")

      if (int) {
        mod_struct1 <- which(mod_struct == 1)
        colnames(B)[1] <- ""
        var_name1 <- which(var_name_vec != var_name)
        var_name2 <- ifelse(int.count == 0,
                            var_name_vec[var_name1],
                            var_name_vec[var_name_struct])
        var_name2 <- var_name2[-1]
        var_name_struct <- mod_struct1[mod_struct1 %in% var_name1]
        B2 <- as.matrix(B[, var_name_struct, drop = FALSE])
        B3_names <- B_names_vec[var_name_struct][-1]

        if (trunc.type == 2) {
          B2a <- matrix(rep(B2[, best.var], 2), ncol = 2)
          B_temp <- cbind(B, B2a*cbind(b1_new, b2_new))
          B_new <- B2a * cbind(b1_new, b2_new)
          var_name3 <- var_name2[best.var]
          colnames(B_new) <- rep(var_name3, 2)
        } else {
          B2b <- matrix(B2[, best.var], ncol = 1)
          B_temp <- cbind(B, B2b * b1_new)  # Interaction model with one truncated basis function (i.e., the positive part).
          B_new <- B2b * b1_new
          colnames(B_new) <- rep(var_name2[1], 1)
          var_name3 <- var_name2[best.var]
          colnames(B_new) <- rep(var_name3, 1)
        }

        B_names <- paste(B3_names[best.var], B_name1, sep = "*")
        B_names <- ifelse(trunc.type == 2,
                          c(B_names, paste(B3_names[best.var], B_name2, sep = "*")),
                          B_names)

        var_name_list1 <- vector("list")
        for (ll in seq(ncol(B_new))) {
          colnames(B_new)[ll] <- paste(var_name, colnames(B_new)[ll], sep = ":")
          var_name_list1 <- c(var_name_list1, list(colnames(B_new)[ll]))
          int.count1 <- int.count1 + 1
        }
      } else {
        var_name_list1 <- vector("list")
        if (trunc.type == 2) {
          B_temp <- cbind(B, b1_new, b2_new) # Additive model with both truncated basis functions.
          B_new <- cbind(b1_new, b2_new)
          B_names <- c(B_name1, B_name2)
          var_name_list1 <- c(var_name_list1, list(var_name))
          var_name_list1 <- c(var_name_list1, list(var_name))  # Repeat it because there are two new truncated basis function in the set.
        } else {
          B_temp <- cbind(B, b1_new) # Additive model with one truncated basis function (i.e., the positive part).
          B_new <- b1_new
          B_names <- B_name1
          var_name_list1 <- c(var_name_list1, list(var_name))
        }
      }

      if (is.gee) {
        meas_model <- score_fun_gee(Y = Y,
                                    N = N,
                                    n_vec = n_vec,
                                    VS.est_list = VS.est_list,
                                    AWA.est_list = AWA.est_list,
                                    J2_list = J2_list,
                                    Sigma2_list = Sigma2_list,
                                    J11.inv = J11.inv,
                                    JSigma11 = JSigma11,
                                    mu.est = mu.est,
                                    V.est = V.est,
                                    B1 = B_temp,
                                    XA = B_new)
      } else {
        meas_model <- score_fun_glm(Y = Y,
                                    VS.est_list = VS.est_list,
                                    A_list = A_list,
                                    B1_list = B1_list,
                                    mu.est = mu.est,
                                    V.est = V.est,
                                    B1 = B_temp,
                                    XA = B_new)
      }

      score2 <- meas_model$score
      meas_model0 <- stat_out(Y = Y,
                              B1 = B_temp,
                              TSS = TSS,
                              GCV.null = GCV.null,
                              pen = pen)
      if (meas_model0$GCVq1 < -10 || round(score2, 4) <= 0) {
        var.mod_temp <- c(var.mod_temp, NA)
        min_knot_vec_temp <- c(min_knot_vec_temp, NA)
        int.count1_temp <- c(int.count1_temp, NA)
        is.int_temp <- c(is.int_temp, int)

        trunc.type_temp <- c(trunc.type_temp, NA)
        X_red_temp <- c(X_red_temp, list(NA))
        B_new_list_temp <- c(B_new_list_temp, list(NA))
        var_name_list1_temp <- c(var_name_list1_temp, list(NA))
        B_names_temp <- c(B_names_temp, list(NA))
        B_temp_list <- c(B_temp_list, list(NA))
        score_term_temp <- c(score_term_temp, NA)
        if (length(var.mod_temp) == q) {
          breakFlag <- TRUE
          break
        } else {
          next
        }
      } else if (meas_model0$GCVq1 >= -10 || round(score2, 4) > 0) {
        score_term_temp <- c(score_term_temp, score2)
      }

      var.mod_temp <- c(var.mod_temp, var_name)
      min_knot_vec_temp <- c(min_knot_vec_temp, min_knot1)
      int.count1_temp <- c(int.count1_temp, int.count1)
      is.int_temp <- c(is.int_temp, int)
      trunc.type_temp <- c(trunc.type_temp, trunc.type)
      B_new_list_temp <- c(B_new_list_temp, list(B_new))
      var_name_list1_temp <- c(var_name_list1_temp, list(var_name_list1))
      B_names_temp <- c(B_names_temp, list(B_names))
      X_red_temp <- c(X_red_temp, list(X_red))
      B_temp_list <- c(B_temp_list, list(B_temp))
    }

    if (breakFlag) {
      break
    }

    best.mod <- which.max(score_term_temp)  # Finds the best model (i.e., the max LOF) from your candidate model/basis set. This becomes the new parent
    score2 <- score_term_temp[best.mod]
    min_knot_vec1 <- min_knot_vec_temp[best.mod]
    int.count1 <- int.count1_temp[best.mod]
    int <- is.int_temp[best.mod]
    trunc.type <- trunc.type_temp[best.mod]
    B_new <- B_new_list_temp[[best.mod]]
    var_name_list1 <- var_name_list1_temp[[best.mod]]
    B_names <- B_names_temp[[best.mod]]
    X_red <- X_red_temp[[best.mod]]
    B_temp <- B_temp_list[[best.mod]]
    score_term <- c(score_term, score2)
    min_knot_vec <- c(min_knot_vec, min_knot_vec1)
    pred.name_vec <- c(pred.name_vec, colnames(B_new)[1])
    cut_vec <- c(cut_vec, X_red[min_knot_vec1])
    trunc.type_vec <- c(trunc.type_vec, trunc.type)

    if (score_term[k + 1] < tols_score) {
      breakFlag <- TRUE
      break
    } else {
      if (int) {
        mod_struct <- c(mod_struct, rep(c(rep(2, int.count1 / trunc.type)), trunc.type))
        int.count <- int.count + 1
      } else {
        mod_struct <- c(mod_struct, rep(1, trunc.type))
      }
      B <- B_temp
      var_name_vec <- c(var_name_vec, colnames(B_new))
      var_name_list <- c(var_name_list, var_name_list1)
      B_names_vec <- c(B_names_vec, B_names)
      k <- k + 1
      m <- m + 2
    }

    # terminate if high dimensionality occurs or maximum # of hinge functions is reached
    if (nrow(B) <= (ncol(B) + 2) || m >= M) {
      ok <- FALSE
    }
  }

  # Algorithm 3 (backward pass) as in Friedman (1991) but for GLM/GEE use WIC.

  colnames(B) <- B_names_vec
  WIC_vec_2 <- NA
  full.wic <- 0
  B_new <- B
  ncol_B <- ncol(B)
  cnames_2 <- list(colnames(B_new))

  wic_mat_2 <- matrix(NA_real_, ncol = ncol_B , nrow = ncol_B)
  colnames(wic_mat_2) <- B_names_vec
  wic_mat_2 <- cbind(wic_mat_2, rep(NA_real_, ncol_B))
  colnames(wic_mat_2)[(ncol_B  + 1)] <- "Forward pass model"

  wic_mat_2[1, (ncol_B + 1)] <- full.wic
  wic1_2 <- backward_sel_WIC(Y = Y, B_new = B_new)
  wic_mat_2[2, 2:(length(wic1_2) + 1)] <- wic1_2
  WIC_2 <- sum(apply(wic_mat_2[1:2, ], 1, min, na.rm = TRUE)) + 2 * ncol(B_new)
  WIC_vec_2 <- c(WIC_vec_2, WIC_2)

  variable.lowest_2 <- as.numeric(which(wic1_2 == min(wic1_2, na.rm = TRUE))[1])
  var.low.vec_2 <- c(colnames(B_new)[variable.lowest_2 + 1])
  B_new_2 <- as.matrix(B_new[, -(variable.lowest_2 + 1)])
  cnames_2 <- c(cnames_2, list(colnames(B_new_2)))
  for (i in 2:(ncol_B - 1)) {
    wic1_2 <- backward_sel_WIC(Y = Y, B_new = B_new_2)
    if (i != (ncol_B - 1)) {
      wic_mat_2[(i + 1), colnames(B_new_2)[-1]] <- wic1_2
      WIC_2 <- sum(apply(wic_mat_2[1:(i + 1), ], 1, min, na.rm = TRUE)) + 2 * ncol(B_new_2)
      WIC_vec_2 <- c(WIC_vec_2, WIC_2)
      variable.lowest_2 <- as.numeric(which(wic1_2 == min(wic1_2, na.rm = TRUE))[1])
      var.low.vec_2 <- c(var.low.vec_2, colnames(B_new_2)[variable.lowest_2 + 1])
      B_new_2 <- as.matrix(B_new_2[, -(variable.lowest_2 + 1)])
    } else {
      wic_mat_2[(i + 1), colnames(B_new_2)[-1]] <- wic1_2
      WIC_2 <- sum(apply(wic_mat_2[1:(ncol_B), ], 1, min, na.rm = TRUE)) + 2 * ncol(B_new_2)
      WIC_vec_2 <- c(WIC_vec_2, WIC_2)
      B_new_2 <- as.matrix(B_new_2[, -(variable.lowest_2)])
      colnames(B_new_2) <- "Intercept"
    }
    cnames_2 <- c(cnames_2, list(colnames(B_new_2)))
  }

  # Some final model output, WIC, GCV etc.
  B_final <- as.matrix(B[, colnames(B) %in% cnames_2[[which.min(WIC_vec_2)]]])
  model_df <- as.data.frame(B_final)
  clean_varnames <- gsub("\\(", "", colnames(model_df)[-1])
  clean_varnames <- gsub("\\)", "", clean_varnames)
  clean_varnames <- gsub("-", "_", clean_varnames)
  clean_varnames <- paste0("h_", clean_varnames)
  colnames(model_df) <- c("Intercept", clean_varnames)
  model_formula <- paste(colnames(model_df), collapse = " + ")
  model_formula <- paste0("Y ~ -1 + ", model_formula)
  if (!is.null(Y.offset)) {
    model_df <- dplyr::mutate(model_df,
                              cell_offset = Y.offset)
    model_formula <- paste0(model_formula, " + ", "offset(log(1 / cell_offset))")  # this is correct i promise
  }
  model_formula <- stats::as.formula(model_formula)
  if (is.gee) {
    final_mod <- geeM::geem(model_formula,
                            data = model_df,
                            id = id.vec,
                            family = MASS::negative.binomial(theta_hat, link = log),
                            corstr = cor.structure,
                            scale.fix = FALSE,
                            sandwich = TRUE)
  } else {
    final_mod <- MASS::glm.nb(model_formula,
                              data = model_df,
                              method = "glm.fit2",
                              link = log,
                              init.theta = 1,
                              y = FALSE,
                              model = FALSE)
  }
  # format results
  if (!is.gee) {
    final_mod <- stripGLM(glm.obj = final_mod)
  }
  res <- list(final_mod = final_mod,
              basis_mtx = NULL,
              WIC_mtx = NULL,
              GCV = NULL,
              model_type = ifelse(is.gee, "GEE", "GLM"),
              coef_names = names(stats::coef(final_mod)),
              marge_coef_names = colnames(B_final))
  if (return.basis) {
    res$basis_mtx <- model_df
  }
  if (return.GCV) {
    df1a <- ncol(B_final) + pen * (ncol(B_final) - 1) / 2  # This matches the {earth} package, SAS and Friedman (1991) penalty
    res$GCV <- sum((Y - stats::fitted(final_mod))^2) / (NN * (1 - df1a / NN)^2)
  }
  if (return.WIC) {
    res$WIC_mtx <- wic_mat_2
  }
  class(res) <- "marge"
  return(res)
}
