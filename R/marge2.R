#' \code{marge2} models for negative binomial data.
#'
#' @name marge2
#' @description MARS fitting function for generalized linear models (GLMs).
#' @importFrom glm2 glm.fit2
#' @param X_pred A matrix of the predictor variables. Note that this matrix should include a column of 1's for the intercept term.
#' @param Y The response variable. A vector of length n by N.
#' @param pen A set penalty used for the GCV (note: MARGE doesn't actually use this). The default is 2.
#' @param tols_score The set tolerance for monitoring the convergence for the difference in score statistics between the parent and candidate model (this is the lack-of-fit criterion used for MARGE). The default is 0.00001
#' @param M A set threshold for the number of basis functions to be used. The default is 7.
#' @param minspan A set minimum span value. The default is \code{minspan = NULL}.
#' @param return.basis Whether the basis model matrix (denoted \code{B_final}) should be returned as part of the \code{marge} model object. Defaults to FALSE since it makes the mdoel object much larger than necessary.
#' @param return.wic Whether the WIC matrix and final WIC value should be returned as part of the \code{marge} model object. Defaults to FALSE because I don't ever use them.
#' @param return.GCV Whether the final GCV value should be returned as part of the \code{marge} model object. Defaults to FALSE because I never use it.
#' @details For further details please look at the \code{mars_ls} function - there are more details on the general MARS algorithm. MARGE will produce output for two penalties: 2 and log(N). A figure is automatically generated plotting WIC against the no. of parameters.
#' @return \code{marge} returns a list of calculated values consisting of:
#' @return \code{B_final}Tthe basis model matrix for the final model fit.
#' @return \code{wic_mat} A matrix of WIC values (with both penalties) for MARGE models given by the forward pass.
#' @return \code{min_wic} The WIC (with both penalties) for the final MARGE model.
#' @return \code{GCV} The GCV for the final selected model.
#' @return \code{y_pred} The fitted values from the final selected model (with both penalties).
#' @return \code{final_mod} The final selected (with both penalties) model matrix.
#' @author Jakub Stoklosa and David I. Warton.
#' @references Friedman, J. (1991). Multivariate adaptive regression splines. \emph{The Annals of Statistics}, \strong{19}, 1--67.
#' @references Stoklosa, J., Gibb, H. and Warton, D.I. (2014). Fast forward selection for generalized estimating equations with a large number of predictor variables. \emph{Biometrics}, \strong{70}, 110--120.
#' @references Stoklosa, J. and Warton, D.I. (2018). A generalized estimating equation approach to multivariate adaptive regression splines. \emph{Journal of Computational and Graphical Statistics}, \strong{27}, 245--253.
#' @seealso \code{\link{mars_ls}} and \code{\link{backward_sel_WIC}}
#' @export
#' @importFrom gamlss gamlss
#' @importFrom mvabund manyglm
#' @importFrom MASS glm.nb
#' @importFrom stats binomial poisson
#' @examples
#' marge(X_pred = pseudotime_df, Y = expr_vec)

marge2 <- function(X_pred = NULL,
                   Y = NULL,
                   pen = 2,
                   tols_score = 0.00001,
                   M = 7,
                   minspan = NULL,
                   return.basis = FALSE,
                   return.wic = FALSE,
                   return.GCV = FALSE) {
  # check inputs
  if (any(unlist(lapply(list(X_pred, Y), is.null)))) stop("Some inputs to marge() are missing.")
  NN <- length(Y)  # Total sample size
  q <- ncol(X_pred)  # Number of predictor variables

  # Algorithm 2 (forward pass) as in Friedman (1991). Uses score statistics instead of RSS, etc.

  B <- as.matrix(rep(1, NN))  # Start with the intercept model.

  colnames(B) <- c("Intercept")
  var_name_vec <- c("Intercept")
  var_name_list <- list("Intercept")
  B_names_vec <- c("Intercept")
  min_knot_vec <- c("Intercept")
  pred.name_vec <- c("Intercept")
  cut_vec <- c("Intercept")
  trunc.type_vec <- c(1)
  is.int_vec <- c("FALSE")
  mod_struct <- c(1)  # Univariate (1) or interaction (2).
  score_term <- c(0)
  TSS <- sum((Y - mean(Y))^2)
  GCV.null <- TSS / (NN * (1 - (1 / NN))^2)

  # Null model setup.

  m <- 1
  k <- 1
  breakFlag <- FALSE
  ok <- TRUE
  int.count <- 0

  while(ok) {
    if (breakFlag) break

    var.mod_temp <- c()
    score_term_temp <- c()
    min_knot_vec_temp <- c()
    int.count1_temp <- c()
    is.int_temp <- c()
    trunc.type_temp <- c()

    B_new_list_temp <- list()
    var_name_list1_temp <- list()
    B_names_temp <- list()
    X_red_temp <- list()
    B_temp_list <- list()

    # Obtain/calculate the null stats here (speeds things up).
    B_null_stats <- stat_out_score_glm_null(Y = Y, B_null = B)
    VS.est_list <- B_null_stats$VS.est_list
    A_list <- B_null_stats$A_list
    B1_list <- B_null_stats$B1_list
    mu.est <- B_null_stats$mu.est
    V.est <- B_null_stats$V.est

    for (v in seq(q)) {
      var_name <- colnames(X_pred)[v]
      X <- round(X_pred[, v], 4)
      #X <- X_pred[, v]

      X_red1 <- min_span(X_red = c(round(X, 4)), q = q, minspan = minspan)  # Reduce the space between knots.
      X_red2 <- max_span(X_pred = c(round(X, 4)), q = q)  # Truncate the ends of data to avoid extreme values.

      X_red <- intersect(X_red1, X_red2)

      score_knot_both_int_mat <- c()
      score_knot_both_add_mat <- c()
      score_knot_one_int_mat <- c()
      score_knot_one_add_mat <- c()

      int.count1 <- 0

      in.set <- ifelse(ncol(B) > 1, sum(!var_name_vec %in% var_name), 0)

      for (t in seq(length(X_red))) {
        b1_new <- matrix(tp1(x = X,  t = X_red[t]), ncol = 1)  # Pairs of truncated functions.
        b2_new <- matrix(tp2(x = X, t = X_red[t]), ncol = 1)

        score_knot_both_int <- c()
        score_knot_both_add <- c()
        score_knot_one_int <- c()
        score_knot_one_add <- c()

        if (in.set == 0) {
          B_new_both_add <- cbind(B, b1_new, b2_new)   # Additive model with both truncated functions.
          B_new_one_add <- cbind(B, b1_new)         # Additive model with one truncated function (positive part).

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

          score_knot_both_add <- c(score_knot_both_add, meas_model_both_add$score)
          score_knot_one_add <- c(score_knot_one_add, meas_model_one_add$score)

          score_knot_both_int <- score_knot_one_int <- -100000   # Interaction set is impossible since there is nothing to interact with, so let the LOF measure be a huge negative number.
        } else {
          var_name_struct <- which(((var_name != var_name_vec) * mod_struct) == 1)
          colnames(B)[1] <- c("")
          B2 <- as.matrix(B[, var_name_struct])
          if (k != 1 & (any(!var_name_vec[-1] %in% var_name))) B2 <- as.matrix(B2[, -1])
          if (ncol(B2) == 0) B2 <- as.matrix(B[, 1])

          for (nn in seq(ncol(B2))) {
            B2a <- matrix(rep(B2[, nn], 2), ncol = 2)
            B2b <- matrix(B2[, nn], ncol = 1)
            B_new_both_int <- cbind(B, B2a * cbind(b1_new, b2_new))
            B_new_one_int <- cbind(B, B2b * b1_new)  # Interaction model with one truncated function (i.e., the positive part).

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

            score_knot_both_int <- c(score_knot_both_int, meas_model_both_int$score)
            score_knot_one_int <- c(score_knot_one_int, meas_model_one_int$score)
          }

          B_new_both_add <- cbind(B, b1_new, b2_new)
          B_new_one_add <- cbind(B, b1_new)

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
            temp <- utils::tail(which(utils::tail(max(round(score_knot, 6), na.rm = T), n = 1) == round(score_knot, 6), arr.ind = TRUE), n = 1)
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
          if (utils::tail(max(score_knot_one_int_mat, na.rm = T), n = 1) > utils::tail(max(score_knot_both_add_mat, na.rm = T), n = 1)) {
            trunc.type <- 1
            if (utils::tail(max(score_knot_one_int_mat, na.rm = T), n = 1) > utils::tail(max(score_knot_one_add_mat, na.rm = T), n = 1)) {
              int <- TRUE
              score_knot <- score_knot_one_int_mat
              temp <- utils::tail(which(utils::tail(max(round(score_knot, 6), na.rm = T), n = 1) == round(score_knot, 6), arr.ind = T), n = 1)
              min_knot1 <- temp[1]
              best.var <- temp[2]
            } else {
              int <- FALSE
              score_knot <- score_knot_one_add_mat
              temp <- utils::tail(which(utils::tail(max(round(score_knot, 6), na.rm = T), n = 1) == round(score_knot, 6), arr.ind = T), n = 1)
              min_knot1 <- temp[1]
              best.var <- temp[2]
            }
          } else {
            int <- FALSE
            if (utils::tail(max(score_knot_both_add_mat, na.rm = T), n = 1) <= utils::tail(max(score_knot_one_add_mat, na.rm = T), n = 1)) {
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
            temp <- utils::tail(which(utils::tail(max(round(score_knot, 6), na.rm = T), n = 1) == round(score_knot, 6), arr.ind = T), n = 1)
            min_knot1 <- temp[1]
            best.var <- temp[2]
          }
          if (utils::tail(max(score_knot_both_int_mat, na.rm = T), n = 1) > utils::tail(max(score_knot_one_add_mat, na.rm = TRUE), n = 1)) {
            int <- TRUE
            trunc.type <- 2
            score_knot <- score_knot_both_int_mat
            temp <- utils::tail(which(utils::tail(max(round(score_knot, 6), na.rm = TRUE), n = 1) == round(score_knot, 6), arr.ind = TRUE), n = 1)
            min_knot1 <- temp[1]
            best.var <- temp[2]
          } else {
            int <- F
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
        colnames(B)[1] <- c("")
        var_name1 <- which(var_name_vec != var_name)
        var_name2 <- ifelse(int.count == 0, var_name_vec[var_name1], var_name_vec[var_name_struct])
        var_name2 <- var_name2[-1]
        var_name_struct <- mod_struct1[mod_struct1 %in% var_name1]
        B2 <- as.matrix(B[, var_name_struct, drop = FALSE])
        #B2 <- as.matrix(B2[, -1, drop = FALSE])
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
        B_names <- ifelse(trunc.type == 2, c(B_names, paste(B3_names[best.var], B_name2, sep = "*")), B_names)

        var_name_list1 <- list()
        for (ll in seq(ncol(B_new))) {
          colnames(B_new)[ll] <- paste(var_name, colnames(B_new)[ll], sep = ":")
          var_name_list1 <- c(var_name_list1, list(colnames(B_new)[ll]))
          int.count1 <- int.count1 + 1
        }
      } else {
        var_name_list1 <- list()
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

      meas_model <- score_fun_glm(Y = Y,
                                  VS.est_list = VS.est_list,
                                  A_list = A_list,
                                  B1_list = B1_list,
                                  mu.est = mu.est,
                                  V.est = V.est,
                                  B1 = B_temp,
                                  XA = B_new)

      score2 <- meas_model$score

      meas_model0 <- stat_out(Y = Y,
                              B1 = B_temp,
                              TSS = TSS,
                              GCV.null = GCV.null,
                              pen = pen)
      GCVq2 <- meas_model0$GCVq1

      if (GCVq2 < (-10) || round(score2, 4) <= 0) {
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
      }

      if (GCVq2 >= (-10) || round(score2, 4) > 0) {
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

    }  # Terminate the for () loop to end v (variables) here.

    if (breakFlag) break

    best.mod <- which.max(score_term_temp)  # Finds the best model (i.e., the max LOF) from your candidate model/basis set. This becomes the new parent.

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
    is.int_vec <- c(is.int_vec, int)

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

    if (nrow(B) <= (ncol(B) + 2)) {  # To avoid the p > N issue!
      ok <- FALSE
    }

    if (m >= M) { # If model exceeds no. of set terms, terminate it.
      ok <- FALSE
    }
  }

  colnames(B) <- B_names_vec
  B2 <- B

  # Algorithm 3 (backward pass) as in Friedman (1991) but for GLM/GEE use WIC.

  WIC_vec_2 <- NA

  full.wic <- 0

  B_new <- B
  ncol_Bnew <- ncol(B_new)
  ncol_B <- ncol(B)
  cnames_2 <- list(colnames(B_new))

  wic_mat_2 <- matrix(NA, ncol = ncol_B , nrow = ncol_B )
  colnames(wic_mat_2) <- colnames(B)
  wic_mat_2 <- cbind(wic_mat_2, rep(NA, ncol_B ))
  colnames(wic_mat_2)[(ncol_B  + 1)] <- "Forward pass model"

  wic_mat_2[1, (ncol_B + 1)] <- full.wic

  wic1_2 <- backward_sel_WIC(Y = Y, B = B_new)

  wic_mat_2[2, 2:(length(wic1_2) + 1)] <- wic1_2

  WIC_2 <- sum(apply(wic_mat_2[1:2, ], 1, min, na.rm = TRUE)) + 2 * ncol_Bnew

  WIC_vec_2 <- c(WIC_vec_2, WIC_2)

  variable.lowest_2 <- as.numeric(which(wic1_2 == min(wic1_2, na.rm = TRUE))[1])
  var.low.vec_2 <- c(colnames(B_new)[variable.lowest_2 + 1])
  B_new_2 <- as.matrix(B_new[, -(variable.lowest_2 + 1)])

  cnames_2 <- c(cnames_2, list(colnames(B_new_2)))

  for (i in 2:(ncol_B - 1)) {
    if (i != (ncol_B - 1)) {
      wic1_2 <- backward_sel_WIC(Y = Y, B_new = B_new_2)

      wic_mat_2[(i + 1), colnames(B_new_2)[-1]] <- wic1_2

      WIC_2 <- sum(apply(wic_mat_2[1:(i + 1), ], 1, min, na.rm = TRUE)) + 2 * ncol(B_new_2)

      WIC_vec_2 <- c(WIC_vec_2, WIC_2)

      variable.lowest_2 <- as.numeric(which(wic1_2 == min(wic1_2, na.rm = TRUE))[1])
      var.low.vec_2 <- c(var.low.vec_2, colnames(B_new_2)[variable.lowest_2 + 1])

      B_new_2 <- as.matrix(B_new_2[, -(variable.lowest_2 + 1)])
    } else {
      full.fit_2 <- gamlss::gamlss(Y ~ B_new_2 - 1, family = "NBI", trace = FALSE)
      sink(tempfile())
      full.wald_2 <- ((as.matrix(summary(full.fit_2))[, 3])[-c(1, nrow(as.matrix(summary(full.fit_2))))])^2
      sink()

      wic1_2 <- full.wald_2
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
  final_mod <- MASS::glm.nb(c(t(Y)) ~ B_final - 1, method = "glm.fit2", init.theta = 1, y = FALSE, model = FALSE)
  p_2 <- ncol(B_final)
  df1a <- p_2 + pen * (p_2 - 1) / 2  # This matches the earth() package, SAS and Friedman (1991) penalty.

  RSS1 <- sum((Y - stats::fitted(final_mod))^2)
  GCV1 <- RSS1 / (NN * (1 - df1a / NN)^2)

  min_wic_own <- min(wic_mat_2, na.rm = TRUE)

  final_mod <- stripGLM(glm.obj = final_mod)

  z <- NULL
  if (return.basis) { z$bx <- B_final }
  if (return.wic) { z$wic_mat <- wic_mat_2; z$min_wic_own <- min_wic_own; }
  if (return.GCV) { z$GCV <- GCV1 }
  z$final_mod <- final_mod
  class(z) <- "marge"

  return(z)
}
