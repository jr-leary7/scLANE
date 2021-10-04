#' MARS model for I.I.D. Gaussian data.
#'
#' @name mars_ls
#' @description MARS estimation function for independent Gaussian response data. If non-Gaussian response data is used and a family function is specified (i.e., a GLM is used) then \code{mars_ls} applies a \code{glm} on the final selected basis function.
#' @param X_pred : a matrix of the predictor variables (note that this matrix excludes the intercept term).
#' @param Y : the response variable.
#' @param pen : a set penalty that is used for the GCV. The default is 2.
#' @param tols : a set tolerance for monitoring the convergence for the RSS between the parent and candidate model (this is the lack-of-fit criterion used for MARS). The default is 0.00001.
#' @param M : a set threshold for the number of basis functions to be used. Defaults to 7.
#' @param minspan : a set minimum span value. The default is \code{minspan = NULL}.
#' @details Note that no argument is provided for an "additive structure only" model but one could just use \code{mgcv::gam()} instead.
#' @return \code{mars_ls} returns a list of calculated values consisting of:
#' @return \code{B_final} : the final selected basis model matrix.
#' @return \code{GCV_mat} : a matrix of GCV values for each fitted model.
#' @return \code{min_GCV_own} : the GCV for the final selected model.
#' @return \code{y_pred} : the fitted values from the final selected model.
#' @return \code{final_mod} : the final selected model matrix.
#' @author Jakub Stoklosa and David I. Warton.
#' @references Friedman, J. (1991). Multivariate adaptive regression splines. \emph{The Annals of Statistics}, \strong{19}, 1--67.
#' @references Milborrow, S. (2017a). Notes on the \code{earth} package. Package vignette. Available at: \url{http://127.0.0.1:31355/library/earth/doc/earth-notes.pdf}.
#' @references Milborrow, S. (2017b). \code{earth}: Multivariate Adaptive Regression Splines. R package version 4.4.7. Available at \url{http://CRAN.R-project.org/package = earth.}
#' @references Stoklosa, J. and Warton, D.I. (2018). A generalized estimating equation approach to multivariate adaptive regression splines. \emph{Journal of Computational and Graphical Statistics}, \strong{27}, 245--253.
#' @importFrom stats binomial poisson
#' @seealso \code{\link{marge}}
#' @examples
#' mars_ls(X_pred = pseudotime_df, Y = expr_vec)

mars_ls <- function(X_pred = NULL,
                    Y = NULL,
                    pen = 2,
                    tols = 0.00001,
                    M = 7,
                    minspan = NULL) {
  # check inputs
  if (any(sapply(c(X_pred, Y), is.null))) stop("Some inputs to mars_ls() are missing.")
  N <- length(Y)    # Sample size.
  q <- ncol(X_pred)  # No. of predictor variables.

  # Algorithm 2 (forward pass) as in Friedman (1991). Build an overfitted model using Rsq as the lack-of-fit (LOF).

  B <- as.matrix(rep(1, N))  # Start with the intercept model.

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

  # Null model setup - find the RSS, GCV and df. So initialize everything!

  B_null <- as.matrix(rep(1, N))
  RSS_term <- c(sum((Y - stats::fitted(stats::lm.fit(B_null, Y)))^2))
  TSS <- sum((Y - mean(Y))^2)
  RSSq_term <- c(1 - RSS_term[1] / TSS)
  df1 <- 1  # Set degrees freedom to 1.
  GCV_term <- c(RSS_term / (N * (1 - df1 / N)^2))
  GCVq_term <- c(0)
  GCV.null <- GCV_term[1]
  m <- 1
  k <- 1
  breakFlag <- FALSE
  ok <- TRUE
  int.count <- 0

  # Starting with loop 1 of algorithm considers terms (M). If model exceeds no. of set M we then terminate the loop.

  while(ok) {
    if (breakFlag) break

    var.mod_temp <- c()
    RSSq_term_temp <- c()
    GCVq_term_temp <- c()
    min_knot_vec_temp <- c()
    int.count1_temp <- c()
    is.int_temp <- c()
    trunc.type_temp <- c()

    B_new_list_temp <- list()
    var_name_list1_temp <- list()
    B_names_temp <- list()
    X_red_temp <- list()
    B_temp_list <- list()

    # Loop 2 of algorithm considers all predictor variables in the data (q). Start with first variable.
    # Note that variables can be repeated in the model, but cannot interact with itself (it doesn't
    # really makes sense for it too anyway)!

    for (v in seq(q)) {
      var_name <- colnames(X_pred)[v]
      X <- round(X_pred[, v], 4)  # Round to 4 digits (to match earth()). X starts the new basis.

      # These next few bits truncate the knot selection space: (1) reduce the space between knots and (2) chop off ends bits.

      X_red1 <- min_span(X_red = c(round(X, 4)), q = q, minspan = minspan)  # Reduce the space between knots.
      X_red2 <- max_span(X_pred = c(round(X, 4)), q = q)  # Truncate the ends of data to avoid extreme values.
      X_red <- intersect(X_red1, X_red2)

      RSSq_knot_both_int_mat <- c()
      RSSq_knot_both_add_mat <- c()
      RSSq_knot_one_int_mat <- c()
      RSSq_knot_one_add_mat <- c()
      GCVq_knot_both_int_mat <- c()
      GCVq_knot_both_add_mat <- c()
      GCVq_knot_one_int_mat <- c()
      GCVq_knot_one_add_mat <- c()

      int.count1 <- 0

      # This next bit is an indicator to check if the new variable can interact with any other variable already in the set. 0 = no and >0, possible.

      if (ncol(B) > 1) {
        in.set <- sum(!var_name_vec%in%var_name)
      } else if (ncol(B) == 1) {
        in.set <- 0
      }

      # Loop 3 of algorithm considers ALL the knot locations (after trim) of the chosen variable in the current loop. The strategy is simple - for
      # each knot find the lack-of-fit measures, and compare these with seperate fitted additive and interactions models.

      for (t in seq(length(X_red))) {
        b1_new <- matrix(tp1(x = X, t = X_red[t]), ncol = 1)  # Pairs of truncated power basis functions, positive and negative functions (first or both must be considered).
        b2_new <- matrix(tp2(x = X, t = X_red[t]), ncol = 1)

        RSSq_knot_both_int <- c()
        RSSq_knot_both_add <- c()
        RSSq_knot_one_int <- c()
        RSSq_knot_one_add <- c()
        GCVq_knot_both_int <- c()
        GCVq_knot_both_add <- c()
        GCVq_knot_one_int <- c()
        GCVq_knot_one_add <- c()

        # Find statistics for additive models only here, in fact this is only used once after initial entry of first variable, and used again if first variable was chosen.

        if (in.set == 0) {
          B_new_both_add <- cbind(B, b1_new, b2_new)  # Additive model with both truncated functions.
          B_new_one_add <- cbind(B, b1_new)  # Additive model with one truncated function (positive part).
          # calculate stats
          meas_model_both_add <- stat_out(Y = Y, B1 = B_new_both_add, TSS = TSS, GCV.null = GCV.null, pen = pen)
          meas_model_one_add <- stat_out(Y = Y, B1 = B_new_one_add, TSS = TSS, GCV.null = GCV.null, pen = pen)
          # save stats
          RSSq_knot_both_add <- c(RSSq_knot_both_add, meas_model_both_add$RSSq1)
          RSSq_knot_one_add <- c(RSSq_knot_one_add, meas_model_one_add$RSSq1)
          GCVq_knot_both_add <- c(GCVq_knot_both_add, meas_model_both_add$GCVq1)
          GCVq_knot_one_add <- c(GCVq_knot_one_add, meas_model_one_add$GCVq1)
          # Since no other variables in the set then no interaction possible ---> let the measures be some huge negative number
          GCVq_knot_both_int <- GCVq_knot_one_int <- RSSq_knot_both_int <- RSSq_knot_one_int <- -10000
        } else {
          # Find LOF measures for all basis with interactions. Need to also fit additive models (they are candidates too)!
          var_name_struct <- which(((var_name != var_name_vec) * mod_struct) == 1)  # This ensures that an interaction with the same variable is not included in the new basis.
          colnames(B)[1] <- c("")
          B2 <- as.matrix(B[, var_name_struct])
          if (k != 1 & (sum(!var_name_vec[-1]%in%var_name) > 0)) B2 <- as.matrix(B2[, -1])
          if (ncol(B2) == 0) B2 <- as.matrix(B[, 1])

          for (nn in seq(ncol(B2))) {
            B2a <- matrix(rep(B2[, nn], 2), ncol = 2)  # This is a basis function for potential parent basis. Need both part (hinges) here.
            B2b <- matrix(B2[, nn], ncol = 1)
            B_new_both_int <- cbind(B, B2a * cbind(b1_new, b2_new))
            B_new_one_int <- cbind(B, B2b * b1_new)  # Interaction model with one truncated function (the positive part).
            # calculate stats
            meas_model_both_int <- stat_out(Y = Y, B1 = B_new_both_int, TSS = TSS, GCV.null = GCV.null, pen = pen)
            meas_model_one_int <- stat_out(Y = Y, B1 = B_new_one_int, TSS = TSS, GCV.null = GCV.null, pen = pen)
            # save stats
            RSSq_knot_both_int <- c(RSSq_knot_both_int, meas_model_both_int$RSSq1)
            RSSq_knot_one_int <- c(RSSq_knot_one_int, meas_model_one_int$RSSq1)
            GCVq_knot_both_int <- c(GCVq_knot_both_int, meas_model_both_int$GCVq1)
            GCVq_knot_one_int <- c(GCVq_knot_one_int, meas_model_one_int$GCVq1)
          }

          B_new_both_add <- cbind(B, b1_new, b2_new)
          B_new_one_add <- cbind(B, b1_new)
          # calculate stats
          meas_model_both_add <- stat_out(Y = Y, B1 = B_new_both_add, TSS = TSS, GCV.null = GCV.null, pen = pen)
          meas_model_one_add <- stat_out(Y = Y, B1 = B_new_one_add, TSS = TSS, GCV.null = GCV.null, pen = pen)
          # save stats
          RSSq_knot_both_add <- c(RSSq_knot_both_add, meas_model_both_add$RSSq1)
          RSSq_knot_one_add <- c(RSSq_knot_one_add, meas_model_one_add$RSSq1)
          GCVq_knot_both_add <- c(GCVq_knot_both_add, meas_model_both_add$GCVq1)
          GCVq_knot_one_add <- c(GCVq_knot_one_add, meas_model_one_add$GCVq1)
        }

        # This next bit combines all lack-of-fit (LOF) measures that were found for each knot (t).
        RSSq_knot_both_int_mat <- rbind(RSSq_knot_both_int_mat, RSSq_knot_both_int)
        RSSq_knot_both_add_mat <- rbind(RSSq_knot_both_add_mat, RSSq_knot_both_add)
        RSSq_knot_one_int_mat <- rbind(RSSq_knot_one_int_mat, RSSq_knot_one_int)
        RSSq_knot_one_add_mat <- rbind(RSSq_knot_one_add_mat, RSSq_knot_one_add)
        GCVq_knot_both_int_mat <- rbind(GCVq_knot_both_int_mat, GCVq_knot_both_int)
        GCVq_knot_both_add_mat <- rbind(GCVq_knot_both_add_mat, GCVq_knot_both_add)
        GCVq_knot_one_int_mat <- rbind(GCVq_knot_one_int_mat, GCVq_knot_one_int)
        GCVq_knot_one_add_mat <- rbind(GCVq_knot_one_add_mat, GCVq_knot_one_add)
      } # Loop for t ends here.

      # These four conditions are used to determine whether we include interactions (first) or additive parts only to the model, Then they look further by checking
      # which truncations to include for each parent basis in the set - i.e., one truncated function or both. All must be tested/checked for FUBAR (NA) results!!!
      # We apply the parsimonious principle here in case of tied measures - i.e., choose: additive model > interaction model, and one variable model > two variable model.

      if (all((apply(RSSq_knot_both_int_mat, 1, is.na))) & all((apply(RSSq_knot_one_int_mat, 1, is.na)))) {  # Look at both truncation types (check for NA). This says interactions were FUBAR because of an NA present in both.
        int <- FALSE
        if (any(!is.na(RSSq_knot_both_add_mat)) & any(!is.na(RSSq_knot_one_add_mat))) {  # This says that an additive model for both one and two trunc. types were OK.
          if (utils::tail(max(RSSq_knot_both_add_mat, na.rm = TRUE), n = 1) > utils::tail(max(RSSq_knot_one_add_mat, na.rm = TRUE), n = 1)) {
            trunc.type <- 2
            RSSq_knot <- RSSq_knot_both_add_mat
            min_knot1 <- utils::tail(which(utils::tail(max(round(RSSq_knot, 6), na.rm = TRUE), n = 1) == round(RSSq_knot, 6), arr.ind = TRUE), n = 1)[1]
          }
          if (utils::tail(max(RSSq_knot_both_add_mat, na.rm = TRUE), n = 1) <= utils::tail(max(RSSq_knot_one_add_mat, na.rm = TRUE), n = 1)) {
            trunc.type <- 1
            RSSq_knot <- RSSq_knot_one_add_mat
            min_knot1 <- utils::tail(which(utils::tail(max(round(RSSq_knot, 6), na.rm = TRUE), n = 1) == round(RSSq_knot, 6), arr.ind = TRUE), n = 1)[1]
          }
        }
        if (all(is.na(RSSq_knot_both_add_mat)) & any(!is.na(RSSq_knot_one_add_mat))) {  # This says additive for one trunc type was FUBAR.
          trunc.type <- 1  # For example, additive and both gave FUBAR results, so use only the one (+) truncated function in the set.
          RSSq_knot <- RSSq_knot_one_add_mat
          min_knot1 <- utils::tail(which.max(round(RSSq_knot, 6)), n = 1)  # We look for max because recall: GCVq1 <- 1 - GCV1/GCV.null and RSSq1 <- 1-RSS1/TSS. Want these to be max. Tail is used in case of ties.
        }
        if (any(!is.na(RSSq_knot_both_add_mat)) & all(is.na(RSSq_knot_one_add_mat))) { # This says additive for both trunc type was OK.
          trunc.type <- 2
          RSSq_knot <- RSSq_knot_one_add_mat
          min_knot1 <- utils::tail(which.max(round(RSSq_knot, 6)), n = 1)
        }
        if (all(is.na(RSSq_knot_both_add_mat)) & all(is.na(RSSq_knot_one_add_mat))) {  # This says everything (both additive trunc type) were FUBAR.
          breakFlag <- TRUE
          break
        }
      }

      if (sum(!(apply(RSSq_knot_both_int_mat, 1, is.na))) == 0 & any(!(apply(RSSq_knot_one_int_mat, 1, is.na)))) {  # This says interactions when using both trunc. type were FUBR.
        if (all(is.na(RSSq_knot_both_add_mat))) {
          trunc.type <- 1
          if (all(is.na(RSSq_knot_one_add_mat))) {
            int <- TRUE
            RSSq_knot <- RSSq_knot_one_int_mat
            temp <- utils::tail(which(utils::tail(max(round(RSSq_knot, 6), na.rm = TRUE), n = 1) == round(RSSq_knot, 6), arr.ind = TRUE), n = 1)
            min_knot1 <- temp[1]
            best.var <- temp[2]
          }
          if (utils::tail(max(RSSq_knot_one_int_mat, na.rm = TRUE), n = 1) > utils::tail(max(RSSq_knot_one_add_mat, na.rm = TRUE), n = 1)) { # This says interactions using one trunc. type were superior over the additive model.
            int <- TRUE
            RSSq_knot <- RSSq_knot_one_int_mat
            temp <- utils::tail(which(utils::tail(max(round(RSSq_knot, 6), na.rm = TRUE), n = 1) == round(RSSq_knot, 6), arr.ind = TRUE), n = 1)
            min_knot1 <- temp[1]
            best.var <- temp[2]
          }
          if (utils::tail(max(RSSq_knot_one_int_mat, na.rm = TRUE), n = 1) <= utils::tail(max(RSSq_knot_one_add_mat, na.rm = TRUE), n = 1)) {
            int <- FALSE
            RSSq_knot <- RSSq_knot_one_add_mat
            min_knot1 <- utils::tail(which.max(round(RSSq_knot, 6)), n = 1)
          }
        }
        if (any(!is.na(RSSq_knot_both_add_mat))) {
          if (utils::tail(max(RSSq_knot_one_int_mat, na.rm = TRUE), n = 1) > utils::tail(max(RSSq_knot_both_add_mat, na.rm = TRUE), n = 1)) {
            trunc.type <- 1
            if (utils::tail(max(RSSq_knot_one_int_mat, na.rm = TRUE), n = 1) > utils::tail(max(RSSq_knot_one_add_mat, na.rm = TRUE), n = 1)) {
              int <- TRUE
              RSSq_knot <- RSSq_knot_one_int_mat
              temp <- utils::tail(which(utils::tail(max(round(RSSq_knot, 6), na.rm = TRUE), n = 1) == round(RSSq_knot, 6), arr.ind = TRUE), n = 1)
              min_knot1 <- temp[1]
              best.var <- temp[2]
            }
            if (utils::tail(max(RSSq_knot_one_int_mat, na.rm = TRUE)) <= utils::tail(max(RSSq_knot_one_add_mat, na.rm = TRUE))) {
              int <- FALSE
              RSSq_knot <- RSSq_knot_one_add_mat
              min_knot1 <- utils::tail(which.max(round(RSSq_knot, 6)), n = 1)
            }
          }
          if (utils::tail(max(RSSq_knot_one_int_mat, na.rm = TRUE), n = 1) <= utils::tail(max(RSSq_knot_both_add_mat, na.rm = TRUE), n = 1)) {
            int <- FALSE
            if (utils::tail(max(RSSq_knot_both_add_mat, na.rm = TRUE), n = 1) > utils::tail(max(RSSq_knot_one_add_mat, na.rm = TRUE), n = 1)) {
              trunc.type <- 2
              RSSq_knot <- RSSq_knot_both_add_mat
              min_knot1 <- utils::tail(which(utils::tail(max(round(RSSq_knot, 6), na.rm = TRUE), n = 1) == round(RSSq_knot, 6), arr.ind = TRUE), n = 1)[1]
            }
            if (utils::tail(max(RSSq_knot_both_add_mat, na.rm = TRUE), n = 1) <= utils::tail(max(RSSq_knot_one_add_mat, na.rm = TRUE), n = 1)) {
              trunc.type <- 1
              RSSq_knot <- RSSq_knot_one_add_mat
              min_knot1 <- utils::tail(which(utils::tail(max(round(RSSq_knot, 6), na.rm = TRUE), n = 1) == round(RSSq_knot, 6), arr.ind = TRUE), n = 1)[1]
            }
          }
        }
      }

      if (any(!(apply(RSSq_knot_both_int_mat, 1, is.na))) & all((apply(RSSq_knot_one_int_mat, 1, is.na)))) { # This says interactions for one trunc. type was FUBAR.
        if (all(is.na(RSSq_knot_both_add_mat))) {  # Check if additive for only one trunc type was OK.
          if (all(is.na(RSSq_knot_one_add_mat))) {
            int <- TRUE
            trunc.type <- 2
            RSSq_knot <- RSSq_knot_both_int_mat
            temp <- utils::tail(which(utils::tail(max(round(RSSq_knot, 6), na.rm = TRUE), n = 1) == round(RSSq_knot, 6), arr.ind = TRUE), n = 1)
            min_knot1 <- temp[1]
            best.var <- temp[2]
          }
          if (utils::tail(max(RSSq_knot_both_int_mat, na.rm = TRUE), n = 1) > utils::tail(max(RSSq_knot_one_add_mat, na.rm = TRUE), n = 1)) {  # Compare interaction (with both) with additive (with one).
            int <- TRUE
            trunc.type <- 2
            RSSq_knot <- RSSq_knot_both_int_mat
            temp <- utils::tail(which(utils::tail(max(round(RSSq_knot, 6), na.rm = TRUE), n = 1) == round(RSSq_knot, 6), arr.ind = TRUE), n = 1)
            min_knot1 <- temp[1]
            best.var <- temp[2]
          }
          if (utils::tail(max(RSSq_knot_both_int_mat, na.rm = TRUE), n = 1) <= utils::tail(max(RSSq_knot_one_add_mat, na.rm = TRUE), n = 1)) {
            int <- FALSE
            trunc.type <- 1
            RSSq_knot <- RSSq_knot_one_add_mat
            min_knot1 <- utils::tail(which.max(round(RSSq_knot, 6)), n = 1)
          }
        }
        if (any(!is.na(RSSq_knot_both_add_mat))) {
          if (utils::tail(max(RSSq_knot_both_int_mat, na.rm = TRUE), n = 1) > utils::tail(max(RSSq_knot_both_add_mat, na.rm = TRUE), n = 1)) {
            trunc.type <- 1
            if (utils::tail(max(RSSq_knot_both_int_mat, na.rm = TRUE), n = 1) > utils::tail(max(RSSq_knot_one_add_mat, na.rm = TRUE), n = 1)) {
              int <- TRUE
              RSSq_knot <- RSSq_knot_both_int_mat
              temp <- utils::tail(which(utils::tail(max(round(RSSq_knot, 6), na.rm = TRUE), n = 1) == round(RSSq_knot, 6), arr.ind = TRUE), n = 1)
              min_knot1 <- temp[1]
              best.var <- temp[2]
            }
            if (utils::tail(max(RSSq_knot_both_int_mat, na.rm = TRUE)) <= utils::tail(max(RSSq_knot_one_add_mat, na.rm = TRUE))) {
              int <- FALSE
              RSSq_knot <- RSSq_knot_one_add_mat
              min_knot1 <- utils::tail(which.max(round(RSSq_knot, 6)), n = 1)
            }
          }
          if (utils::tail(max(RSSq_knot_both_int_mat, na.rm = TRUE), n = 1) <= utils::tail(max(RSSq_knot_both_add_mat, na.rm = TRUE), n = 1)) {
            int <- FALSE
            if (utils::tail(max(RSSq_knot_both_add_mat, na.rm = TRUE), n = 1) > utils::tail(max(RSSq_knot_one_add_mat, na.rm = TRUE), n = 1)) {
              trunc.type <- 2
              RSSq_knot <- RSSq_knot_both_add_mat
              min_knot1 <- utils::tail(which.max(round(RSSq_knot, 6)), n = 1)
            }
            if (utils::tail(max(RSSq_knot_both_add_mat, na.rm = TRUE), n = 1) <= utils::tail(max(RSSq_knot_one_add_mat, na.rm = TRUE), n = 1)) {
              trunc.type <- 1
              RSSq_knot <- RSSq_knot_one_add_mat
              min_knot1 <- utils::tail(which.max(round(RSSq_knot, 6)), n = 1)
            }
          }
        }
      }

      if (any(!(apply(RSSq_knot_both_int_mat, 1, is.na))) & any(!(apply(RSSq_knot_one_int_mat, 1, is.na)))) { # This says interactions for one trunc. and both types were OK.
        if (utils::tail(max(RSSq_knot_both_int_mat, na.rm = TRUE), n = 1) > utils::tail(max(RSSq_knot_one_int_mat, na.rm = TRUE), n = 1)) {
          if (any(!is.na(RSSq_knot_both_add_mat))) {
            if (utils::tail(max(RSSq_knot_both_add_mat, na.rm = TRUE), n = 1) >= utils::tail(max(RSSq_knot_both_int_mat, na.rm = TRUE), n = 1)) {
              int <- FALSE
              if (utils::tail(max(RSSq_knot_both_add_mat, na.rm = TRUE), n = 1) > utils::tail(max(RSSq_knot_one_add_mat, na.rm = TRUE), n = 1)) {
                trunc.type <- 2
                RSSq_knot <- RSSq_knot_both_add_mat
                min_knot1 <- utils::tail(which.max(round(RSSq_knot, 6)), n = 1)
              }
              if (utils::tail(max(RSSq_knot_both_add_mat, na.rm = TRUE), n = 1) <= utils::tail(max(RSSq_knot_one_add_mat, na.rm = TRUE), n = 1)) {
                trunc.type <- 1
                RSSq_knot <- RSSq_knot_one_add_mat
                min_knot1 <- utils::tail(which.max(round(RSSq_knot, 6)))
              }
            }
            if (utils::tail(max(RSSq_knot_both_add_mat, na.rm = TRUE), n = 1) < utils::tail(max(RSSq_knot_both_int_mat, na.rm = TRUE), n = 1)) {
              if (utils::tail(max(RSSq_knot_both_int_mat, na.rm = TRUE), n = 1) > utils::tail(max(RSSq_knot_one_add_mat, na.rm = TRUE), n = 1)) {
                int <- TRUE
                trunc.type <- 2
                RSSq_knot <- RSSq_knot_both_int_mat
                temp <- utils::tail(which(utils::tail(max(round(RSSq_knot, 6), na.rm = TRUE), n = 1) == round(RSSq_knot, 6), arr.ind = TRUE), n = 1)
                min_knot1 <- temp[1]
                best.var <- temp[2]
              }
              if (utils::tail(max(RSSq_knot_both_int_mat, na.rm = TRUE), n = 1) <= utils::tail(max(RSSq_knot_one_add_mat, na.rm = TRUE), n = 1)) {
                int <- FALSE
                trunc.type <- 1
                RSSq_knot <- RSSq_knot_one_add_mat
                min_knot1 <- utils::tail(which.max(round(RSSq_knot, 6)), n = 1)
              }
            }
          }
          if (all(is.na(RSSq_knot_both_add_mat))) {
            if (all(is.na(RSSq_knot_one_add_mat))) {
              int <- TRUE
              trunc.type <- 2
              RSSq_knot <- RSSq_knot_both_int_mat
              temp <- utils::tail(which(utils::tail(max(round(RSSq_knot, 6), na.rm = TRUE), n = 1) == round(RSSq_knot, 6), arr.ind = TRUE), n = 1)
              min_knot1 <- temp[1]
              best.var <- temp[2]
            }
            if (utils::tail(max(RSSq_knot_one_add_mat, na.rm = TRUE), n = 1) >= utils::tail(max(RSSq_knot_both_int_mat, na.rm = TRUE), n = 1)) {
              int <- FALSE
              trunc.type <- 1
              RSSq_knot <- RSSq_knot_one_add_mat
              min_knot1 <- utils::tail(which.max(round(RSSq_knot, 6)), n = 1)
            }
            if (utils::tail(max(RSSq_knot_one_add_mat, na.rm = TRUE), n = 1) < utils::tail(max(RSSq_knot_both_int_mat, na.rm = TRUE), n = 1)) {
              int <- TRUE
              trunc.type <- 2
              RSSq_knot <- RSSq_knot_both_int_mat
              temp <- utils::tail(which(utils::tail(max(round(RSSq_knot, 6), na.rm = TRUE), n = 1) == round(RSSq_knot, 6), arr.ind = TRUE), n = 1)
              min_knot1 <- temp[1]
              best.var <- temp[2]
            }
          }
        }
        if (utils::tail(max(RSSq_knot_both_int_mat, na.rm = TRUE), n = 1) <= utils::tail(max(RSSq_knot_one_int_mat, na.rm = TRUE), n = 1)) {
          if (any(!is.na(RSSq_knot_both_add_mat))) {
            if (utils::tail(max(RSSq_knot_both_add_mat, na.rm = TRUE), n = 1) >= utils::tail(max(RSSq_knot_one_int_mat, na.rm = TRUE), n = 1)) {
              int <- FALSE
              if (utils::tail(max(RSSq_knot_both_add_mat, na.rm = TRUE), n = 1) > utils::tail(max(RSSq_knot_one_add_mat, na.rm = TRUE), n = 1)) {
                trunc.type <- 2
                RSSq_knot <- RSSq_knot_both_add_mat
                min_knot1 <- utils::tail(which.max(round(RSSq_knot, 6)), n = 1)
              }
              if (utils::tail(max(RSSq_knot_both_add_mat, na.rm = TRUE), n = 1) <= utils::tail(max(RSSq_knot_one_add_mat, na.rm = TRUE), n = 1)) {
                trunc.type <- 1
                RSSq_knot <- RSSq_knot_one_add_mat
                min_knot1 <- utils::tail(which.max(round(RSSq_knot, 6)), n = 1)
              }
            }
            if (utils::tail(max(RSSq_knot_both_add_mat, na.rm = TRUE), n = 1) < utils::tail(max(RSSq_knot_one_int_mat, na.rm = TRUE), n = 1)) {
              trunc.type <- 1
              if (utils::tail(max(RSSq_knot_one_int_mat, na.rm = TRUE), n = 1) > utils::tail(max(RSSq_knot_one_add_mat, na.rm = TRUE), n = 1)) {
                int <- TRUE
                RSSq_knot <- RSSq_knot_one_int_mat
                temp <- utils::tail(which(utils::tail(max(round(RSSq_knot, 6), na.rm = TRUE), n = 1) == round(RSSq_knot, 6), arr.ind = TRUE), n = 1)
                min_knot1 <- temp[1]
                best.var <- temp[2]
              }
              if (utils::tail(max(RSSq_knot_one_int_mat, na.rm = TRUE), n = 1) <= utils::tail(max(RSSq_knot_one_add_mat, na.rm = TRUE), n = 1)) {
                int <- FALSE
                RSSq_knot <- RSSq_knot_one_add_mat
                min_knot1 <- utils::tail(which.max(round(RSSq_knot, 6)), n = 1)
              }
            }
          }
          if (all(is.na(RSSq_knot_both_add_mat))) {
            if (all(is.na(RSSq_knot_one_add_mat))) {
              int <- TRUE
              trunc.type <- 1
              RSSq_knot <- RSSq_knot_one_int_mat
              temp <- utils::tail(which(utils::tail(max(round(RSSq_knot, 6), na.rm = TRUE), n = 1) == round(RSSq_knot, 6), arr.ind = TRUE), n = 1)
              min_knot1 <- temp[1]
              best.var <- temp[2]
            }
            if (any(!is.na(RSSq_knot_one_add_mat))) {
              trunc.type <- 1
              if (utils::tail(max(RSSq_knot_one_add_mat, na.rm = TRUE), n = 1) >= utils::tail(max(RSSq_knot_one_int_mat, na.rm = TRUE), n = 1)) {
                int <- FALSE
                RSSq_knot <- RSSq_knot_one_add_mat
                min_knot1 <- utils::tail(which.max(round(RSSq_knot, 6)))
              }
              if (utils::tail(max(RSSq_knot_one_add_mat, na.rm = TRUE), n = 1) < utils::tail(max(RSSq_knot_one_int_mat, na.rm = TRUE), n = 1)) {
                int <- TRUE
                RSSq_knot <- RSSq_knot_one_int_mat
                temp <- utils::tail(which(utils::tail(max(round(RSSq_knot, 6), na.rm = TRUE), n = 1) == round(RSSq_knot, 6), arr.ind = TRUE), n = 1)
                min_knot1 <- temp[1]
                best.var <- temp[2]
              }
            }
          }
        }
      }

      b1_new <- matrix(tp1(x = X, t = X_red[min_knot1]), ncol = 1)
      b2_new <- matrix(tp2(x = X, t = X_red[min_knot1]), ncol = 1)
      colnames(b1_new) <- var_name
      colnames(b2_new) <- var_name

      B_name1 <- paste("(", var_name, "-", signif(X_red[min_knot1], 4), ")", sep = "")
      B_name2 <- paste("(", signif(X_red[min_knot1], 4), "-", var_name, ")", sep = "")

      if (int) {
        mod_struct1 <- which(mod_struct == 1)
        colnames(B)[1] <- c("")
        var_name1 <- which(var_name_vec != var_name)
        if (int.count == 0) {
          var_name2 <- var_name_vec[var_name1]
        } else if (int.count > 0) {
          var_name2 <- var_name_vec[var_name_struct]
        }
        var_name_struct <- mod_struct1[mod_struct1 %in% var_name1]
        B2 <- as.matrix(B[, var_name_struct])
        B3_names <- B_names_vec[var_name_struct]
        B3_names <- B3_names[-1]
        B2 <- as.matrix(B2[, -1])
        var_name2 <- var_name2[-1]

        if (trunc.type == 2) {
          B2a <- matrix(rep(B2[, best.var], 2), ncol = 2)
          B_temp <- cbind(B, B2a * cbind(b1_new, b2_new))
          B_new <- B2a * cbind(b1_new, b2_new)
          var_name3 <- var_name2[best.var]
          colnames(B_new) <- rep(var_name3, 2)
        } else if (trunc.type == 1) {
          B2b <- matrix(B2[, best.var], ncol = 1)
          B_temp <- cbind(B, B2b * b1_new)  # Interaction model with one truncated function (i.e., the positive part).
          B_new <- B2b * b1_new
          colnames(B_new) <- rep(var_name2[1], 1)
          var_name3 <- var_name2[best.var]
          colnames(B_new) <- rep(var_name3, 1)
        }

        B_names <- paste(B3_names[best.var], B_name1, sep = "*")
        if (trunc.type == 2) {
          B_names <- c(B_names, paste(B3_names[best.var], B_name2, sep = "*"))
        } else if (trunc.type == 1) {
          B_names <- B_names
        }

        var_name_list1 <- list()
        for (ll in seq(ncol(B_new))) {
          colnames(B_new)[ll] <- paste(var_name, colnames(B_new)[ll], sep = ":")
          var_name_list1 <- c(var_name_list1, list(colnames(B_new)[ll]))
          int.count1 <- int.count1 + 1
        }

        pp <- ncol(B_temp)
        if (trunc.type == 2) {
          colnames(B_temp)[((pp - 1):pp)] <- var_name_list1[[1]]
        } else if (trunc.type == 1) {
          colnames(B_temp)[pp] <- var_name_list1[[1]]
        }
      } else {
        var_name_list1 <- list()
        if (trunc.type == 2) {
          B_temp <- cbind(B, b1_new, b2_new)  # Additive model with both truncated functions.
          B_new <- cbind(b1_new, b2_new)
          B_names <- c(B_name1, B_name2)
          var_name_list1 <- c(var_name_list1, list(var_name))
          var_name_list1 <- c(var_name_list1, list(var_name))  # Repeat it because there are new basis function selected.
        } else if (trunc.type == 1) {
          B_temp <- cbind(B, b1_new)  # Additive model with one truncated function (positive part).
          B_new <- b1_new
          B_names <- B_name1
          var_name_list1 <- c(var_name_list1, list(var_name))
        }
      }

      meas_model <- stat_out(Y = Y, B1 = B_temp, TSS = TSS, GCV.null = GCV.null, pen = pen)
      GCVq2 <- meas_model$GCVq1
      RSSq2 <- meas_model$RSSq1

      if (GCVq2 < (-10) | (round(RSSq2, 4) > (1 - tols))) {
        # writeLines("MARS Tolerance criteria met 1. \n")
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

        RSSq_term_temp <- c(RSSq_term_temp, NA)
        GCVq_term_temp <- c(GCVq_term_temp, NA)

        if (length(var.mod_temp) == q) {
          breakFlag <- TRUE
          break
        } else {
          next
        }
      } else if (GCVq2 >= (-10) | (round(RSSq2, 4) <= (1 - tols))) {
        RSSq_term_temp <- c(RSSq_term_temp, RSSq2)
        GCVq_term_temp <- c(GCVq_term_temp, GCVq2)
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

    } # End the for () loop here.

    if (breakFlag) break

    best.mod <- which.max(RSSq_term_temp)  # Finds the best model (max LOF) from your candidate model/basis set. This will become the new parent.

    RSSq2 <- RSSq_term_temp[best.mod]
    GCVq2 <- GCVq_term_temp[best.mod]
    min_knot_vec1 <- min_knot_vec_temp[best.mod]
    int.count1 <- int.count1_temp[best.mod]
    int <- is.int_temp[best.mod]
    trunc.type <- trunc.type_temp[best.mod]

    B_new <- B_new_list_temp[[best.mod]]
    var_name_list1 <- var_name_list1_temp[[best.mod]]
    B_names <- B_names_temp[[best.mod]]
    X_red <- X_red_temp[[best.mod]]
    B_temp <- B_temp_list[[best.mod]]

    RSSq_term <- c(RSSq_term, RSSq2)
    GCVq_term <- c(GCVq_term, GCVq2)
    min_knot_vec <- c(min_knot_vec, min_knot_vec1)
    pred.name_vec <- c(pred.name_vec, colnames(B_new)[1])
    cut_vec <- c(cut_vec, X_red[min_knot_vec1])
    trunc.type_vec <- c(trunc.type_vec, trunc.type)
    is.int_vec <- c(is.int_vec, int)

    if (abs(RSSq_term[k] - RSSq_term[k + 1]) < tols) {
      breakFlag <- TRUE
      break
    }

    if (abs(RSSq_term[k] - RSSq_term[k + 1]) >= tols) {
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

    if (nrow(B) <= (ncol(B) + 2)) { # To avoid p>N issues!
      ok <- FALSE
    }

    if (m >= M) {   # If model exceeds no. of set terms, terminate it.
      ok <- FALSE
    }

    stat_out(Y = Y, B1 = B, TSS = TSS, GCV.null = GCV.null, pen = pen)
  }  # The term (m) for () loop ends here.

  colnames(B) <- B_names_vec
  ncol_B <- ncol(B)

  B2 <- B # Set this for final model to be the same as model selected from the forward pass.

  # Algorithm 3 (backward pass), as in Friedman (1991). This uses GCV as the selection criterion.

  p <- ncol_B + pen * (ncol_B - 1) / 2  # This matches the earth() package, SAS and Friedman (1991) penalty.

  full_RSS <- sum((Y - stats::fitted(stats::lm.fit(B - 1, Y)))^2)
  full_GCV <- full_RSS / (N * (1 - p / N)^2)
  full_GCV <- 1 - (full_GCV / GCV.null)

  B_new <- B
  GCV_mat <- matrix(NA, ncol = ncol_B, nrow = ncol_B)
  colnames(GCV_mat) <- colnames(B)
  GCV_mat <- cbind(GCV_mat, rep(NA, ncol_B))
  colnames(GCV_mat)[(ncol_B + 1)] <- "Forward pass model"

  GCV_mat[1, (ncol_B + 1)] <- full_GCV
  GCV1 <- backward_sel(Y = Y, B_new = B_new, pen = pen, GCV.null = GCV.null)
  GCV_mat[2, 2:(length(GCV1) + 1)] <- GCV1
  variable.lowest <- utils::tail(which(GCV1 == max(GCV1, na.rm = TRUE)), n = 1)
  var.low.vec <- c(colnames(B_new)[variable.lowest + 1])
  B_new <- as.matrix(B_new[, -(variable.lowest + 1)])

  for (i in 2:(ncol_B - 1)) {
    GCV1 <- backward_sel(Y = Y, B_new = B_new, pen = pen, GCV.null = GCV.null)
    variable.lowest <- utils::tail(which(GCV1 == max(GCV1, na.rm = TRUE)), n = 1)
    var.low.vec <- c(var.low.vec, colnames(B_new)[variable.lowest + 1])
    if (i != (ncol_B - 1)) {
      GCV_mat[(i + 1), colnames(B_new)[-1]] <- GCV1
    } else {
      GCV_mat[(i + 1), 1] <- GCV1
    }
    B_new <- as.matrix(B_new[, -(variable.lowest + 1)])
  }

  number.rm1 <- which(max(GCV_mat, na.rm = TRUE) == GCV_mat, arr.ind = TRUE)[1] - 1

  if (number.rm1 == 0) {
    # if (print.disp == TRUE) writeLines("\n ** Forward pass model was chosen after pruning/backward selection for MARS** \n")
    B_final <- B
  } else if (number.rm1 == (nrow(GCV_mat) - 1)) {
    B_final <- as.matrix(B[, 1])
  }

  if (!(number.rm1 == 0) & !(number.rm1 == (nrow(GCV_mat) - 1))) {
    number.rm2 <- rep(NA, length(var.low.vec_red))
    var.low.vec_red <- var.low.vec[1:number.rm1]
    for (k in seq(length(var.low.vec_red))) {
      number.rm2[k] <- which(var.low.vec_red[k] == colnames(B))
    }
    B_final <- B[, -c(number.rm2)]
  }

  if (ncol(B_final) == 1) colnames(B_final) <- "Intercept"

  # Some measures to be reported, GCV, RSq, etc. for the final model fit.
  est.many <- mvabund::manyglm(Y ~ B_final - 1,
                               family = "negative.binomial",
                               maxiter = 1000,
                               maxiter2 = 100)
  final_mod <- MASS::glm.nb(c(t(Y)) ~ B_final - 1,
                            method = "glm.fit2",
                            init.theta = est.many$theta)

  final_stats <- stat_out(Y = Y, B1 = B_final, TSS = TSS, GCV.null = GCV.null, pen = pen)

  res <- NULL
  res$bx <- B_final
  res$GCV_mat <- GCV_mat
  res$min_GCV_own <- final_stats$GCV1
  res$y_pred <- stats::predict(final_mod)
  res$final_mod <- final_mod

  class(z) <- "marge"
  return(z)
}
