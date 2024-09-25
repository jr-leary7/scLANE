#' Build an NB GLMM using truncated power basis functions.
#'
#' @name fitGLMM
#' @author Jack Leary
#' @description Fits a negative binomial generalized linear mixed model using truncated power basis function splines as input. The basis matrix can be created adaptively using subject-specific estimation of optimal knots using \code{\link{marge2}}, or basis functions can be evenly space across quantiles. The resulting model can output subject-specific and population-level fitted values.
#' @importFrom purrr map_dfr pmap_dfc
#' @importFrom dplyr mutate select if_else
#' @importFrom stats quantile logLik fitted as.formula
#' @importFrom glmmTMB glmmTMB nbinom2
#' @param X_pred A matrix with one column containing cell ordering values. Defaults to NULL.
#' @param Y A vector of raw single cell counts. Defaults to NULL.
#' @param Y.offset (Optional) An offset to be included in the final model fit. Defaults to NULL.
#' @param id.vec A vector of subject IDs. Defaults to NULL.
#' @param adaptive Should basis functions be chosen adaptively? Defaults to TRUE.
#' @param approx.knot Should knot approximation be used in the calls to \code{\link{marge2}}? This speeds up computation somewhat. Defaults to TRUE.
#' @param M.glm The number of possible basis functions to use in the calls to \code{\link{marge2}} when choosing basis functions adaptively.
#' @param return.basis (Optional) Whether the basis model matrix (denoted \code{B_final}) should be returned as part of the \code{marge} model object. Defaults to FALSE.
#' @param return.GCV (Optional) Whether the final GCV value should be returned as part of the \code{marge} model object. Defaults to FALSE.
#' @param verbose (Optional) Should intermediate output be printed to the console? Defaults to FALSE.
#' @return An object of class \code{marge} containing the fitted model & other optional quantities of interest (basis function matrix, GCV, etc.).
#' @seealso \code{\link[glmmTMB]{glmmTMB}}
#' @seealso \code{\link{testDynamic}}
#' @seealso \code{\link{modelLRT}}
#' @export
#' @examples
#' data(sim_counts)
#' data(sim_pseudotime)
#' cell_offset <- createCellOffset(sim_counts)
#' glmm_mod <- fitGLMM(X_pred = sim_pseudotime,
#'                     Y = BiocGenerics::counts(sim_counts)[4, ],
#'                     Y.offset = cell_offset,
#'                     id.vec = sim_counts$subject)

fitGLMM <- function(X_pred = NULL,
                    Y = NULL,
                    Y.offset = NULL,
                    id.vec = NULL,
                    adaptive = TRUE,
                    approx.knot = TRUE,
                    M.glm = 3,
                    return.basis = FALSE,
                    return.GCV = FALSE,
                    verbose = FALSE) {
  # check inputs
  if (is.null(X_pred) || is.null(Y) || is.null(id.vec)) { stop("You forgot some inputs to fitGLMM().") }
  if (is.unsorted(id.vec)) { stop("Your data must be ordered by subject, please do so before running fitGLMM().") }
  pen <- 2  # GCV criterion penalty -- matches marge2()
  # fit NB GLMM
  if (adaptive) {
      glm_marge_knots <- purrr::map(unique(id.vec),
                                          function(x) {
                                            marge_mod_sub <- marge2(X_pred = X_pred[which(id.vec == x), , drop = FALSE],
                                                                    Y = Y[which(id.vec == x)],
                                                                    Y.offset = Y.offset[which(id.vec == x)],
                                                                    approx.knot = approx.knot,
                                                                    M = M.glm,
                                                                    return.basis = TRUE)
                                            })
        keepmods <- which(sapply(glm_marge_knots, function(x) length(x$coef_names)) > 1)
        glm_marge_knots <- glm_marge_knots[keepmods]

        allKnots <- lapply(glm_marge_knots, function(x) extractBreakpoints(x)$Breakpoint)
        allCoef <- lapply(glm_marge_knots, function(x) names(coef(x$final_mod)[-1]))
        allOldCoef <- lapply(glm_marge_knots, function(x) paste0("B_final", x$marge_coef_names[-1]))
        tp_fun <- lapply(allCoef, function(x) dplyr::if_else(grepl("h_[0-9]", x), "tp2", "tp1"))

        glm_marge_knots <- data.frame(knot = do.call(c, allKnots),
                                     coef = do.call(c, allCoef),
                                     old_coef = do.call(c, allOldCoef),
                                     tp_fun = do.call(c, tp_fun))

        glmm_basis_df <- purrr::pmap_dfc(list(glm_marge_knots$knot,
                                          glm_marge_knots$tp_fun,
                                          seq_len(nrow(glm_marge_knots))),
                                     function(k, f, i) {
                                       if (f == "tp1") {
                                         basis <- tp1(x = X_pred[, 1], t = k)
                                       } else {
                                         basis <- tp2(x = X_pred[, 1], t = k)
                                       }
                                       basis_df <- data.frame(basis)
                                       colnames(basis_df) <- paste0("X", i)
                                       return(basis_df)
                                     })
        if (ncol(glmm_basis_df) == 1) {
            mod_formula <- stats::as.formula(paste0("Y ~ ",
                                                    paste(colnames(glmm_basis_df), collapse = " + "),
                                                    " + (1 + ", paste(colnames(glmm_basis_df), collapse = " + "),
                                                    " | subject)"))
            glmm_basis_df_new <- dplyr::mutate(glmm_basis_df,
                                               Y = Y,
                                               subject = id.vec,
                                               .before = 1)
            nonzero_coefs <- 1
        } else {
            marge_style_names <- glm_marge_knots$old_coef
            coef_names <- glm_marge_knots$coef
            if (verbose) {
              message(paste0("Generated a total of ",
                             length(coef_names),
                             " basis functions across ",
                             length(unique(id.vec)),
                             " subjects."))
            }
            # run NB LASSO with all possible basis functions
            lasso_formula <- stats::as.formula(paste0("Y ~ ", paste0(colnames(glmm_basis_df), collapse = " + ")))
            if (is.null(Y.offset)) {
              pruned_model <- mpath::glmregNB(lasso_formula,
                                              data = glmm_basis_df,
                                              parallel = FALSE,
                                              nlambda = 50,
                                              penalty = "snet",
                                              alpha = .5,
                                              standardize = TRUE,
                                              trace = FALSE,
                                              maxit.theta = 1,
                                              link = log)
            } else {
              pruned_model <- mpath::glmregNB(lasso_formula,
                                              data = glmm_basis_df,
                                              offset = log(1 / Y.offset),
                                              parallel = FALSE,
                                              nlambda = 50,
                                              penalty = "snet",
                                              alpha = .5,
                                              standardize = TRUE,
                                              trace = FALSE,
                                              maxit.theta = 1,
                                              link = log)
            }
            # identify nonzero basis functions in minimum AIC model
            nonzero_coefs <- which(as.numeric(pruned_model$beta[, which.min(pruned_model$aic)]) != 0)
            # build formula automatically
            mod_formula <- stats::as.formula(paste0("Y ~ ",
                                                    paste(colnames(glmm_basis_df)[nonzero_coefs], collapse = " + "),
                                                    " + (1 + ", paste(colnames(glmm_basis_df)[nonzero_coefs], collapse = " + "),
                                                    " | subject)"))
            glmm_basis_df_new <- dplyr::mutate(glmm_basis_df,
                                           Y = Y,
                                           subject = id.vec,
                                           .before = 1)
       }
    if (is.null(Y.offset)) {
      glmm_mod <- glmmTMB::glmmTMB(mod_formula,
                                   data = glmm_basis_df_new,
                                   family = glmmTMB::nbinom2(link = "log"),
                                   se = TRUE,
                                   REML = FALSE)
    } else {
      glmm_mod <- glmmTMB::glmmTMB(mod_formula,
                                   data = glmm_basis_df_new,
                                   offset = log(1 / Y.offset),
                                   family = glmmTMB::nbinom2(link = "log"),
                                   se = TRUE,
                                   REML = FALSE)
    }
  } else {
    glmm_basis_df_new <- data.frame(X1 = tp1(X_pred[, 1], t = round(as.numeric(stats::quantile(X_pred[, 1], 1/3)), 4)),
                                    X2 = tp1(X_pred[, 1], t = round(as.numeric(stats::quantile(X_pred[, 1], 2/3)), 4)),
                                    X3 = tp2(X_pred[, 1], t = round(as.numeric(stats::quantile(X_pred[, 1], 1/3)), 4)),
                                    X4 = tp2(X_pred[, 1], t = round(as.numeric(stats::quantile(X_pred[, 1], 2/3)), 4)),
                                    Y = Y,
                                    subject = id.vec)
    marge_style_names <- c(paste0("B_final(", round(as.numeric(stats::quantile(X_pred[, 1], 1/3)), 4), "-PT)"),
                           paste0("B_final(PT-", round(as.numeric(stats::quantile(X_pred[, 1], 1/3)), 4), ")"),
                           paste0("B_final(", round(as.numeric(stats::quantile(X_pred[, 1], 2/3)), 4), "-PT)"),
                           paste0("B_final(PT-", round(as.numeric(stats::quantile(X_pred[, 1], 2/3)), 4), ")"))
    coef_names <- c(paste0("h_", round(as.numeric(stats::quantile(X_pred[, 1], 1/3)), 4), "_PT"),
                    paste0("h_PT_", round(as.numeric(stats::quantile(X_pred[, 1], 1/3)), 4)),
                    paste0("h_", round(as.numeric(stats::quantile(X_pred[, 1], 2/3)), 4), "_PT"),
                    paste0("h_PT_", round(as.numeric(stats::quantile(X_pred[, 1], 2/3)), 4)))
    nonzero_coefs <- seq(4)
    if (is.null(Y.offset)) {
      glmm_mod <- glmmTMB::glmmTMB(Y ~ X1 + X2 + X3 + X4 + (1 + X1 + X2 + X3 + X4 | subject),
                                   data = glmm_basis_df_new,
                                   family = glmmTMB::nbinom2(link = "log"),
                                   se = TRUE,
                                   REML = FALSE)
    } else {
      glmm_mod <- glmmTMB::glmmTMB(Y ~ X1 + X2 + X3 + X4 + (1 + X1 + X2 + X3 + X4 | subject),
                                   data = glmm_basis_df_new,
                                   offset = log(1 / Y.offset),
                                   family = glmmTMB::nbinom2(link = "log"),
                                   se = TRUE,
                                   REML = FALSE)
    }
  }
  # set up results
  marge_style_names <- c("B_finalIntercept", marge_style_names[nonzero_coefs])
  coef_names <- c("Intercept", coef_names[nonzero_coefs])
  res <- list(final_mod = glmm_mod,
              basis_mtx = NULL,
              WIC_mtx = NULL,  # only used in GLM / GEE models
              GCV = NULL,
              model_type = "GLMM",
              coef_names = coef_names,
              marge_coef_names = marge_style_names)
  if (return.basis) {
    res$basis_mtx <- dplyr::select(glmm_basis_df_new, -c(subject, Y))
  }
  if (return.GCV) {
    p_2 <- attributes(stats::logLik(glmm_mod))$df  # degrees of freedom, differs from how it's done in GLM / GEE marge2() -- might need to change to approximation method later idk
    df1a <- p_2 + pen * (p_2 - 1) / 2
    RSS1 <- sum((Y - stats::fitted(glmm_mod))^2)
    res$GCV <- RSS1 / (length(Y) * (1 - df1a / length(Y))^2)  # Note there are some problems w/ GCV in correlated data, so use this as a metric of GoF w/ caution
  }
  class(res) <- "marge"
  return(res)
}
