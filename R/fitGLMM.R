#' Build an NB GLMM using truncated power basis functions.
#'
#' @name fitGLMM
#' @author Jack Leary
#' @description Fits a negative binomial generalized linear mixed model using truncated power basis function splines as input. The basis matrix can be created adaptively using subject-specific estimation of optimal knots using \code{\link{marge2}}, or basis functions can be evenly space across quantiles. The resulting model can output subject-specific and population-level fitted values.
#' @importFrom purrr map_dfr pmap_dfc
#' @importFrom dplyr mutate select
#' @importFrom stats quantile logLik fitted
#' @importFrom glmmTMB glmmTMB nbinom1
#' @param X_pred A matrix with one column containing cell ordering values. Defaults to NULL.
#' @param Y A vector of raw single cell counts. Defaults to NULL.
#' @param id.vec A vector of subject IDs. Defaults to NULL.
#' @param adaptive Should basis functions be chosen adaptively? Defaults to FALSE.
#' @param approx.knot Should knot approximation be used in the calls to \code{\link{marge2}}? This speeds up computation somewhat. Defaults to TRUE.
#' @param M.glm The number of possible basis functions to use in the calls to \code{\link{marge2}} when choosing basis functions adaptively.
#' @param return.basis (Optional) Whether the basis model matrix (denoted \code{B_final}) should be returned as part of the \code{marge} model object. Defaults to FALSE.
#' @param return.GCV (Optional) Whether the final GCV value should be returned as part of the \code{marge} model object. Defaults to FALSE.
#' @return An object of class \code{marge} containing the fitted model & other optional quantities of interest (basis function matrix, GCV, etc.).
#' @seealso \code{\link[glmmTMB]{glmmTMB}}
#' @seealso \code{\link{testDynamic}}
#' @seealso \code{\link{modelLRT}}
#' @export
#' @examples
#' \dontrun{
#' fitGLMM(X_pred = pt_df,
#'         Y = raw_counts,
#'         id.vec = subject_vec)
#' }

fitGLMM <- function(X_pred = NULL,
                    Y = NULL,
                    id.vec = NULL,
                    adaptive = FALSE,
                    approx.knot = TRUE,
                    M.glm = 3,
                    return.basis = FALSE,
                    return.GCV = FALSE) {
  # check inputs
  if (is.null(X_pred) || is.null(Y) || is.null(id.vec)) { stop("You forgot some inputs to fitGLMM().") }
  if (is.unsorted(id.vec)) { stop("Your data must be ordered by subject, please do so before running fitGLMM().") }

  # fit NB GLMM
  if (adaptive) {
    glm_marge_knots <- purrr::map_dfr(unique(id.vec),
                                      function(x) {
                                        marge_mod_sub <- marge2(X_pred = X_pred[which(id.vec == x), , drop = FALSE],
                                                                Y = Y[which(id.vec == x)],
                                                                approx.knot = approx.knot,
                                                                M = M.glm,
                                                                return.basis = TRUE)
                                        knot_df <- data.frame(knot = extractBreakpoints(marge_mod_sub)$Breakpoint,
                                                              coef = names(coef(marge_mod_sub$final_mod)[-1])) %>%
                                                   mutate(tp_fun = case_when(grepl("-[0-9]", coef) ~ "tp1",
                                                                             TRUE ~ "tp2"))
                                        return(knot_df)
                                      })
    glmm_basis_df <- purrr::pmap_dfc(list(glm_marge_knots$knot,
                                          glm_marge_knots$tp_fun,
                                          1:nrow(glm_marge_knots)),
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
    marge_style_names <- glm_marge_knots$coef
    # build formula automatically
    mod_formula <- as.formula(paste0("Y ~ ",
                                     paste(colnames(glmm_basis_df), collapse = " + "),
                                     " + (1 | subject) + ",
                                     "(1 + ", paste(colnames(glmm_basis_df), collapse = " + "),
                                     " | subject)"))
    glmm_basis_df <- dplyr::mutate(glmm_basis_df,
                                   Y = Y,
                                   subject = id.vec,
                                   .before = 1)
    glmm_mod <- glmmTMB::glmmTMB(mod_formula,
                                 data = glmm_basis_df,
                                 family = glmmTMB::nbinom1(),
                                 se = TRUE,
                                 REML = FALSE)
  } else {
    glmm_basis_df <- data.frame(X1 = tp1(X_pred[, 1], t = round(as.numeric(stats::quantile(X_pred[, 1], 1/3)), 2)),
                                X2 = tp1(X_pred[, 1], t = round(as.numeric(stats::quantile(X_pred[, 1], 2/3)), 2)),
                                X3 = tp2(X_pred[, 1], t = round(as.numeric(stats::quantile(X_pred[, 1], 1/3)), 2)),
                                X4 = tp2(X_pred[, 1], t = round(as.numeric(stats::quantile(X_pred[, 1], 2/3)), 2)),
                                Y = Y,
                                subject = id.vec)
    marge_style_names <- c(paste0("B_final(", round(as.numeric(stats::quantile(X_pred[, 1], 1/3)), 2), "-Lineage_A)"),
                           paste0("B_final(LineageA-", round(as.numeric(stats::quantile(X_pred[, 1], 1/3)), 2), ")"),
                           paste0("B_final(", round(as.numeric(stats::quantile(X_pred[, 1], 2/3)), 2), "-Lineage_A)"),
                           paste0("B_final(LineageA-", round(as.numeric(stats::quantile(X_pred[, 1], 2/3)), 2), ")"))
    glmm_mod <- glmmTMB::glmmTMB(Y ~ X1 + X2 + X3 + X4 + (1 | subject) + (1 + X1 + X2 + X3 + X4 | subject),
                                 data = glmm_basis_df,
                                 family = glmmTMB::nbinom1(),
                                 se = TRUE,
                                 REML = FALSE)
  }
  # set up results
  marge_style_names <- c("B_finalIntercept", marge_style_names)
  res <- list(final_mod = glmm_mod,
              bx = NULL,
              wic_mat = NULL,  # only used in GLM / GEE models
              min_wic_own = NULL,  # only used in GLM / GEE models
              GCV = NULL,
              model_type = "GLMM",
              marge_coef_names = marge_style_names)
  if (return.basis) {
    res$bx <- dplyr::select(glmm_basis_df, -c(subject, Y)) %>%
              as.matrix()
  }
  if (return.GCV) {
    NN <- length(Y)
    p_2 <- attributes(stats::logLik(glmm_mod))$df  # degrees of freedom, differs from how it's done in GLM / GEE marge2() -- might need to change to approximation method later idk
    df1a <- p_2 + 2 * (p_2 - 1) / 2  # This matches the {earth} package, SAS and Friedman (1991) penalty & thus pen = 2 as in marge2()
    RSS1 <- sum((Y - stats::fitted(glmm_mod))^2)
    GCV1 <- RSS1 / (NN * (1 - df1a / NN)^2)  # Note there are some problems w/ GCV in correlated data, so use this as a metric of GoF w/ caution
    res$GCV <- GCV1
  }
  class(res) <- "marge"
  return(res)
}