pull.null.sumy <- function(mod.obj, is.gee, is.glmm) {
  # check inputs
  if (inherits(mod.obj, "try-error")) {
    res <- list(null_sumy_df = NA,
                null_pred_df = NA,
                null_ll = NA_real_,
                null_dev = NA_real_)
  } else {
    # pull model summary
    if (is.gee) {
      null_sumy_df <- try({
        null_gee_summary <- summary(mod.obj)
        data.frame(term = null_gee_summary$coefnames,
                   estimate = unname(null_gee_summary$beta),
                   `std.error` = unname(null_gee_summary$se.robust),
                   statistic = unname(null_gee_summary$wald.test),
                   `p.value` = unname(null_gee_summary$p))
      }, silent = TRUE)
      null_pred_df <- try({
        robust_vcov_mat <- as.matrix(mod.obj$var)
        data.frame(null_link_fit = predict(mod.obj),
                   null_link_se = sqrt(apply((tcrossprod(mod.obj$X, robust_vcov_mat)) * mod.obj$X, 1, sum)))
      }, silent = TRUE)
    } else if (is.glmm) {
      null_sumy_df <- try({
        null_glmm_summary <- as.data.frame(broom.mixed::tidy(mod.obj, effects = "fixed"))
        data.frame(term = null_glmm_summary$term,
                   estimate = null_glmm_summary$estimate[1],
                   `std.error` = null_glmm_summary$std.error[1],
                   statistic = null_glmm_summary$statistic[1],
                   `p.value` = null_glmm_summary$p.value[1])
      }, silent = TRUE)
      null_pred_df <- try({
        data.frame(predict(mod.obj, type = "link", se.fit = TRUE)[1:2]) %>%
          dplyr::rename(null_link_fit = fit, null_link_se = se.fit)
      }, silent = TRUE)
    }  else {
      null_sumy_df <- try({
        as.data.frame(broom.mixed::tidy(mod.obj)) # saves a few bytes by converting from tibble
      }, silent = TRUE)
      null_pred_df <- try({
        data.frame(stats::predict(mod.obj, type = "link", se.fit = TRUE)[1:2]) %>%
          dplyr::rename(null_link_fit = fit, null_link_se = se.fit)
      }, silent = TRUE)
    }

    null_ll <- ifelse(is.gee, NA_real_, as.numeric(stats::logLik(mod.obj)))
    null_dev <- ifelse((is.gee || is.glmm), NA_real_, stats::deviance(mod.obj))

    res <- list(null_sumy_df = null_sumy_df,
                null_pred_df = null_pred_df,
                null_ll = null_ll,
                null_dev = null_dev)
  }
 return(res)
}





pull.marge.sumy <- function(mod.obj,
                            is.gee,
                            is.glmm) {
  # check inputs
  if (inherits(mod.obj, "try-error")) {
    res <- list(marge_pred_df = NA,
                marge_sumy_df = NA,
                ll_marge = NA_real_,
                marge_fit_notes = NA_character_)
  } else {
    # pull model summary
    if (is.gee) {
      robust_vcov_mat <- as.matrix(mod.obj$final_mod$var)
      marge_pred_df <- try({
        data.frame(marge_link_fit = predict(mod.obj$final_mod),
                   marge_link_se = sqrt(apply((tcrossprod(mod.obj$final_mod$X, robust_vcov_mat)) * mod.obj$final_mod$X, 1, sum)))
      }, silent = TRUE)
      marge_sumy_df <- try({
        marge_gee_summary <- summary(mod.obj$final_mod)
        data.frame(term = marge_gee_summary$coefnames,
                   estimate = unname(marge_gee_summary$beta),
                   `std.error` = unname(marge_gee_summary$se.robust),
                   statistic = unname(marge_gee_summary$wald.test),
                   `p.value` = unname(marge_gee_summary$p))
      }, silent = TRUE)

    } else if (is.glmm) {
      marge_pred_df <- try({
        data.frame(predict(mod.obj$final_mod, type = "link", se.fit = TRUE)[1:2]) %>%
          dplyr::rename(marge_link_fit = fit, marge_link_se = se.fit)
      }, silent = TRUE)
      marge_sumy_df <- try({
        marge_glmm_summary <- broom.mixed::tidy(mod.obj$final_mod, effects = "fixed")
        data.frame(term = mod.obj$marge_coef_names,  # use marge-style hinge function names instead of X1, X2, ..., XP (makes interpretability easier)
                   estimate = marge_glmm_summary$estimate,
                   `std.error` = marge_glmm_summary$std.error,
                   statistic = marge_glmm_summary$statistic,
                   `p.value` = marge_glmm_summary$p.value)
      }, silent = TRUE)

    } else {
      marge_pred_df <- try({
        data.frame(stats::predict(mod.obj$final_mod, type = "link", se.fit = TRUE)[1:2]) %>%
          dplyr::rename(marge_link_fit = fit, marge_link_se = se.fit)
      }, silent = TRUE)
      marge_sumy_df <- try({
        as.data.frame(broom.mixed::tidy(mod.obj$final_mod)) %>%
          lapply(unname) %>%
          as.data.frame()
      }, silent = TRUE)
    }

    # get log-likelihood for GLMM or GLM cases
    if (is.glmm) {
      ll_marge <- -mod.obj$final_mod$fit$objective
    } else if (!(is.gee || is.glmm)) {
      ll_marge <- as.numeric(stats::logLik(mod.obj$final_mod))
    } else {
      ll_marge <- NA_real_
    }
    # check positive-definiteness of Hessian for GLMM -- might have an effect on LRT stat / accompanying p-value
    if (is.glmm) {
      if (!mod.obj$final_mod$sdr$pdHess) {
        marge_fit_notes <- "Non-positive definite Hessian in GLMM, probably due to shallow log-likelihood. Be careful!"
      } else {
        marge_fit_notes <- NA_character_
      }
    } else {
      marge_fit_notes <- NA_character_
    }
    marge_dev <- ifelse((is.gee || is.glmm), NA_real_, stats::deviance(mod.obj$final_mod))
    res <- list(marge_pred_df = marge_pred_df,
                marge_sumy_df = marge_sumy_df,
                ll_marge = ll_marge,
                marge_dev = marge_dev,
                marge_fit_notes = marge_fit_notes)
  }
 return(res)
}
