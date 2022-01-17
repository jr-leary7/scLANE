#' Test whether a gene is dynamic over pseudotime.
#'
#' @name testDynamic
#' @author Jack Leary
#' @description This function tests whether a NB \code{marge} model is better than a null (intercept-only) NB GLM using the Likelihood Ratio Test. In effect, the test tells us whether a gene's expression changes (in any way) over pseudotime.
#' @importFrom bigstatsr as_FBM
#' @importFrom foreach foreach %dopar% registerDoSEQ
#' @importFrom doParallel registerDoParallel
#' @importFrom parallel makeCluster detectCores stopCluster clusterEvalQ clusterExport
#' @importFrom MASS glm.nb
#' @importFrom glm2 glm.fit2
#' @importFrom dplyr rename mutate relocate
#' @importFrom broom tidy
#' @importFrom stats deviance
#' @param expr.mat A dense matrix of integer-valued counts. Defaults to NULL.
#' @param pt A data.frame containing a single column - the pseudotime or latent time estimates for each cell. Defaults to NULL.
#' @param genes A character vector of genes to model. If not provided, defaults to all genes in \code{expr.mat}. Defaults to NULL.
#' @param parallel.exec A boolean indicating whether a parallel \code{foreach} loop should be used to generate results more quickly. Defaults to FALSE.
#' @param n.cores (Optional) If running in parallel, how many cores should be used? Defaults to 2.
#' @param n.potential.basis.fns The number of possible basis functions. See the parameter "M" in \code{\link{marge2}}. Defaults to 5.
#' @param track.time A boolean indicating whether the amount of time the function takes to run should be tracked and printed to the console. Useful for debugging. Defaults to FALSE.
#' @return A list of lists; each sublist contains a gene name, default \code{marge} vs. null model test results, model statistics, and fitted values. Use \code{\link{getResultsDE}} to tidy the results.
#' @seealso \code{\link{getResultsDE}}
#' @export
#' @examples
#' \dontrun{testDynamic(expr.mat = raw_counts, pt = pseudotime_df, parallel.exec = false)}
#' \dontrun{testDynamic(expr.mat = raw_counts,
#'                      pt = pseudotime_df,
#'                      parallel.exec = TRUE,
#'                      n.cores = 8,
#'                      n.potential.basis.fns = 7)}

testDynamic <- function(expr.mat = NULL,
                        pt = NULL,
                        genes = NULL,
                        parallel.exec = FALSE,
                        n.cores = 2,
                        n.potential.basis.fns = 5,
                        track.time = FALSE) {
  # check inputs
  if (is.null(expr.mat) || is.null(pt)) { stop("You forgot some inputs to testDynamic().") }
  if (any(class(expr.mat) != c("matrix", "array"))) stop("Input expr.mat must be a matrix of integer counts.")
  if (class(pt) != "data.frame") stop("Input pt must be of class data.frame.")
  if (track.time) { start_time <- Sys.time() }
  if (is.null(genes)) { genes <- colnames(expr.mat) }
  if (parallel.exec) {
    cl <- parallel::makeCluster(n.cores)
    doParallel::registerDoParallel(cl)
    clusterEvalQ(cl, expr = {
      set.seed(312)
      rm(list = ls())
    })
  } else {
    cl <- foreach::registerDoSEQ()
    set.seed(312)
  }
  expr.mat <- bigstatsr::as_FBM(expr.mat)
  # build list of objects to prevent from being sent to workers
  if (any(ls(envir = .GlobalEnv) %in% c("expr.mat", "genes", "pt", "n.potential.basis.fns"))) {
    no_export <-  c(ls(envir = .GlobalEnv)[-which(ls(envir = .GlobalEnv) %in% c("expr.mat", "genes", "pt", "n.potential.basis.fns"))],
                    ls()[-which(ls() %in% c("expr.mat", "genes", "pt", "n.potential.basis.fns"))])
  } else {
    no_export <- c(ls(envir = .GlobalEnv),
                   ls()[-which(ls() %in% c("expr.mat", "genes", "pt", "n.potential.basis.fns"))])
  }
  no_export <- unique(no_export)
  # iterate over genes & build models
  test_stats <- foreach::foreach(i = seq_along(genes),
                                 .combine = "list",
                                 .multicombine = ifelse(length(genes) > 1, TRUE, FALSE),
                                 .maxcombine = ifelse(length(genes) > 1, length(genes), 2),
                                 .packages = c("glm2", "scLANE", "MASS",  "bigstatsr", "broom", "dplyr", "stats"),
                                 .noexport = no_export,
                                 .verbose = FALSE) %dopar% {
    # run original MARGE model
    marge_mod <- tryCatch(
      scLANE::marge2(X_pred = pt,
                     Y = expr.mat[, i],
                     M = n.potential.basis.fns),
      error = function(e) "Model error"
    )
    # fit null model for comparison - must use MASS::glm.nb() because log-likelihood differs when using lm()
    null_mod <- tryCatch(
      MASS::glm.nb(expr.mat[, i] ~ 1,
                   method = "glm.fit2",
                   y = FALSE,
                   model = FALSE),
      error = function(e) "Model error"
    )
    if (all(class(null_mod) != "character")) {
      null_mod <- scLANE::stripGLM(glm.obj = null_mod)
    }
    # prepare results if there were errors in either null or MARGE model
    if (unname(marge_mod == "Model error" & all(class(null_mod) == "character"))) {
      res_list <- list(Gene = genes[i],
                       LRT_Stat = NA,
                       P_Val = NA,
                       LogLik_MARGE = NA,
                       LogLik_Null = NA,
                       Dev_MARGE = NA,
                       Dev_Null = NA,
                       Model_Status = "Null and MARGE model errors",
                       MARGE_Summary = NA,
                       Null_Summary = NA,
                       MARGE_Preds = NA,
                       Null_Preds = NA,
                       MARGE_Slope_Data = NA)
    } else if (unname(marge_mod == "Model error" & !all(class(null_mod) == "character"))) {
      # generate null model summary table
      null_sumy_df <- broom::tidy(null_mod) %>% as.data.frame()  # saves a few bytes
      # compute fitted values + SE for null model
      null_pred_df <- data.frame(stats::predict(null_mod, type = "link", se.fit = TRUE)[1:2]) %>%
                      dplyr::rename(null_link_fit = fit, null_link_se = se.fit)
      res_list <- list(Gene = genes[i],
                       LRT_Stat = NA,
                       P_Val = NA,
                       LogLik_MARGE = NA,
                       LogLik_Null = as.numeric(stats::logLik(null_mod)),
                       Dev_MARGE = NA,
                       Dev_Null = stats::deviance(null_mod),
                       Model_Status = "MARGE model error, null model OK",
                       MARGE_Summary = NA,
                       Null_Summary = null_sumy_df,
                       MARGE_Preds = NA,
                       Null_Preds = null_pred_df,
                       MARGE_Slope_Data = NA)
    } else if (unname(marge_mod != "Model error" & !all(class(null_mod) == "character"))) {
      # compute fitted values + SE for both models
      null_pred_df <- tryCatch(
        data.frame(stats::predict(null_mod, type = "link", se.fit = TRUE)[1:2]) %>%
        dplyr::rename(null_link_fit = fit, null_link_se = se.fit),
        error = function(e) "Null model prediction error"
      )
      marge_pred_df <- tryCatch(
        data.frame(stats::predict(marge_mod$final_mod, type = "link", se.fit = TRUE)[1:2]) %>%
        dplyr::rename(marge_link_fit = fit, marge_link_se = se.fit),
        error = function(e) "MARGE model prediction error"
      )
      # generate data for slope test
      marge_slope_df <- scLANE:::createSlopeTestData(marge.model = marge_mod, pt = pt) %>%
                        dplyr::mutate(Gene = genes[i]) %>%
                        dplyr::relocate(Gene)
      # generate model summary tables
      null_sumy_df <- broom::tidy(null_mod) %>% as.data.frame()
      marge_sumy_df <- broom::tidy(marge_mod$final_mod) %>%
                       as.data.frame() %>%
                       lapply(FUN = unname) %>%
                       as.data.frame()
      # compute LRT stat using asymptotic Chi-squared approximation
      lrt_res <- scLANE::modelLRT(mod.1 = marge_mod$final_mod, mod.0 = null_mod)
      # prepare results
      if (class(null_pred_df) == "character" & class(marge_pred_df) == "character") {
        mod_status <- "MARGE & null model prediction errors"
      } else if (class(null_pred_df) == "character" & class(marge_pred_df) == "data.frame") {
        mod_status <- "Null model prediction error"
      } else if (class(null_pred_df) == "data.frame" & class(marge_pred_df) == "character") {
        mod_status <- "MARGE model prediction error"
      } else {
        mod_status <- "All clear"
      }
      res_list <- list(Gene = genes[i],
                       LRT_Stat = lrt_res$LRT_Stat,
                       P_Val = lrt_res$P_Val,
                       LogLik_MARGE = lrt_res$Alt_Mod_LL,
                       LogLik_Null = lrt_res$Null_Mod_LL,
                       Dev_MARGE = stats::deviance(marge_mod$final_mod),
                       Dev_Null = stats::deviance(null_mod),
                       Model_Status = mod_status,
                       MARGE_Summary = marge_sumy_df,
                       Null_Summary = null_sumy_df,
                       MARGE_Preds = marge_pred_df,
                       Null_Preds = null_pred_df,
                       MARGE_Slope_Data = marge_slope_df)
    } else {
      stop(sprintf("Conditions for marge & null model convergence incorrect in testDynamic() - iteration %s.", i))
    }
    rm(marge_mod, null_mod)  # a vain attempt to conserve memory
    res_list
  }
  # end parallelization & clean up
  sink(tempfile())
  if (parallel.exec) {
    parallel::clusterEvalQ(cl, expr = {
      rm(list = ls(all.names = TRUE)); gc(verbose = FALSE, full = TRUE)
    })
    parallel::stopCluster(cl)
  }
  rm(cl)
  gc(verbose = FALSE, full = TRUE)
  sink()
  # prepare results
  names(test_stats) <- genes
  if (track.time) {
    end_time <- Sys.time()
    total_time <- end_time - start_time
    total_time_units <- attributes(total_time)$units
    total_time_numeric <- as.numeric(total_time)
    print(sprintf("testDynamic evaluated %s genes in %s %s",
                  length(genes),
                  round(total_time_numeric, 3),
                  total_time_units))
  }
  return(test_stats)
}
