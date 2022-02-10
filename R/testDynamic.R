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
#' @return A list of list, where each element is a gene and each gene contains sublists for each element. Each gene-lineage sublist contains a gene name, lineage number, default \code{marge} vs. null model test results, model statistics, and fitted values. Use \code{\link{getResultsDE}} to tidy the results.
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
  n_lineages <- ncol(pt)
  # build models per-lineage per-gene, parallelize over genes
  test_stats <- foreach::foreach(i = seq_along(genes),
                                 .combine = "list",
                                 .multicombine = ifelse(length(genes) > 1, TRUE, FALSE),
                                 .maxcombine = ifelse(length(genes) > 1, length(genes), 2),
                                 .packages = c("glm2", "scLANE", "MASS",  "bigstatsr", "broom", "dplyr", "stats"),
                                 .noexport = no_export,
                                 .verbose = FALSE) %dopar% {
    lineage_list <- vector("list", n_lineages)
    for (j in seq(n_lineages)) {
      lineage_cells <- which(!is.na(pt[, j]))
      # run original MARGE model
      marge_mod <- try(
        { scLANE::marge2(X_pred = pt[lineage_cells, j, drop = FALSE],
                         Y = expr.mat[lineage_cells, i, drop = FALSE],
                         M = n.potential.basis.fns) },
        silent = TRUE
      )
      # fit null model for comparison - must use MASS::glm.nb() because log-likelihood differs when using lm()
      null_mod <- try(
        { MASS::glm.nb(expr.mat[lineage_cells, i, drop = FALSE] ~ 1,
                       method = "glm.fit2",
                       y = FALSE,
                       model = FALSE) },
        silent = TRUE
      )
      if (all(class(null_mod) != "try-error")) {
        null_mod <- scLANE::stripGLM(glm.obj = null_mod)
      }
      # prepare results if there were errors in either null or MARGE model
      if (all(class(marge_mod) == "try-error") & all(class(null_mod) == "try-error")) {
        mod_status <- "MARGE & null model errors"
        # generate empty dataframe for slope test
        slope_data_error <- data.frame(Gene = genes[i],
                                       Lineage = j,
                                       Breakpoint = NA_real_,
                                       Rounded_Breakpoint = NA_real_,
                                       Direction = NA_character_,
                                       P_Val = NA_real_,
                                       Notes = mod_status)
        # create lineage result list
        lineage_list[[j]] <- list(Gene = genes[i],
                                  Lineage = j,
                                  LRT_Stat = NA_real_,
                                  P_Val = NA_real_,
                                  LogLik_MARGE = NA_real_,
                                  LogLik_Null = NA_real_,
                                  Dev_MARGE = NA_real_,
                                  Dev_Null = NA_real_,
                                  Model_Status = mod_status,
                                  MARGE_Summary = NA,
                                  Null_Summary = NA,
                                  MARGE_Preds = NA,
                                  Null_Preds = NA,
                                  MARGE_Slope_Data = slope_data_error)
      # prepare results if there were errors in MARGE model but not in null model
      } else if (all(class(marge_mod) == "try-error") & !all(class(null_mod) == "try-error")) {
        mod_status <- "MARGE model error, null model OK"
        # generate empty dataframe for slope test
        slope_data_error <- data.frame(Gene = genes[i],
                                       Lineage = j,
                                       Breakpoint = NA_real_,
                                       Rounded_Breakpoint = NA_real_,
                                       Direction = NA_character_,
                                       P_Val = NA_real_,
                                       Notes = mod_status)
        # generate null model summary table
        null_sumy_df <- broom::tidy(null_mod) %>% as.data.frame()  # saves a few bytes by converting from tibble
        # compute fitted values + SE for null model
        null_pred_df <- data.frame(stats::predict(null_mod, type = "link", se.fit = TRUE)[1:2]) %>%
                        dplyr::rename(null_link_fit = fit, null_link_se = se.fit)
        # create lineage result list
        lineage_list[[j]] <- list(Gene = genes[i],
                                  Lineage = j,
                                  LRT_Stat = NA_real_,
                                  P_Val = NA_real_,
                                  LogLik_MARGE = NA_real_,
                                  LogLik_Null = as.numeric(stats::logLik(null_mod)),
                                  Dev_MARGE = NA_real_,
                                  Dev_Null = stats::deviance(null_mod),
                                  Model_Status = mod_status,
                                  MARGE_Summary = NA,
                                  Null_Summary = null_sumy_df,
                                  MARGE_Preds = NA,
                                  Null_Preds = null_pred_df,
                                  MARGE_Slope_Data = slope_data_error)
      # prepare results if there were errors in null model but not in marge model
      } else if (!all(class(marge_mod) == "try-error") & all(class(null_mod) == "try-error")) {
        mod_status <- "MARGE model OK, null model error"
        # compute fitted values + SE for marge model if possible
        marge_pred_df <- try(
          { data.frame(stats::predict(marge_mod$final_mod, type = "link", se.fit = TRUE)[1:2]) %>%
              dplyr::rename(marge_link_fit = fit, marge_link_se = se.fit) },
          silent = TRUE
        )
        # generate data for slope test
        marge_slope_df <- scLANE:::createSlopeTestData(marge.model = marge_mod, pt = pt[lineage_cells, j, drop = FALSE]) %>%
                          dplyr::mutate(Gene = genes[i], Lineage = j) %>%
                          dplyr::relocate(Gene, Lineage)
        # create marge summary table
        marge_sumy_df <- broom::tidy(marge_mod$final_mod) %>%
                         as.data.frame() %>%
                         lapply(unname) %>%
                         as.data.frame()
        # create lineage result list
        lineage_list[[j]] <- list(Gene = genes[i],
                                  Lineage = j,
                                  LRT_Stat = NA_real_,
                                  P_Val = NA_real_,
                                  LogLik_MARGE = as.numeric(stats::logLik(marge_mod$final_mod)),
                                  LogLik_Null = NA_real_,
                                  Dev_MARGE = stats::deviance(marge_mod$final_mod),
                                  Dev_Null = NA_real_,
                                  Model_Status = mod_status,
                                  MARGE_Summary = marge_sumy_df,
                                  Null_Summary = NA,
                                  MARGE_Preds = marge_pred_df,
                                  Null_Preds = NA,
                                  MARGE_Slope_Data = marge_slope_df)

      # prepare results if neither model had errors
      } else if (!all(class(marge_mod) == "try-error") & !all(class(null_mod) == "try-error")) {
        mod_status <- "MARGE & null model OK"
        # compute fitted values + SE for both models (if possible - this very rarely won't work)
        null_pred_df <- try(
          { data.frame(stats::predict(null_mod, type = "link", se.fit = TRUE)[1:2]) %>%
            dplyr::rename(null_link_fit = fit, null_link_se = se.fit) },
          silent = TRUE
        )
        marge_pred_df <- try(
          { data.frame(stats::predict(marge_mod$final_mod, type = "link", se.fit = TRUE)[1:2]) %>%
            dplyr::rename(marge_link_fit = fit, marge_link_se = se.fit) },
          silent = TRUE
        )
        # generate data for slope test
        marge_slope_df <- scLANE:::createSlopeTestData(marge.model = marge_mod, pt = pt[lineage_cells, j, drop = FALSE]) %>%
                          dplyr::mutate(Gene = genes[i], Lineage = j) %>%
                          dplyr::relocate(Gene, Lineage)
        # generate model summary tables
        null_sumy_df <- broom::tidy(null_mod) %>% as.data.frame()
        marge_sumy_df <- broom::tidy(marge_mod$final_mod) %>%
                         as.data.frame() %>%
                         lapply(unname) %>%
                         as.data.frame()
        # compute LRT stat using asymptotic Chi-squared approximation
        lrt_res <- scLANE::modelLRT(mod.1 = marge_mod$final_mod, mod.0 = null_mod)
        # create lineage result list
        lineage_list[[j]] <- list(Gene = genes[i],
                                  Lineage = j,
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
        stop(sprintf("Conditions for marge or null model fits not met for gene %s on lineage %s.", genes[i], j))
      }
    }
    names(lineage_list) <- paste0("Lineage_", seq(n_lineages))
    rm(marge_mod, null_mod)  # a vain attempt to conserve memory
    lineage_list
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
    print(sprintf("testDynamic evaluated %s genes with %s lineages apiece in %s %s",
                  length(genes),
                  n_lineages,
                  round(total_time_numeric, 3),
                  total_time_units))
  }
  class(test_stats) <- "scLANE"
  return(test_stats)
}
