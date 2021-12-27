#' Test whether a gene is dynamic over pseudotime.
#'
#' @name testDynamic
#' @author Jack Leary
#' @description This function tests whether a NB \code{marge} model is better than a null (intercept-only) NB GLM using the Likelihood Ratio Test. In effect, the test tells us whether a gene's expression changes (in any way) over pseudotime.
#' @importFrom bigmemory as.big.matrix describe attach.big.matrix
#' @importFrom foreach foreach %dopar%
#' @importFrom doParallel registerDoParallel
#' @importFrom parallel makeCluster detectCores stopCluster clusterEvalQ clusterExport
#' @importFrom MASS glm.nb
#' @importFrom glm2 glm.fit2
#' @importFrom stats deviance
#' @param expr A sparse (this is imporant for memory reasons) matrix of integer-valued counts. Defaults to NULL.
#' @param pt A data.frame containing a single column - the pseudotime or latent time estimates for each cell. Defaults to NULL.
#' @param genes A character vector of genes to model. If not provided, defaults to all genes in \code{expr}. Defaults to NULL.
#' @param parallel A boolean indicating whether a parallel \code{foreach} loop should be used to generate results more quickly. Defaults to FALSE.
#' @param n.cores (Optional) If running in parallel, how many cores should be used? Defaults to 2.
#' @param n.potential.basis.fns The number of possible basis functions. See the parameter "M" in \code{\link{marge2}}. Defaults to 5.
#' @param track.time A boolean indicating whether the amount of time the function takes to run should be tracked and printed to the console. Useful for debugging. Defaults to FALSE.
#' @return A list of lists; each sublist contains a gene name, default \code{marge} vs. null model test results & model statistics, and a \code{ggplot} of the models' fitted values. Use \code{\link{getResultsDE}} to tidy the results.
#' @export
#' @examples
#' \dontrun{testDynamic(expr = raw_counts, pt = pseudotime_df)}
#' \dontrun{testDynamic(expr = raw_counts,
#'                      pt = pseudotime_df,
#'                      parallel = TRUE,
#'                      n.cores = 8,
#'                      n.potential.basis.fns = 7)}

testDynamic <- function(expr = NULL,
                        pt = NULL,
                        genes = NULL,
                        parallel = FALSE,
                        n.cores = 2,
                        n.potential.basis.fns = 5,
                        track.time = FALSE) {
  # check inputs
  if (is.null(expr) || is.null(pt)) { stop("You forgot some inputs to testDynamic().") }
  if (any(class(expr) != c("matrix", "array"))) stop("Input expr must be a matrix of integer counts.")
  if (class(pt) != "data.frame") stop("Input pt must be of class data.frame.")
  if (track.time) { start_time <- Sys.time() }
  if (is.null(genes)) { genes <- colnames(expr) }
  # set up parallel operation if desired, otherwise just set n.cores = 1 (this causes a very slight setup time slowdown, but it's negligible compared to total runtime)
  cl <- parallel::makeCluster(ifelse(parallel, n.cores, 1))
  doParallel::registerDoParallel(cl)
  sink(tempfile())
  clusterEvalQ(cl, expr = {
    set.seed(312)
    rm(list = ls()); gc()
    NULL
  })
  expr <- bigmemory::as.big.matrix(expr, shared = TRUE)
  big_mat_desc <- bigmemory::describe(expr)
  clusterEvalQ(cl, expr = {
    library(bigmemory)
    NULL
  })
  clusterExport(cl, "big_mat_desc")
  clusterEvalQ(cl, expr = {
    big_mat2 <- bigmemory::attach.big.matrix(big_mat_desc)
    NULL
  })
  sink()
  # build list of objects to prevent from being sent to workers
  if (any(ls(envir = .GlobalEnv) %in% c("big_mat_desc", "genes", "pt", "n.potential.basis.fns"))) {
    no_export <-  c(ls(envir = .GlobalEnv)[-which(ls(envir = .GlobalEnv) %in% c("big_mat_desc", "genes", "pt", "n.potential.basis.fns"))],
                    ls()[-which(ls() %in% c("big_mat_desc", "genes", "pt", "n.potential.basis.fns"))])
  } else {
    no_export <- c(ls(envir = .GlobalEnv),
                   ls()[-which(ls() %in% c("big_mat_desc", "genes", "pt", "n.potential.basis.fns"))])
  }
  no_export <- unique(no_export)
  # iterate over genes & build models
  test_stats <- foreach::foreach(i = seq_along(genes),
                                 .combine = "list",
                                 .multicombine = ifelse(length(genes) > 1, TRUE, FALSE),
                                 .maxcombine = ifelse(length(genes) > 1, length(genes), 2),
                                 .packages = c("glm2", "scLANE", "MASS",  "bigmemory", "stats"),
                                 .export = c("genes", "pt", "big_mat_desc", "n.potential.basis.fns"),
                                 .noexport = no_export,
                                 .verbose = FALSE) %dopar% {
    # run original MARGE model
    marge_mod <- tryCatch(
      scLANE::marge2(X_pred = pt,
                     Y = big_mat2[, i],
                     M = n.potential.basis.fns),
      error = function(e) "Model error"
    )

    # fit null model for comparison - must use MASS::glm.nb() because log-likelihood differs when using lm()
    null_mod <- tryCatch(
      MASS::glm.nb(big_mat2[, i] ~ 1,
                   method = "glm.fit2",
                   y = FALSE,
                   model = FALSE),
      error = function(e) "Model error"
    )
    if (!any(null_mod == "Model error")) {
      null_mod <- scLANE::stripGLM(glm.obj = null_mod)
    }

    # prepare results if there were errors in either null or MARGE model
    if (any(marge_mod == "Model error")) {
      if (any(null_mod == "Model error")) {
        res_list <- list(Gene = genes[i],
                         LRT_Stat = NA,
                         P_Val = NA,
                         LogLik_MARGE = NA,
                         LogLik_Null = NA,
                         Dev_MARGE = NA,
                         Dev_Null = NA,
                         Model_Status = "Null and MARGE model errors",
                         MARGE_Model = NA,
                         Null_Model = null_mod)
      } else {
        res_list <- list(Gene = genes[i],
                         LRT_Stat = NA,
                         P_Val = NA,
                         LogLik_MARGE = NA,
                         LogLik_Null = as.numeric(stats::logLik(null_mod)),
                         Dev_MARGE = NA,
                         Dev_Null = deviance(null_mod),
                         Model_Status = "Original MARGE model error",
                         MARGE_Model = NA,
                         Null_Model = null_mod)
      }
    } else {
      # compute LRT stat using asymptotic Chi-squared approximation
      lrt_res <- scLANE::modelLRT(mod.1 = marge_mod$final_mod, mod.0 = null_mod)
      marge_ll <- lrt_res$Alt_Mod_LL
      null_ll <- lrt_res$Null_Mod_LL
      lrt_stat <- lrt_res$LRT_Stat
      p_val <- lrt_res$P_Val
      # prepare results
      marge_dev <- stats::deviance(marge_mod$final_mod)
      null_dev <- stats::deviance(null_mod)
      res_list <- list(Gene = genes[i],
                       LRT_Stat = lrt_stat,
                       P_Val = p_val,
                       LogLik_MARGE = marge_ll,
                       LogLik_Null = null_ll,
                       Dev_MARGE = marge_dev,
                       Dev_Null = null_dev,
                       Model_Status = "All clear",
                       MARGE_Model = marge_mod,
                       Null_Model = null_mod)
    }
    res_list
  }
  # end parallelization
  sink(tempfile())
  parallel::clusterEvalQ(cl, expr = {
    rm(list = ls(all.names = TRUE)); gc(verbose = FALSE, full = TRUE)
  })
  parallel::stopCluster(cl)
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
