#' Test whether a gene is dynamic over pseudotime.
#'
#' @name testDynamic
#' @description This function tests whether a NB \code{marge} model is better than a null (intercept-only) NB GLM using the Likelihood Ratio Test. In effect, the test tells us whether a gene's expression changes (in any way) over pseudotime.
#' @importFrom foreach foreach
#' @importFrom doParallel registerDoParallel
#' @importFrom parallel makeCluster detectCores stopCluster clusterEvalQ
#' @importFrom MASS glm.nb
#' @param expr A matrix of integer-valued counts. Defaults to NULL.
#' @param pt A data.frame containing a single column - the pseudotime or latent time estimates for each cell. Defaults to NULL.
#' @param genes A character vector of genes to model. If not provided, defaults to all genes in \code{expr}. Defaults to NULL.
#' @param parallel A boolean indicating whether a parallel \code{foreach} loop should be used to generate results more quickly. Defaults to FALSE.
#' @param n.cores (Optional) If running in parallel, how many cores should be used? Defaults to 2.
#' @param track.time A boolean indicating whether the amount of time the function takes to run should be tracked and printed to the console. Useful for debugging. Defaults to FALSE.
#' @return A list of lists; each sublist contains a gene name, default \code{marge} vs. null model test results & model statistics, and a \code{ggplot} of the models' fitted values. Use \code{\link{getResultsDE}} to tidy the results.
#' @export
#' @examples
#' testDynamic(expr = raw_counts, pt = pseudotime_df)
#' testDynamic(expr = raw_counts, pt = pseudotime_df, parallel = TRUE, n.cores = 8)

testDynamic <- function(expr = NULL,
                        pt = NULL,
                        genes = NULL,
                        parallel = FALSE,
                        n.cores = 2,
                        track.time = FALSE) {
  # check inputs
  if (is.null(expr) || is.null(pt)) { stop("You forgot some inputs to testDynamic().") }
  if (any(class(expr) != c("matrix", "array"))) stop("Input expr must be a matrix of integer counts.")
  if (class(pt) != "data.frame") stop("Input pt must be a data.frame.")
  if (is.null(genes)) { genes <- colnames(expr) }
  if (track.time) { start_time <- Sys.time() }
  # set up parallel operation if desired
  `%op%` <- ifelse(parallel, foreach::`%dopar%`, foreach::`%do%`)
  if (parallel) {
    cl <- parallel::makeCluster(ifelse(is.null(n.cores), parallel::detectCores() - 1, n.cores))
    doParallel::registerDoParallel(cl)
    sink(tempfile())
    parallel::clusterEvalQ(cl, {
      set.seed(312)
      NULL
    })
    sink()
  }
  # iterate over genes & build models
  test_stats <- foreach::foreach(i = seq_along(genes),
                                 .combine = "list",
                                 .multicombine = ifelse(length(genes) > 1, TRUE, FALSE),
                                 .maxcombine = ifelse(length(genes) > 1, length(genes), 2),
                                 .packages = c("glm2", "scLANE", "ggplot2", "MASS"),
                                 .export = c("expr", "genes", "pt")) %op% {
    # run original MARGE model
    gene_data <- expr[, genes[i]]
    marge_mod <- tryCatch (
      scLANE::marge2(X_pred = pt,
                     Y = gene_data,
                     M = 7),
      error = function(e) "Model error"
    )

    # fit null model for comparison - must use MASS::glm.nb() because log-likelihood differs when using lm()
    null_mod <- MASS::glm.nb(gene_data ~ 1, method = "glm.fit2", y = FALSE, model = FALSE)
    null_mod <- stripGLM(glm.obj = null_mod)
    if (any(marge_mod == "Model error")) {
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
    } else {
      # compute LRT stat using asymptotic Chi-squared approximation
      lrt_res <- scLANE::modelLRT(mod.1 = marge_mod$final_mod, mod.0 = null_mod)
      marge_ll <- lrt_res$Alt_Mod_LL
      null_ll <- lrt_res$Null_Mod_LL
      lrt_stat <- lrt_res$LRT_Stat
      p_val <- lrt_res$P_Val
      # prepare results
      marge_dev <- deviance(marge_mod$final_mod)
      null_dev <- deviance(null_mod)
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
  if (parallel) { stopCluster(cl) }
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
