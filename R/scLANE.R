#' Identify differentially expressed genes over pseudotime.
#'
#' @name scLANE
#' @description This function uses adaptive NB GLM spline models as implemented in the \code{marge} package to identify non-linear patterns in gene expression across pseudotime.
#' @importFrom foreach foreach
#' @importFrom doParallel registerDoParallel
#' @importFrom parallel makeCluster detectCores stopCluster clusterEvalQ
#' @param expr A matrix of integer-valued counts. Defaults to NULL.
#' @param pt A data.frame containing a single column - the pseudotime or latent time estimates for each cell. Defaults to NULL.
#' @param genes A character vector of genes to model. Defaults to NULL.
#' @param n.cores The number of CPUs to use when fitting models in parallel. Defaults to \code{parallel::detectCores() - 1}.
#' @param n.boot The number of bootstrap samples to use in calculating the \code{marge} model's \emph{p}-value. Defaults to 200.
#' @param track.time A boolean indicating whether the time it takes to execute \code{scLANE()} should be tracked & printed at the end. Defaults to TRUE.
#' @return A list of lists; each sublist contains a gene name, default \code{marge} vs. null model test results & model statistics, bootstrapped \code{marge} vs. null model test results & model statistics, and a \code{ggplot} of the models' fitted values.
#' @export
#' @examples
#' scLANE(expr = raw_counts, pt = psuedotime_df, genes = my_gene_vec)

scLANE <- function(expr = NULL,
                   pt = NULL,
                   genes = NULL,
                   n.cores = NULL,
                   n.boot = 200,
                   track.time = TRUE) {
  # check inputs
  if (any(sapply(c(expr, pt, genes), is.null))) stop("You forgot one of the three main inputs.")
  if (any(class(expr) != c("matrix", "array"))) stop("`expr` must be a matrix of integer counts.")
  if (class(pt) != "data.frame") stop("`pt` must be a data.frame.")
  # set up parallelization
  if (track.time) { start_time <- Sys.time() }
  cl <- parallel::makeCluster(ifelse(is.null(n.cores), parallel::detectCores() - 1, n.cores))
  doParallel::registerDoParallel(cl)
  sink(tempfile())
  parallel::clusterEvalQ(cl, {
    set.seed(312)
    NULL
  })
  sink()
  test_stats <- foreach::foreach(i = seq_along(genes),
                                 .combine = "list",
                                 .multicombine = TRUE,
                                 .maxcombine = length(genes),
                                 .packages = c("glm2", "scLANE", "ggplot2"),
                                 .export = c("expr", "genes", "pt", "n.boot")) %dopar% {
    # run original MARGE model
    gene_data <- expr[, genes[i]]
    marge_mod <- tryCatch (
      scLANE::marge2(X_pred = pt,
                     Y = gene_data,
                     M = 7),
      error = function(e) "Model error"
    )
    if (any(marge_mod == "Model error")) {
      res_list <- list(Gene = genes[i],
                       LRT_Stat = NA,
                       P_Val = NA,
                       Boot_P_Val = NA,
                       Boot_P_Val_SE = NA,
                       LogLik_MARGE = NA,
                       LogLik_Null = NA,
                       Dev_MARGE = NA,
                       Dev_Null = NA,
                       Model_Status = "Original MARGE model error",
                       MARGE_model = NA,
                       Plot = NA)
    } else {
      # compare MARGE to null model (must use NB GLM b/c later we simulate from it's distribution, lm doesn't work correctly)
      null_mod <- MASS::glm.nb(gene_data ~ 1)
      model_plot <- scLANE::PlotMARGE(model = marge_mod$final_mod,
                                      gene.counts = gene_data,
                                      pt = pt,
                                      gene = genes[i],
                                      null.mod = null_mod)
      # compute LRT stat using asymptotic \Chi^2 approximation (likely incorrect)
      marge_ll <- logLik(marge_mod$final_mod)
      null_ll <- logLik(null_mod)
      lrt_stat <- as.numeric(-2 * (marge_ll - null_ll))
      dgr_free <- length(coef(marge_mod$final_mod)) - 1
      if (dgr_free == 0) { dgr_free <- 1 }
      p_val <- pchisq(lrt_stat, dgr_free, lower = FALSE)
      # estimate p-value w/ parametric bootstrap
      boot_stats <- rep(NA, n.boot)
      for(j in seq(n.boot)){
        sim_counts <- unlist(simulate(null_mod))
        sim_null <- MASS::glm.nb(sim_counts ~ 1)
        sim_marge <- tryCatch(
          scLANE::marge2(X_pred = pt,
                         Y = sim_counts,
                         M = 7),
          error = function(e) "Model error"
        )
        if (any(sim_marge == "Model error")) {
          boot_stats[j] <- NA
          next
        }
        # simulated LRT statistic
        sim_marge_ll <- logLik(sim_marge$final_mod)
        sim_null_ll <- logLik(sim_null)
        boot_stats[j] <- as.numeric(2 * (sim_marge_ll - sim_null_ll))
      }
      # prepare results
      boot_stats <- na.omit(boot_stats)
      if (length(boot_stats) < .9 * n.boot) {
        res_list <- list(Gene = genes[i],
                         LRT_Stat = lrt_stat,
                         P_Val = p_val,
                         Boot_P_Val = "Bootstrap MARGE model error",
                         Boot_P_Val_SE = "Bootstrap MARGE model error",
                         LogLik_MARGE = as.numeric(logLik(marge_mod$final_mod)),
                         LogLik_Null = as.numeric(logLik(null_mod)),
                         Dev_MARGE = deviance(marge_mod$final_mod),
                         Dev_Null = deviance(null_mod),
                         Model_Status = "Original MARGE model OK, errors during bootstrap",
                         MARGE_model = marge_mod,
                         Plot = model_plot)
      } else {
        boot_p_val <- mean(boot_stats > lrt_stat)
        boot_se <- sqrt(boot_p_val * (1 - boot_p_val) / n.boot)
        res_list <- list(Gene = genes[i],
                         LRT_Stat = lrt_stat,
                         P_Val = p_val,
                         Boot_P_Val = boot_p_val,
                         Boot_P_Val_SE = boot_se,
                         LogLik_MARGE = as.numeric(logLik(marge_mod$final_mod)),
                         LogLik_Null = as.numeric(logLik(null_mod)),
                         Dev_MARGE = deviance(marge_mod$final_mod),
                         Dev_Null = deviance(null_mod),
                         Model_Status = "All clear",
                         MARGE_model = marge_mod,
                         Plot = model_plot)
      }
    }
    res_list
  }
  # end parallelization
  stopCluster(cl)
  names(test_stats) <- genes
  end_time <- Sys.time()
  if (track.time) {
    total_time <- end_time - start_time
    total_time_units <- attributes(total_time)$units
    total_time_numeric <- as.numeric(total_time)
    print(sprintf("scLANE evaluated %s genes in %s %s",
                  length(genes),
                  round(total_time_numeric, 3),
                  total_time_units))
  }
  return(test_stats)
}
