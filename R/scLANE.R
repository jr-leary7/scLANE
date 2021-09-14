#' Identify differentially expressed genes over pseudotime.
#'
#' @name scLANE
#' @description This function uses GLM spline models as implemented in the \code{marge} package to identify non-linear patterns in gene expression across pseudotime.
#' @import marge
#' @import foreach
#' @import doParallel
#' @param expr A data.frame of integer-valued counts. Defaults to NULL.
#' @param pt A data.frame containing a single column - the pseudotime or latent time estimates for each cell. Defaults to NULL.
#' @param genes A character vector of genes to model. Defaults to NULL.
#' @param n.boot The number of bootstrap samples to use in calculating the \code{marge} model's \emph{p}-value. Defaults to 200.
#' @param track.time A boolean indicating whether the time it takes to execute \code{scLANE()} should be tracked & printed at the end. Defaults to TRUE.
#' @return A list of lists; each sublist contains a gene name, default \code{marge} vs. null model test results & model statistics, bootstrapped \code{marge} vs. null model test results & model statistics, and a \code{ggplot} of the models' fitted values.
#' @export
#' @example
#' scLANE(expr = raw_counts, pt = pt_df, genes = my_gene_vec)

scLANE <- function(expr = NULL,
                   pt = NULL,
                   genes = NULL,
                   n.boot = 200,
                   track.time = TRUE) {
  # check inputs
  if (any(sapply(c(expr, pt, genes), is.null))) stop("You forgot one of the three main inputs.")
  # set up parallelization
  start_time <- Sys.time()
  cl <- makeCluster(detectCores() - 1)
  registerDoParallel(cl)
  set.seed(312)
  test_stats <- foreach(i = seq_along(genes),
                        .combine = "list",
                        .multicombine = TRUE,
                        .maxcombine = length(genes),
                        .packages = c("glm2", "marge", "ggplot2"),
                        .export = c("expr", "genes", "pt", "n.boot", "fit_marge", "plot_marge")) %dopar% {
                          # run original MARGE model
                          gene_data <- expr[, genes[i]]
                          n <- 1
                          N <- nrow(pt)
                          id <- rep(1:N, each = n)
                          marge_mod <- tryCatch (
                            fit_marge(X = pt,
                                      Y = gene_data,
                                      N = N,
                                      n = n,
                                      id = id,
                                      family = "poisson",
                                      tol = 0.0001,
                                      M = 5,
                                      nb = TRUE),
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
                                             Model_Status = "Original model error",
                                             Plot = NA)
                          } else {
                            null_mod <- MASS::glm.nb(gene_data ~ 1)
                            model_plot <- PlotMARGE(model = marge_mod$final_mod[[1]],
                                                    gene.counts = gene_data,
                                                    pt = pt,
                                                    gene = genes[i])
                            marge_ll <- logLik(marge_mod$final_mod[[1]])
                            null_ll <- logLik(null_mod)
                            lrt_stat <- as.numeric(2 * (marge_ll - null_ll))
                            dgr_free <- length(coef(marge_mod$final_mod[[1]])) - 1
                            if (dgr_free == 0) dgr_free <- 1
                            p_val <- pchisq(lrt_stat, dgr_free, lower = FALSE)
                            # estimate p-value w/ parametric bootstrap
                            boot_stats <- rep(NA, n.boot)
                            set.seed(629)
                            for(j in seq(n.boot)){
                              sim_counts <- unlist(simulate(null_mod))
                              sim_null <- MASS::glm.nb(sim_counts ~ 1)
                              sim_marge <- tryCatch(
                                fit_marge(X = pt,
                                          Y = sim_counts,
                                          N = N,
                                          n = n,
                                          id = id,
                                          family = "poisson",
                                          tol = 0.0001,
                                          M = 5,
                                          nb = TRUE),
                                error = function(e) "Model error"
                              )
                              if (any(sim_marge == "Model error")) {
                                boot_stats[j] <- NA
                                next
                              }
                              # simulated LRT statistic
                              sim_marge_ll <- logLik(sim_marge$final_mod[[1]])
                              sim_null_ll <- logLik(sim_null)
                              boot_stats[j] <- as.numeric(2 * (sim_marge_ll - sim_null_ll))
                            }
                            # prepare results
                            boot_stats <- na.omit(boot_stats)
                            if (length(boot_stats) < .9 * n.boot) {
                              res_list <- list(Gene = genes[i],
                                               LRT_Stat = lrt_stat,
                                               P_Val = p_val,
                                               Boot_P_Val = "Bootstrap model error",
                                               Boot_P_Val_SE = "Bootstrap model error",
                                               LogLik_MARGE = as.numeric(logLik(marge_mod$final_mod[[1]])),
                                               LogLik_Null = as.numeric(logLik(null_mod)),
                                               Dev_MARGE = deviance(marge_mod$final_mod[[1]]),
                                               Dev_Null = deviance(null_mod),
                                               Model_Status = "Original model OK, errors in boostrap",
                                               Plot = model_plot)
                            } else {
                              boot_p_val <- mean(boot_stats > lrt_stat)
                              boot_se <- sqrt(boot_p_val * (1 - boot_p_val) / n.boot)
                              res_list <- list(Gene = genes[i],
                                               LRT_Stat = lrt_stat,
                                               P_Val = p_val,
                                               Boot_P_Val = boot_p_val,
                                               Boot_P_Val_SE = boot_se,
                                               LogLik_MARGE = as.numeric(logLik(marge_mod$final_mod[[1]])),
                                               LogLik_Null = as.numeric(logLik(null_mod)),
                                               Dev_MARGE = deviance(marge_mod$final_mod[[1]]),
                                               Dev_Null = deviance(null_mod),
                                               Model_Status = "All clear",
                                               Plot = model_plot)
                            }
                          }
                          res_list
                        }
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
