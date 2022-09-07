#' Test whether a gene is dynamic over pseudotime.
#'
#' @name testDynamic
#' @author Jack Leary
#' @description This function tests whether a NB \code{marge} model is better than a null (intercept-only) NB GLM using the Likelihood Ratio Test. In effect, the test tells us whether a gene's expression changes (in any way) over pseudotime.
#' @importFrom BiocGenerics counts
#' @importFrom slingshot slingPseudotime
#' @importFrom bigstatsr as_FBM
#' @importFrom foreach foreach %dopar% registerDoSEQ
#' @importFrom doParallel registerDoParallel
#' @importFrom parallel makeCluster detectCores stopCluster clusterEvalQ clusterExport
#' @importFrom MASS glm.nb negative.binomial
#' @importFrom glm2 glm.fit2
#' @importFrom dplyr rename mutate relocate
#' @importFrom broom tidy
#' @importFrom stats deviance
#' @importFrom geeM geem
#' @param expr.mat Either a \code{SingleCellExperiment} object from which counts can be extracted, or a dense matrix of integer-valued counts. Defaults to NULL.
#' @param pt Either the output from \code{\link[slingshot]{SlingshotDataSet}} object from which pseudotime can be generated, or a data.frame containing the pseudotime or latent time estimates for each cell (can be multiple columns / lineages). Defaults to NULL.
#' @param genes A character vector of genes to model. If not provided, defaults to all genes in \code{expr.mat}. Defaults to NULL.
#' @param is.gee Should a GEE framework be used instead of the default GLM? Defaults to FALSE.
#' @param id.vec If the GEE framework is being used, a vector of subject IDs to use as input to \code{\link[geeM]{geem}}. Defaults to NULL.
#' @param cor.structure If the GEE framework is used, specifies the desired working correlation structure. Must be one of "ar1", "independence", "unstructured", or "exchangeable". Defaults to "exchangeable".
#' @param parallel.exec A boolean indicating whether a parallel \code{foreach} loop should be used to generate results more quickly. Defaults to FALSE.
#' @param n.cores (Optional) If running in parallel, how many cores should be used? Defaults to 2.
#' @param n.potential.basis.fns (Optional) The number of possible basis functions. See the parameter \code{M} in \code{\link{marge2}}. Defaults to 5.
#' @param track.time (Optional) A boolean indicating whether the amount of time the function takes to run should be tracked and printed to the console. Useful for debugging. Defaults to FALSE.
#' @param log.file (Optional) A boolean indicating whether iteration tracking should be printed to \code{"log.txt"}. Can be useful for debugging. Defaults to FALSE.
#' @param log.iter (Optional) If logging is enabled, how often should iterations be printed to the logfile. Defaults to 1, or every iteration.
#' @return A list of lists, where each element is a gene and each gene contains sublists for each element. Each gene-lineage sublist contains a gene name, lineage number, default \code{marge} vs. null model test results, model statistics, and fitted values. Use \code{\link{getResultsDE}} to tidy the results.
#' @seealso \code{\link{getResultsDE}}
#' @seealso \code{\link{testSlope}}
#' @export
#' @examples
#' \dontrun{testDynamic(expr.mat = raw_counts, pt = pseudotime_df, parallel.exec = false)}
#' \dontrun{testDynamic(expr.mat = sce_obj, pt = slingshot_obj, genes = rownames(sce_obj)[1:100])}
#' \dontrun{testDynamic(expr.mat = raw_counts,
#'                      pt = pseudotime_df,
#'                      parallel.exec = TRUE,
#'                      n.cores = 8,
#'                      n.potential.basis.fns = 7)}
#' \dontrun{testDynamic(expr.mat = raw_counts,
#'                      pt = pseudotime_df,
#'                      is.gee = TRUE,
#'                      id.vec = my_subject_ids,
#'                      cor.structure = "independence",
#'                      parallel.exec = TRUE,
#'                      n.cores = 8,
#'                      n.potential.basis.fns = 7)}

testDynamic <- function(expr.mat = NULL,
                        pt = NULL,
                        genes = NULL,
                        is.gee = FALSE,
                        id.vec = NULL,
                        cor.structure = "exchangeable",
                        parallel.exec = FALSE,
                        n.cores = 2,
                        n.potential.basis.fns = 5,
                        track.time = FALSE,
                        log.file = FALSE,
                        log.iter = 1) {
  # check inputs
  if (is.null(expr.mat) || is.null(pt)) { stop("You forgot some inputs to testDynamic().") }
  if (any(class(expr.mat) == "SingleCellExperiment")) {
    expr.mat <- as.matrix(t(BiocGenerics::counts(expr.mat)))
  }
  if (any(class(expr.mat) != c("matrix", "array"))) { stop("Input expr.mat must be a matrix of integer counts.") }
  if (class(pt) == "SlingshotDataSet") {
    pt <- as.data.frame(slingshot::slingPseudotime(pt))
  }
  if (class(pt) != "data.frame") { stop("pt must be of class data.frame.") }
  if (is.gee && is.null(id.vec)) { stop("You must provide a vector of IDs if you're using the GEE framework.") }
  if (is.gee && is.unsorted(id.vec)) { stop("Your data must be ordered by subject, please do so before running testDynamic() with is.gee = TRUE.") }
  cor.structure <- tolower(cor.structure)
  if (is.gee && !(cor.structure %in% c("ar1", "independence", "exchangeable"))) { stop("GEE models require a specified correlation structure.") }
  if (track.time) { start_time <- Sys.time() }
  if (is.null(genes)) { genes <- colnames(expr.mat) }
  # set column names automatically to prevent user error
  n_lineages <- ncol(pt)
  colnames(pt) <- paste0("Lineage_", LETTERS[1:n_lineages])
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
  if (any(ls(envir = .GlobalEnv) %in% c("expr.mat", "genes", "pt", "n.potential.basis.fns", "n_lineages", "id.vec", "cor.structure", "is.gee", "log.file", "log.iter"))) {
    no_export <-  c(ls(envir = .GlobalEnv)[-which(ls(envir = .GlobalEnv) %in% c("expr.mat", "genes", "pt", "n.potential.basis.fns", "n_lineages", "id.vec", "cor.structure", "is.gee", "log.file", "log.iter"))],
                    ls()[-which(ls() %in% c("expr.mat", "genes", "pt", "n.potential.basis.fns", "n_lineages", "id.vec", "cor.structure", "is.gee", "log.file", "log.iter"))])
  } else {
    no_export <- c(ls(envir = .GlobalEnv),
                   ls()[-which(ls() %in% c("expr.mat", "genes", "pt", "n.potential.basis.fns", "n_lineages", "id.vec", "cor.structure", "is.gee", "log.file", "log.iter"))])
  }
  no_export <- unique(no_export)
  # build models per-lineage per-gene, parallelize over genes
  if (log.file) {
    writeLines(c(""), "log.txt")
  }
  test_stats <- foreach::foreach(i = seq_along(genes),
                                 .combine = "list",
                                 .multicombine = ifelse(length(genes) > 1, TRUE, FALSE),
                                 .maxcombine = ifelse(length(genes) > 1, length(genes), 2),
                                 .packages = c("glm2", "scLANE", "MASS",  "bigstatsr", "broom", "dplyr", "stats", "geeM"),
                                 .noexport = no_export,
                                 .verbose = FALSE) %dopar% {
    if (log.file) {
      print_nums <- seq(0, length(genes), log.iter)[-1]
      if (i %in% print_nums) {
        sink("log.txt", append = TRUE)
        cat(paste("Starting iteration:", i, "\n"))
        sink()
      }
    }
    lineage_list <- vector("list", n_lineages)
    for (j in seq(n_lineages)) {
      lineage_cells <- which(!is.na(pt[, j]))
      # run original MARGE model - GEE or not
      if (is.gee) {
        marge_mod <- try(
          { scLANE::marge2(X_pred = pt[lineage_cells, j, drop = FALSE],
                           Y = expr.mat[lineage_cells, i],
                           is.gee = is.gee,
                           id.vec = id.vec[lineage_cells],
                           cor.structure = cor.structure,
                           M = n.potential.basis.fns) },
          silent = TRUE
        )
      } else {
        marge_mod <- try(
          { scLANE::marge2(X_pred = pt[lineage_cells, j, drop = FALSE],
                           Y = expr.mat[lineage_cells, i, drop = FALSE],
                           M = n.potential.basis.fns) },
          silent = TRUE
        )
      }
      # fit null model for comparison - must use MASS::glm.nb() because log-likelihood differs when using lm() if not using a GEE model
      if (is.gee) {
        null_mod <- try(
          { geeM::geem(expr.mat[lineage_cells, i, drop = FALSE] ~ 1,
                       id = id.vec[lineage_cells],
                       family = MASS::negative.binomial(1),
                       corstr = cor.structure) },
          silent = TRUE
        )
      } else {
        null_mod <- try(
          { MASS::glm.nb(expr.mat[lineage_cells, i, drop = FALSE] ~ 1,
                         method = "glm.fit2",
                         y = FALSE,
                         model = FALSE) },
          silent = TRUE
        )
      }
      # slim down GLM object if not a GEE model (which are much smaller for some reason)
      if (!is.gee & all(class(null_mod) != "try-error")) {
        null_mod <- scLANE::stripGLM(glm.obj = null_mod)
      }
      # prepare results if there were errors in either null or MARGE model
      if (all(class(marge_mod) == "try-error") & all(class(null_mod) == "try-error")) {
        mod_status <- "MARGE model error, null model error"
        # generate empty dataframe for slope test
        slope_data_error <- data.frame(Gene = genes[i],
                                       Lineage = LETTERS[j],
                                       Breakpoint = NA_real_,
                                       Rounded_Breakpoint = NA_real_,
                                       Direction = NA_character_,
                                       P_Val = NA_real_,
                                       Notes = mod_status)
        # create lineage result list
        lineage_list[[j]] <- list(Gene = genes[i],
                                  Lineage = LETTERS[j],
                                  Test_Stat = NA_real_,
                                  Test_Stat_Type = ifelse(is.gee, "Wald", "LRT"),
                                  Test_Stat_Note = NA_character_,
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
        # generate null model summary table & fitted values w/ standard errors (need to do so manually for GEE)
        if (is.gee) {
          null_sumy_df <- try({
            null_gee_summary <- summary(null_mod)
            data.frame(term = null_gee_summary$coefnames,
                       estimate = unname(null_gee_summary$beta),
                       `std.error` = unname(null_gee_summary$se.robust),
                       statistic = unname(null_gee_summary$wald.test),
                       `p.value` = unname(null_gee_summary$p))
          }, silent = TRUE)
          null_pred_df <- try({
            robust_vcov_mat <- as.matrix(null_mod$var)
            data.frame(null_link_fit = predict(null_mod),
                       null_link_se = sqrt(apply((tcrossprod(null_mod$X, robust_vcov_mat)) * null_mod$X, 1, sum)))  # wow I love math
          }, silent = TRUE)
        } else {
          null_sumy_df <- try({
            broom::tidy(null_mod) %>% as.data.frame()  # saves a few bytes by converting from tibble
          }, silent = TRUE)
          null_pred_df <- try({
            data.frame(stats::predict(null_mod, type = "link", se.fit = TRUE)[1:2]) %>%
              dplyr::rename(null_link_fit = fit, null_link_se = se.fit)
          }, silent = TRUE)
        }
        # create lineage result list
        lineage_list[[j]] <- list(Gene = genes[i],
                                  Lineage = LETTERS[j],
                                  Test_Stat = NA_real_,
                                  Test_Stat_Type = ifelse(is.gee, "Wald", "LRT"),
                                  Test_Stat_Note = NA_character_,
                                  P_Val = NA_real_,
                                  LogLik_MARGE = NA_real_,
                                  LogLik_Null = ifelse(is.gee, NA_real_, as.numeric(stats::logLik(null_mod))),
                                  Dev_MARGE = NA_real_,
                                  Dev_Null = ifelse(is.gee, NA_real_, stats::deviance(null_mod)),
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
        if (is.gee) {
          marge_pred_df <- try({
            robust_vcov_mat <- as.matrix(marge_mod$final_mod$var)
            data.frame(marge_link_fit = predict(marge_mod$final_mod),
                       marge_link_se = sqrt(apply((tcrossprod(marge_mod$final_mod$X, robust_vcov_mat)) * marge_mod$final_mod$X, 1, sum)))
              }, silent = TRUE)
        } else {
          marge_pred_df <- try({
            data.frame(stats::predict(marge_mod$final_mod, type = "link", se.fit = TRUE)[1:2]) %>%
              dplyr::rename(marge_link_fit = fit, marge_link_se = se.fit)
            }, silent = TRUE)
        }
        # compute summary stat table
        if (is.gee) {
          marge_sumy_df <- try({
            marge_gee_summary <- summary(marge_mod$final_mod)
            data.frame(term = marge_gee_summary$coefnames,
                       estimate = unname(marge_gee_summary$beta),
                       `std.error` = unname(marge_gee_summary$se.robust),
                       statistic = unname(marge_gee_summary$wald.test),
                       `p.value` = unname(marge_gee_summary$p))
          }, silent = TRUE)
        } else {
          marge_sumy_df <- try({
            broom::tidy(marge_mod$final_mod) %>%
              as.data.frame() %>%
              lapply(unname) %>%
              as.data.frame()
          }, silent = TRUE)
        }
        # create slope test input table
        marge_slope_df <- scLANE:::createSlopeTestData(marge.model = marge_mod,
                                                       pt = pt[lineage_cells, j, drop = FALSE],
                                                       is.gee = is.gee) %>%
                          dplyr::mutate(Gene = genes[i], Lineage = LETTERS[j]) %>%
                          dplyr::relocate(Gene, Lineage)
        # create lineage result list
        lineage_list[[j]] <- list(Gene = genes[i],
                                  Lineage = LETTERS[j],
                                  Test_Stat = NA_real_,
                                  Test_Stat_Type = ifelse(is.gee, "Wald", "LRT"),
                                  Test_Stat_Note = NA_character_,
                                  P_Val = NA_real_,
                                  LogLik_MARGE = ifelse(is.gee, NA_real_, as.numeric(stats::logLik(marge_mod$final_mod))),
                                  LogLik_Null = NA_real_,
                                  Dev_MARGE = ifelse(is.gee, NA_real_, stats::deviance(marge_mod$final_mod)),
                                  Dev_Null = NA_real_,
                                  Model_Status = mod_status,
                                  MARGE_Summary = marge_sumy_df,
                                  Null_Summary = NA,
                                  MARGE_Preds = marge_pred_df,
                                  Null_Preds = NA,
                                  MARGE_Slope_Data = marge_slope_df)

      # prepare results if neither model had errors
      } else if (!all(class(marge_mod) == "try-error") & !all(class(null_mod) == "try-error")) {
        mod_status <- "MARGE model OK, null model OK"
        # generate null model summary table & fitted values w/ standard errors (need to do so manually for GEE)
        if (is.gee) {
          null_sumy_df <- try({
            null_gee_summary <- summary(null_mod)
            data.frame(term = null_gee_summary$coefnames,
                       estimate = unname(null_gee_summary$beta),
                       `std.error` = unname(null_gee_summary$se.robust),
                       statistic = unname(null_gee_summary$wald.test),
                       `p.value` = unname(null_gee_summary$p))
          }, silent = TRUE)
          null_pred_df <- try({
            robust_vcov_mat <- as.matrix(null_mod$var)
            data.frame(null_link_fit = predict(null_mod),
                       # can maybe replace below line with null_link_se = unname(summary(null_mod)$se.robust) for a little more speed w/ big N ?
                       null_link_se = sqrt(apply((tcrossprod(null_mod$X, robust_vcov_mat)) * null_mod$X, 1, sum)))  # wow I love math
          }, silent = TRUE)
        } else {
          null_sumy_df <- try({
            broom::tidy(null_mod) %>% as.data.frame()  # saves a few bytes by converting from tibble
          }, silent = TRUE)
          null_pred_df <- try({
            data.frame(stats::predict(null_mod, type = "link", se.fit = TRUE)[1:2]) %>%
              dplyr::rename(null_link_fit = fit, null_link_se = se.fit)
          }, silent = TRUE)
        }
        # generate marge model summary table & fitted values w/ standard errors (need to do so manually for GEE)
        if (is.gee) {
          robust_vcov_mat <- as.matrix(marge_mod$final_mod$var)
          marge_pred_df <- try({
            data.frame(marge_link_fit = predict(marge_mod$final_mod),
                       marge_link_se = sqrt(apply((tcrossprod(marge_mod$final_mod$X, robust_vcov_mat)) * marge_mod$final_mod$X, 1, sum)))
          }, silent = TRUE)
        } else {
          marge_pred_df <- try({
            data.frame(stats::predict(marge_mod$final_mod, type = "link", se.fit = TRUE)[1:2]) %>%
              dplyr::rename(marge_link_fit = fit, marge_link_se = se.fit)
          }, silent = TRUE)
        }
        # compute summary stat table
        if (is.gee) {
          marge_sumy_df <- try({
            marge_gee_summary <- summary(marge_mod$final_mod)
            data.frame(term = marge_gee_summary$coefnames,
                       estimate = unname(marge_gee_summary$beta),
                       `std.error` = unname(marge_gee_summary$se.robust),
                       statistic = unname(marge_gee_summary$wald.test),
                       `p.value` = unname(marge_gee_summary$p))
          }, silent = TRUE)
        } else {
          marge_sumy_df <- try({
            broom::tidy(marge_mod$final_mod) %>%
              as.data.frame() %>%
              lapply(unname) %>%
              as.data.frame()
          }, silent = TRUE)
        }
        # create slope test input table
        marge_slope_df <- scLANE:::createSlopeTestData(marge.model = marge_mod,
                                                       pt = pt[lineage_cells, j, drop = FALSE],
                                                       is.gee = is.gee) %>%
                          dplyr::mutate(Gene = genes[i], Lineage = LETTERS[j]) %>%
                          dplyr::relocate(Gene, Lineage)
        # compute LRT stat using asymptotic Chi-squared approximation
        if (is.gee) {
          test_res <- scLANE::waldTestGEE(mod.1 = marge_mod$final_mod, mod.0 = null_mod)
        } else {
          test_res <- scLANE::modelLRT(mod.1 = marge_mod$final_mod, mod.0 = null_mod)
        }
        # create lineage result list
        lineage_list[[j]] <- list(Gene = genes[i],
                                  Lineage = LETTERS[j],
                                  Test_Stat = ifelse(is.gee, test_res$Wald_Stat, test_res$LRT_Stat),
                                  Test_Stat_Type = ifelse(is.gee, "Wald", "LRT"),
                                  Test_Stat_Note = test_res$Notes,
                                  P_Val = test_res$P_Val,
                                  LogLik_MARGE = ifelse(is.gee, NA_real_, as.numeric(stats::logLik(marge_mod$final_mod))),
                                  LogLik_Null = ifelse(is.gee, NA_real_, as.numeric(stats::logLik(null_mod))),
                                  Dev_MARGE = ifelse(is.gee, NA_real_, stats::deviance(marge_mod$final_mod)),
                                  Dev_Null = ifelse(is.gee, NA_real_, stats::deviance(null_mod)),
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
    names(lineage_list) <- paste0("Lineage_", LETTERS[1:n_lineages])
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
