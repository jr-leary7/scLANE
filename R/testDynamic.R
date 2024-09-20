#' Test whether a gene is dynamic over pseudotime.
#'
#' @name testDynamic
#' @author Jack Leary
#' @description This function tests whether a NB \code{marge} model is better than a null (intercept-only) model using the Likelihood Ratio Test. In effect, the test tells us whether a gene's expression changes (in any way) over pseudotime.
#' @import glm2
#' @import magrittr
#' @importFrom Matrix t
#' @importFrom bigstatsr as_FBM
#' @importFrom utils txtProgressBar setTxtProgressBar
#' @importFrom foreach foreach %dopar% registerDoSEQ
#' @importFrom doSNOW registerDoSNOW
#' @importFrom parallel makeCluster stopCluster clusterEvalQ clusterExport clusterSetRNGStream
#' @importFrom withr with_output_sink
#' @importFrom MASS glm.nb negative.binomial theta.mm
#' @importFrom dplyr rename mutate relocate
#' @importFrom broom.mixed tidy
#' @importFrom purrr imap reduce
#' @importFrom stats predict logLik deviance offset
#' @importFrom geeM geem
#' @importFrom glmmTMB glmmTMB nbinom2
#' @param expr.mat Either a \code{SingleCellExperiment} or \code{Seurat} object from which counts can be extracted, or a matrix of integer-valued counts with genes as rows & cells as columns. Defaults to NULL.
#' @param pt Either the output from \code{\link[slingshot]{SlingshotDataSet}} object from which pseudotime can be generated, or a data.frame containing the pseudotime or latent time estimates for each cell (can be multiple columns / lineages). Defaults to NULL.
#' @param genes A character vector of genes to model. If not provided, defaults to all genes in \code{expr.mat}. Defaults to NULL.
#' @param n.potential.basis.fns (Optional) The maximum number of possible basis functions. See the parameter \code{M} in \code{\link{marge2}}. Defaults to 5.
#' @param size.factor.offset (Optional) An offset to be included in the final model fit. Can be generated easily with \code{\link{createCellOffset}}. Defaults to NULL.
#' @param is.gee Should a GEE framework be used instead of the default GLM? Defaults to FALSE.
#' @param cor.structure If the GEE framework is used, specifies the desired working correlation structure. Must be one of "ar1", "independence", or "exchangeable". Defaults to "ar1".
#' @param id.vec If a GEE or GLMM framework is being used, a vector of subject IDs to use as input to \code{\link[geeM]{geem}} or \code{\link[glmmTMB]{glmmTMB}}. Defaults to NULL.
#' @param is.glmm Should a GLMM framework be used instead of the default GLM? Defaults to FALSE.
#' @param parallel.exec A boolean indicating whether a parallel \code{\link[foreach]{foreach}} loop should be used to generate results more quickly. Defaults to TRUE.
#' @param n.cores (Optional) If running in parallel, how many cores should be used? Defaults to 2.
#' @param approx.knot (Optional) Should the knot space be reduced in order to improve computation time? Defaults to TRUE.
#' @param glmm.adaptive (Optional) Should the basis functions for the GLMM be chosen adaptively? If not, uses 4 evenly spaced knots. Defaults to FALSE.
#' @param verbose (Optional) A boolean indicating whether a progress bar should be printed to the console. Defaults to TRUE.
#' @param random.seed (Optional) The random seed used to initialize RNG streams in parallel. Defaults to 312.
#' @details
#' \itemize{
#' \item If \code{expr.mat} is a \code{Seurat} object, counts will be extracted from the output of \code{\link[SeuratObject]{DefaultAssay}}. If using this functionality, check to ensure the specified assay is correct before running the function. If the input is a \code{SingleCellExperiment} object, the raw counts will be extracted with \code{\link[BiocGenerics]{counts}}.
#' \item If using the GEE or GLMM model architectures, ensure that the observations are sorted by subject ID (this is assumed by the underlying fit implementations). If they are not, the models will error out.
#' }
#' @return A list of lists, where each element is a gene and each gene contains sublists for each element. Each gene-lineage sublist contains a gene name, lineage number, default \code{marge} vs. null model test results, model statistics, and fitted values. Use \code{\link{getResultsDE}} to tidy the results.
#' @seealso \code{\link{getResultsDE}}
#' @seealso \code{\link{testSlope}}
#' @seealso \code{\link{marge2}}
#' @seealso \code{\link[MASS]{glm.nb}}
#' @seealso \code{\link[geeM]{geem}}
#' @seealso \code{\link[glmmTMB]{glmmTMB}}
#' @export
#' @examples
#' data(sim_counts)
#' data(sim_pseudotime)
#' cell_offset <- createCellOffset(sim_counts)
#' scLANE_models <- testDynamic(sim_counts,
#'                              pt = sim_pseudotime,
#'                              size.factor.offset = cell_offset)

testDynamic <- function(expr.mat = NULL,
                        pt = NULL,
                        genes = NULL,
                        n.potential.basis.fns = 5,
                        size.factor.offset = NULL,
                        is.gee = FALSE,
                        cor.structure = "ar1",
                        is.glmm = FALSE,
                        glmm.adaptive = FALSE,
                        id.vec = NULL,
                        parallel.exec = TRUE,
                        n.cores = 2,
                        approx.knot = TRUE,
                        verbose = TRUE,
                        random.seed = 312) {
  # check inputs
  if (is.null(expr.mat) || is.null(pt)) { stop("You forgot some inputs to testDynamic().") }

  # get raw counts from SingleCellExperiment or Seurat object & transpose to cell x gene dense matrix
  if (is.null(genes)) {
    genes <- rownames(expr.mat)
  }
  if (inherits(expr.mat, "SingleCellExperiment")) {
    expr.mat <- BiocGenerics::counts(expr.mat)[genes, ]
  } else if (inherits(expr.mat, "Seurat")) {
    expr.mat <- Seurat::GetAssayData(expr.mat,
                                     slot = "counts",
                                     assay = Seurat::DefaultAssay(expr.mat))
    expr.mat <- expr.mat[genes, ]
  } else if (inherits(expr.mat, "dgCMatrix")) {
    expr.mat <- expr.mat[genes, ]
  } else {
    expr.mat <- expr.mat[genes, ]
  }
  expr.mat <- as.matrix(Matrix::t(expr.mat))  # transpose to dense cell x gene matrix
  if (!(inherits(expr.mat, "matrix") || inherits(expr.mat, "array"))) { stop("Input expr.mat must be coerceable to a matrix of integer counts.") }

  # extract pseudotime dataframe if input is results from Slingshot
  if (inherits(pt, "SlingshotDataSet")) {
    pt <- as.data.frame(slingshot::slingPseudotime(pt))
  }
  if (inherits(pt, "tbl") || inherits(pt, "tbl_df")) {
    pt <- as.data.frame(pt)
  }
  if (!inherits(pt, "data.frame")) { stop("pt must be of class data.frame.") }

  # set pseudotime lineage column names automatically to prevent user error (uses e.g., "Lineage_A", "Lineage_B", etc.)
  n_lineages <- ncol(pt)
  colnames(pt) <- paste0("Lineage_", LETTERS[seq(n_lineages)])

  # ensure subject ID vector meets criteria for GEE / GLMM fitting
  if ((is.gee || is.glmm) && is.null(id.vec)) { stop("You must provide a vector of IDs if you're using GEE / GLMM backends.") }
  if ((is.gee || is.glmm) && is.unsorted(id.vec)) { stop("Your data must be ordered by subject, please do so before running testDynamic() with the GEE / GLMM backends.") }
  cor.structure <- tolower(cor.structure)
  if (is.gee && !(cor.structure %in% c("ar1", "independence", "exchangeable"))) { stop("GEE models require a specified correlation structure.") }

  # set up time tracking
  start_time <- Sys.time()

  # set up progress bar
  if (verbose) {
    withr::with_output_sink(tempfile(), {
      pb <- utils::txtProgressBar(0, length(genes), style = 3)
    })
    progress_fun <- function(n) utils::setTxtProgressBar(pb, n)
    snow_opts <- list(progress = progress_fun)
  } else {
    snow_opts <- list()
  }

  # set up parallel processing
  if (parallel.exec) {
    cl <- parallel::makeCluster(n.cores)
    doSNOW::registerDoSNOW(cl)
    parallel::clusterSetRNGStream(cl, iseed = random.seed)
  } else {
    cl <- foreach::registerDoSEQ()
  }

  # convert dense counts matrix to file-backed matrix
  expr.mat <- bigstatsr::as_FBM(expr.mat,
                                type = "integer",
                                is_read_only = TRUE)

  # build list of objects to prevent from being sent to parallel workers
  necessary_vars <- c("expr.mat", "genes", "pt", "n.potential.basis.fns", "approx.knot", "is.glmm",
                      "n_lineages", "id.vec", "cor.structure", "is.gee", "glmm.adaptive", "size.factor.offset")
  if (any(ls(envir = .GlobalEnv) %in% necessary_vars)) {
    no_export <- c(ls(envir = .GlobalEnv)[-which(ls(envir = .GlobalEnv) %in% necessary_vars)],
                   ls()[-which(ls() %in% necessary_vars)])
  } else {
    no_export <- c(ls(envir = .GlobalEnv),
                   ls()[-which(ls() %in% necessary_vars)])
  }
  no_export <- unique(no_export)
  package_list <- c("glm2", "scLANE", "MASS",  "bigstatsr", "broom.mixed", "dplyr", "stats")
  if (is.gee) {
    package_list <- c(package_list, "geeM")
  } else if (is.glmm) {
    package_list <- c(package_list, "glmmTMB")
  }

  # build models per-lineage per-gene, parallelize over genes
  test_stats <- foreach::foreach(i = seq_along(genes),
                                 .combine = "list",
                                 .multicombine = ifelse(length(genes) > 1, TRUE, FALSE),
                                 .maxcombine = ifelse(length(genes) > 1, length(genes), 2),
                                 .packages = package_list,
                                 .noexport = no_export,
                                 .errorhandling = "pass",
                                 .inorder = TRUE,
                                 .verbose = FALSE,
                                 .options.snow = snow_opts) %dopar% {
    lineage_list <- vector("list", n_lineages)
    for (j in seq(n_lineages)) {
      # pull cells assigned to lineage j
      lineage_cells <- which(!is.na(pt[, j]))

      gene_time_start <- Sys.time()
      # run MARGE model using one of GLM, GEE, or GLMM backends
      if (!is.glmm) {
        marge_mod <- try({
          marge2(X_pred = pt[lineage_cells, j, drop = FALSE],
                 Y = expr.mat[lineage_cells, i],
                 Y.offset = size.factor.offset[lineage_cells],
                 is.gee = is.gee,
                 id.vec = id.vec[lineage_cells],
                 cor.structure = cor.structure,
                 M = n.potential.basis.fns,
                 approx.knot = approx.knot,
                 return.basis = TRUE)
        }, silent = TRUE)
      } else if (is.glmm) {
        marge_mod <- try({
          fitGLMM(X_pred = pt[lineage_cells, j, drop = FALSE],
                  Y = expr.mat[lineage_cells, i],
                  Y.offset = size.factor.offset[lineage_cells],
                  id.vec = id.vec[lineage_cells],
                  M.glm = n.potential.basis.fns,
                  approx.knot = approx.knot,
                  adaptive = glmm.adaptive,
                  return.basis = TRUE)
        }, silent = TRUE)
      }
      gene_time_end <- Sys.time() - gene_time_start
      gene_time_end_numeric <- as.numeric(gene_time_end)
      if (attributes(gene_time_end)$units == "mins") {
        gene_time_end_numeric <- gene_time_end_numeric * 60
      }

      # build formula for null model
      null_mod_df <- data.frame(Y_null = expr.mat[lineage_cells, i],
                                Intercept = 1)
      if (!is.null(id.vec)) {
        null_mod_df <- dplyr::mutate(null_mod_df,
                                     subject = id.vec[lineage_cells])
      }
      if (is.glmm) {
        null_mod_formula <- "Y_null ~ (1 | subject)"
      } else {
        null_mod_formula <- "Y_null ~ -1 + Intercept"
      }
      if (!is.null(size.factor.offset)) {
        null_mod_df <- dplyr::mutate(null_mod_df,
                                     n_offset = size.factor.offset[lineage_cells])
        null_mod_formula <- paste0(null_mod_formula, " + offset(log(1 / n_offset))")
      }
      null_mod_formula <- stats::as.formula(null_mod_formula)

      # fit null model for comparison via Wald or LR test
      if (is.gee) {
        theta_hat <- MASS::theta.mm(y = null_mod_df$Y_null,
                                    mu = mean(null_mod_df$Y_null),
                                    dfr = length(null_mod_df$subject) - 1)
        null_mod <- try({
          geeM::geem(null_mod_formula,
                     id = null_mod_df$subject,
                     data = null_mod_df,
                     family = MASS::negative.binomial(theta_hat),
                     corstr = cor.structure)
        }, silent = TRUE)
      } else if (is.glmm) {
        null_mod <- try({
          glmmTMB::glmmTMB(null_mod_formula,
                           data = null_mod_df,
                           family = glmmTMB::nbinom2(link = "log"),
                           se = TRUE)
        }, silent = TRUE)
      } else {
        null_mod <- try({
          MASS::glm.nb(null_mod_formula,
                       data = null_mod_df,
                       method = "glm.fit2",
                       y = FALSE,
                       model = FALSE,
                       init.theta = 1,
                       link = log)
        }, silent = TRUE)
      }

      # slim down GLM object if not a GEE / GLMM model (which are much smaller for some reason)
      if (!(is.gee || is.glmm)) {
        null_mod <- stripGLM(glm.obj = null_mod)
      }

      # record model fit status
      if (inherits(marge_mod, "try-error")) {
        if (inherits(null_mod, "try-error")) {
          mod_status <- "MARGE model error, null model error"
        } else {
        mod_status <- "MARGE model error, null model OK"
        }
      } else {
        if (inherits(null_mod, "try-error")) {
          mod_status <- "MARGE model OK, null model error"
        } else {
          mod_status <- "MARGE model OK, null model OK"
        }
      }

     # summarize hinge function coefficients
     null_sumy <- pull.null.sumy(null_mod, is.gee, is.glmm)
     marge_sumy <- pull.marge.sumy(marge_mod, is.gee, is.glmm)

     # perform slope test
     marge_slope_df <- createSlopeTestData(marge.model = marge_mod,
                                           pt = pt[lineage_cells, j, drop = FALSE],
                                           is.gee = is.gee,
                                           is.glmm = is.glmm)
     marge_slope_df <-  dplyr::mutate(marge_slope_df,
                                      Gene = genes[i],
                                      Lineage = LETTERS[j],
                                      .before = 1)
     # solve values of slopes across pseudotime intervals
     marge_dynamic_df <- summarizeModel(marge.model = marge_mod, 
                                        pt = pt[lineage_cells, j, drop = FALSE], 
                                        is.glmm = is.glmm)
     if (is.glmm) {
       marge_dynamic_df <- purrr::imap(marge_dynamic_df, \(x, y) {
         data.frame(t(unlist(x))) %>%
         dplyr::mutate(Subject = y, 
                       Gene = genes[i],
                       Lineage = LETTERS[j],
                       .before = 1)
       })
       marge_dynamic_df <- purrr::reduce(marge_dynamic_df, rbind)
     } else {
       marge_dynamic_df <- data.frame(t(unlist(marge_dynamic_df))) %>%
                           dplyr::mutate(Gene = genes[i],
                                         Lineage = LETTERS[j],
                                         .before = 1)
     }

     # format results list
     lineage_list[[j]] <- list(Gene = genes[i],
                               Lineage = LETTERS[j],
                               Test_Stat = NA_real_,
                               Test_Stat_Type = ifelse(is.gee, "Wald", "LRT"),
                               Test_Stat_Note = NA_character_,
                               P_Val = NA_real_,
                               LogLik_MARGE = marge_sumy$ll_marge,
                               LogLik_Null = null_sumy$null_ll,
                               Dev_MARGE = marge_sumy$marge_dev,
                               Dev_Null = null_sumy$null_dev,
                               Model_Status = mod_status,
                               MARGE_Fit_Notes = marge_sumy$marge_fit_notes,
                               Gene_Time = gene_time_end_numeric,
                               MARGE_Summary = marge_sumy$marge_sumy_df,
                               Null_Summary = null_sumy$null_sumy_df,
                               MARGE_Preds = marge_sumy$marge_pred_df,
                               Null_Preds = null_sumy$null_pred_df,
                               MARGE_Slope_Data = marge_slope_df,
                               Gene_Dynamics = marge_dynamic_df)

     # compute test stat using asymptotic Chi-squared approximation
     if (is.gee) {
       test_res <- waldTestGEE(mod.1 = marge_mod, 
                               mod.0 = null_mod, 
                               bias.correct = ifelse(length(unique(id.vec)) < 30, TRUE, FALSE))
     } else {
       test_res <- modelLRT(mod.1 = marge_mod,
                            mod.0 = null_mod,
                            is.glmm = is.glmm)
     }
     # add test stats to result list
     lineage_list[[j]]$Test_Stat <- ifelse(is.gee, test_res$Wald_Stat, test_res$LRT_Stat)
     lineage_list[[j]]$Test_Stat_Note <- test_res$Notes
     lineage_list[[j]]$P_Val <- test_res$P_Val
    }
    names(lineage_list) <- paste0("Lineage_", LETTERS[seq(n_lineages)])
    lineage_list
  }

  # end parallelization & clean up each worker node
  withr::with_output_sink(tempfile(), {
    if (parallel.exec) {
      parallel::clusterEvalQ(cl, expr = {
        gc(verbose = FALSE, full = TRUE)
      })
      parallel::stopCluster(cl)
    }
  })
  # clean up errors w/ proper formatting
  names(test_stats) <- genes
  test_stats <- purrr::imap(test_stats, \(x, y) {
    if (inherits(x, "simpleError")) {
      reformatted_results <- purrr::map(seq(n_lineages), \(z) {
        list(Gene = y,
             Lineage = LETTERS[z],
             Test_Stat = NA_real_,
             Test_Stat_Type = ifelse(is.gee, "Wald", "LRT"),
             Test_Stat_Note = NA_character_,
             P_Val = NA_real_,
             LogLik_MARGE = NA_real_,
             LogLik_Null = NA_real_,
             Dev_MARGE = NA_real_,
             Dev_Null = NA_real_,
             Model_Status = "MARGE model error, null model error",
             Gene_Time = NA_real_,
             MARGE_Fit_Notes = NA_character_,
             MARGE_Summary = NULL,
             Null_Summary = NULL,
             MARGE_Preds = NULL,
             Null_Preds = NULL,
             MARGE_Slope_Data = NULL,
             Gene_Dynamics = NULL)
      })
      names(reformatted_results) <- paste0("Lineage_", LETTERS[seq(n_lineages)])
      return(reformatted_results)
    } else {
      return(x)
    }
  })

  # finalize time tracking
  end_time <- Sys.time()
  total_time <- end_time - start_time
  total_time_units <- attributes(total_time)$units
  total_time_numeric <- as.numeric(total_time)
  time_message <- paste0("scLANE testing completed for ",
                         length(genes),
                         " genes across ",
                         n_lineages,
                         " ",
                         ifelse(n_lineages == 1, "lineage ", "lineages "),
                         "in ",
                         round(total_time_numeric, 3),
                         " ",
                         total_time_units)
  if (verbose) {
    time_message <- paste0("\n", time_message)
  }
  message(time_message)

  # return results
  class(test_stats) <- "scLANE"
  return(test_stats)
}
