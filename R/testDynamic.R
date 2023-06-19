#' Test whether a gene is dynamic over pseudotime.
#'
#' @name testDynamic
#' @author Jack Leary
#' @description This function tests whether a NB \code{marge} model is better than a null (intercept-only) NB GLM using the Likelihood Ratio Test. In effect, the test tells us whether a gene's expression changes (in any way) over pseudotime.
#' @import glm2
#' @import magrittr
#' @importFrom bigstatsr as_FBM
#' @importFrom foreach foreach %dopar% registerDoSEQ
#' @importFrom doParallel registerDoParallel
#' @importFrom parallel makeCluster detectCores stopCluster clusterEvalQ clusterExport clusterSetRNGStream
#' @importFrom withr with_output_sink
#' @importFrom MASS glm.nb negative.binomial theta.mm
#' @importFrom dplyr rename mutate relocate
#' @importFrom broom.mixed tidy
#' @importFrom stats predict logLik deviance offset
#' @importFrom geeM geem
#' @importFrom glmmTMB glmmTMB nbinom2
#' @param expr.mat Either a \code{SingleCellExperiment} or \code{Seurat} object from which counts can be extracted, or a dense matrix of integer-valued counts. Defaults to NULL.
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
#' @param track.time (Optional) A boolean indicating whether the amount of time the function takes to run should be tracked and printed to the console. Useful for debugging. Defaults to FALSE.
#' @param log.file (Optional) A string indicating a \code{.txt} file to which iteration tracking should be written. Can be useful for debugging. Defaults to NULL.
#' @param log.iter (Optional) If logging is enabled, how often should iterations be printed to the logfile. Defaults to 100.
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
#' \dontrun{
#' testDynamic(expr.mat = raw_counts,
#'             pt = pseudotime_df,
#'             parallel.exec = TRUE)
#' testDynamic(expr.mat = sce_obj,
#'             pt = slingshot_obj,
#'             size.factor.offset = sizeFactors(sce_obj),
#'             genes = rownames(sce_obj)[1:100])
#' testDynamic(expr.mat = raw_counts,
#'             pt = pseudotime_df,
#'             parallel.exec = TRUE,
#'             n.cores = 8,
#'             n.potential.basis.fns = 7)
#' testDynamic(expr.mat = raw_counts,
#'             pt = pseudotime_df,
#'             is.gee = TRUE,
#'             id.vec = my_subject_ids,
#'             cor.structure = "ar1",
#'             parallel.exec = TRUE,
#'             n.cores = 8,
#'             n.potential.basis.fns = 7)
#' testDynamic(expr.mat = raw_counts,
#'             pt = pseudotime_df,
#'             parallel.exec = TRUE,
#'             n.cores = 8,
#'             is.glmm = TRUE,
#'             id.vec = my_subject_ids,
#'             log.file = "scLANE_log.txt",
#'             log.iter = 10)
#' }

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
                        track.time = FALSE,
                        log.file = NULL,
                        log.iter = 100) {
  # check inputs
  if (is.null(expr.mat) || is.null(pt)) { stop("You forgot some inputs to testDynamic().") }
  # get raw counts from SingleCellExperiment or Seurat object & transpose to cell x gene dense matrix
  if (inherits(expr.mat, "SingleCellExperiment")) {
    expr.mat <- as.matrix(t(BiocGenerics::counts(expr.mat)))
  } else if (inherits(expr.mat, "Seurat")) {
    expr.mat <- as.matrix(t(Seurat::GetAssayData(expr.mat,
                                                 slot = "counts",
                                                 assay = Seurat::DefaultAssay(expr.mat))))
  }
  if (!(inherits(expr.mat, "matrix") || inherits(expr.mat, "array"))) { stop("Input expr.mat must be coerceable to a matrix of integer counts.") }
  # extract pseudotime dataframe if input is results from Slingshot
  if (inherits(pt, "SlingshotDataSet")) {
    pt <- as.data.frame(slingshot::slingPseudotime(pt))
  }
  if (!inherits(pt, "data.frame")) { stop("pt must be of class data.frame.") }
  # set pseudotime lineage column names automatically to prevent user error (uses e.g., "Lineage_A", "Lineage_B", etc.)
  n_lineages <- ncol(pt)
  colnames(pt) <- paste0("Lineage_", LETTERS[1:n_lineages])
  # ensure subject ID vector meets criteria for GEE / GLMM fitting
  if ((is.gee || is.glmm) && is.null(id.vec)) { stop("You must provide a vector of IDs if you're using GEE / GLMM backends.") }
  if ((is.gee || is.glmm) && is.unsorted(id.vec)) { stop("Your data must be ordered by subject, please do so before running testDynamic() with the GEE / GLMM backends.") }
  cor.structure <- tolower(cor.structure)
  if (is.gee && !(cor.structure %in% c("ar1", "independence", "exchangeable"))) { stop("GEE models require a specified correlation structure.") }
  # fit models for all genes if otherwise unspecified
  if (is.null(genes)) {
    genes <- colnames(expr.mat)
  }
  # set up parallel execution & time tracking
  if (track.time) {
    start_time <- Sys.time()
  }
  if (parallel.exec) {
    cl <- parallel::makeCluster(n.cores)
    doParallel::registerDoParallel(cl)
    parallel::clusterSetRNGStream(cl, iseed = 312)
  } else {
    cl <- foreach::registerDoSEQ()
    set.seed(312)
  }
  # convert dense counts matrix to file-backed matrix
  expr.mat <- bigstatsr::as_FBM(expr.mat,
                                type = "integer",
                                is_read_only = TRUE)
  # set up logging to .txt file if desired
  print_nums <- seq(0, length(genes), log.iter)[-1]
  if (!is.null(log.file)) {
    if (!substr(log.file, nchar(log.file) - 3, nchar(log.file)) == ".txt") {
      log.file <- paste0(log.file, ".txt")
    }
  }
  # build list of objects to prevent from being sent to parallel workers
  necessary_vars <- c("expr.mat", "genes", "pt", "n.potential.basis.fns", "approx.knot", "is.glmm", "print_nums",
                      "n_lineages", "id.vec", "cor.structure", "is.gee", "log.file", "log.iter", "glmm.adaptive", "size.factor.offset")
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
  }
  if (is.glmm) {
    package_list <- c(package_list, "glmmTMB")
  }
  # build models per-lineage per-gene, parallelize over genes
  test_stats <- foreach::foreach(i = seq_along(genes),
                                 .combine = "list",
                                 .multicombine = ifelse(length(genes) > 1, TRUE, FALSE),
                                 .maxcombine = ifelse(length(genes) > 1, length(genes), 2),
                                 .packages = package_list,
                                 .noexport = no_export,
                                 .verbose = FALSE) %dopar% {
    if (!is.null(log.file) && i %in% print_nums) {
      withr::with_output_sink(log.file,
                              code = {
                                cat(paste(Sys.time(), "- Starting iteration:", 20, "\n"))
                              },
                              append = TRUE)
    }
    lineage_list <- vector("list", n_lineages)
    for (j in seq(n_lineages)) {
      lineage_cells <- which(!is.na(pt[, j]))
      # run original MARGE model using one of GLM, GEE, or GLMM backends
      if (is.gee) {
        marge_mod <- try(
          { scLANE::marge2(X_pred = pt[lineage_cells, j, drop = FALSE],
                           Y = expr.mat[lineage_cells, i],
                           Y.offset = size.factor.offset[lineage_cells],
                           is.gee = is.gee,
                           id.vec = id.vec[lineage_cells],
                           cor.structure = cor.structure,
                           M = n.potential.basis.fns,
                           approx.knot = approx.knot) },
          silent = TRUE
        )
      } else if (is.glmm) {
        if (glmm.adaptive) {
          marge_mod <- try(
            { scLANE::fitGLMM(X_pred = pt[lineage_cells, j, drop = FALSE],
                              Y = expr.mat[lineage_cells, i],
                              Y.offset = size.factor.offset[lineage_cells],
                              id.vec = id.vec[lineage_cells],
                              M.glm = n.potential.basis.fns,
                              approx.knot = approx.knot,
                              adaptive = glmm.adaptive) },
            silent = TRUE
          )
        } else {
          marge_mod <- try(
            { scLANE::fitGLMM(X_pred = pt[lineage_cells, j, drop = FALSE],
                              Y = expr.mat[lineage_cells, i],
                              Y.offset = size.factor.offset[lineage_cells],
                              id.vec = id.vec[lineage_cells],
                              M.glm = n.potential.basis.fns,
                              approx.knot = approx.knot,
                              adaptive = glmm.adaptive,
                              return.basis = TRUE) },
            silent = TRUE
          )
        }
      } else {
        marge_mod <- try(
          { scLANE::marge2(X_pred = pt[lineage_cells, j, drop = FALSE],
                           Y = expr.mat[lineage_cells, i, drop = FALSE],
                           Y.offset = size.factor.offset[lineage_cells],
                           M = n.potential.basis.fns,
                           approx.knot = approx.knot,
                           return.basis = TRUE) },
          silent = TRUE
        )
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
        null_mod <- try(
          { geeM::geem(null_mod_formula,
                       id = null_mod_df$subject,
                       data = null_mod_df,
                       family = MASS::negative.binomial(theta_hat),
                       corstr = cor.structure) },
          silent = TRUE
        )
      } else if (is.glmm) {
        null_mod <- try(
          { glmmTMB::glmmTMB(null_mod_formula,
                             data = null_mod_df,
                             family = glmmTMB::nbinom2(link = "log"),
                             se = TRUE) },
          silent = TRUE
        )
      } else {
        null_mod <- try(
          { MASS::glm.nb(null_mod_formula,
                         data = null_mod_df,
                         method = "glm.fit2",
                         y = FALSE,
                         model = FALSE,
                         init.theta = 1,
                         link = log) },
          silent = TRUE
        )
      }
      # slim down GLM object if not a GEE / GLMM model (which are much smaller for some reason)
      if (!(is.gee || is.glmm) && !inherits(null_mod, "try-error")) {
        null_mod <- scLANE::stripGLM(glm.obj = null_mod)
      }
      # prepare results if there were errors in both null and MARGE models
      if (inherits(marge_mod, "try-error") && inherits(null_mod, "try-error")) {
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
                                  MARGE_Fit_Notes = NA_character_,
                                  MARGE_Summary = NA,
                                  Null_Summary = NA,
                                  MARGE_Preds = NA,
                                  Null_Preds = NA,
                                  MARGE_Slope_Data = slope_data_error)
      # prepare results if there were errors in MARGE model but not in null model
      } else if (inherits(marge_mod, "try-error") && !inherits(null_mod, "try-error")) {
        mod_status <- "MARGE model error, null model OK"
        # generate empty dataframe for slope test
        slope_data_error <- data.frame(Gene = genes[i],
                                       Lineage = LETTERS[j],
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
        } else if (is.glmm) {
          null_sumy_df <- try({
            null_glmm_summary <- as.data.frame(broom.mixed::tidy(null_mod, effects = "fixed"))
            data.frame(term = null_glmm_summary$term,
                       estimate = null_glmm_summary$estimate[1],
                       `std.error` = null_glmm_summary$std.error[1],
                       statistic = null_glmm_summary$statistic[1],
                       `p.value` = null_glmm_summary$p.value[1])
          }, silent = TRUE)
          null_pred_df <- try({
            data.frame(predict(null_mod, type = "link", se.fit = TRUE)[1:2]) %>%
              dplyr::rename(null_link_fit = fit, null_link_se = se.fit)
          }, silent = TRUE)
        } else {
          null_sumy_df <- try({
            as.data.frame(broom.mixed::tidy(null_mod)) # saves a few bytes by converting from tibble
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
                                  Dev_Null = ifelse((is.gee || is.glmm), NA_real_, stats::deviance(null_mod)),
                                  Model_Status = mod_status,
                                  MARGE_Fit_Notes = NA_character_,
                                  MARGE_Summary = NA,
                                  Null_Summary = null_sumy_df,
                                  MARGE_Preds = NA,
                                  Null_Preds = null_pred_df,
                                  MARGE_Slope_Data = slope_data_error)
      # prepare results if there were errors in null model but not in marge model
      } else if (!inherits(marge_mod, "try-error") && inherits(null_mod, "try-error")) {
        mod_status <- "MARGE model OK, null model error"
        # generate marge model summary stat table & fitted values w/ standard errors (need to do so manually for GEE)
        if (is.gee) {
          robust_vcov_mat <- as.matrix(marge_mod$final_mod$var)
          marge_pred_df <- try({
            data.frame(marge_link_fit = predict(marge_mod$final_mod),
                       marge_link_se = sqrt(apply((tcrossprod(marge_mod$final_mod$X, robust_vcov_mat)) * marge_mod$final_mod$X, 1, sum)))
          }, silent = TRUE)
          marge_sumy_df <- try({
            marge_gee_summary <- summary(marge_mod$final_mod)
            data.frame(term = marge_gee_summary$coefnames,
                       estimate = unname(marge_gee_summary$beta),
                       `std.error` = unname(marge_gee_summary$se.robust),
                       statistic = unname(marge_gee_summary$wald.test),
                       `p.value` = unname(marge_gee_summary$p))
          }, silent = TRUE)
        } else if (is.glmm) {
          marge_pred_df <- try({
            data.frame(predict(marge_mod$final_mod, type = "link", se.fit = TRUE)[1:2]) %>%
              dplyr::rename(marge_link_fit = fit, marge_link_se = se.fit)
          }, silent = TRUE)
          marge_sumy_df <- try({
            marge_glmm_summary <- broom.mixed::tidy(marge_mod$final_mod, effects = "fixed")
            data.frame(term = marge_mod$marge_coef_names,  # use marge-style hinge function names instead of X1, X2, ..., XP (makes interpretability easier)
                       estimate = marge_glmm_summary$estimate,
                       `std.error` = marge_glmm_summary$std.error,
                       statistic = marge_glmm_summary$statistic,
                       `p.value` = marge_glmm_summary$p.value)
          }, silent = TRUE)
        } else {
          marge_pred_df <- try({
            data.frame(stats::predict(marge_mod$final_mod, type = "link", se.fit = TRUE)[1:2]) %>%
              dplyr::rename(marge_link_fit = fit, marge_link_se = se.fit)
          }, silent = TRUE)
          marge_sumy_df <- try({
            as.data.frame(broom.mixed::tidy(marge_mod$final_mod)) %>%
              lapply(unname) %>%
              as.data.frame()
          }, silent = TRUE)
        }
        # create slope test input table
        marge_slope_df <- scLANE:::createSlopeTestData(marge.model = marge_mod,
                                                       pt = pt[lineage_cells, j, drop = FALSE],
                                                       is.gee = is.gee,
                                                       is.glmm = is.glmm) %>%
                          dplyr::mutate(Gene = genes[i], Lineage = LETTERS[j]) %>%
                          dplyr::relocate(Gene, Lineage)
        # get log-likelihood for GLMM or GLM cases
        if (is.glmm) {
          ll_marge <- -marge_mod$final_mod$fit$objective
        } else if (!(is.gee || is.glmm)) {
          ll_marge <- as.numeric(stats::logLik(marge_mod$final_mod))
        } else {
          ll_marge <- NA_real_
        }
        # check positive-definiteness of Hessian for GLMM -- might have an effect on LRT stat / accompanying p-value
        if (is.glmm) {
          if (!marge_mod$final_mod$sdr$pdHess) {
            marge_fit_notes <- "Non-positive definite Hessian in GLMM, probably due to shallow log-likelihood. Be careful!"
          } else {
            marge_fit_notes <- NA_character_
          }
        } else {
          marge_fit_notes <- NA_character_
        }
        # create lineage result list
        lineage_list[[j]] <- list(Gene = genes[i],
                                  Lineage = LETTERS[j],
                                  Test_Stat = NA_real_,
                                  Test_Stat_Type = ifelse(is.gee, "Wald", "LRT"),
                                  Test_Stat_Note = NA_character_,
                                  P_Val = NA_real_,
                                  LogLik_MARGE = ll_marge,
                                  LogLik_Null = NA_real_,
                                  Dev_MARGE = ifelse((is.gee || is.glmm), NA_real_, stats::deviance(marge_mod$final_mod)),
                                  Dev_Null = NA_real_,
                                  Model_Status = mod_status,
                                  MARGE_Fit_Notes = marge_fit_notes,
                                  MARGE_Summary = marge_sumy_df,
                                  Null_Summary = NA,
                                  MARGE_Preds = marge_pred_df,
                                  Null_Preds = NA,
                                  MARGE_Slope_Data = marge_slope_df)

      # prepare results if neither model had errors
      } else if (!inherits(marge_mod, "try-error") && !inherits(null_mod, "try-error")) {
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
        } else if (is.glmm) {
          null_sumy_df <- try({
            null_glmm_summary <- as.data.frame(broom.mixed::tidy(null_mod, effects = "fixed"))
            data.frame(term = null_glmm_summary$term,
                       estimate = null_glmm_summary$estimate[1],
                       `std.error` = null_glmm_summary$std.error[1],
                       statistic = null_glmm_summary$statistic[1],
                       `p.value` = null_glmm_summary$p.value[1])
          }, silent = TRUE)
          null_pred_df <- try({
            data.frame(predict(null_mod, type = "link", se.fit = TRUE)[1:2]) %>%
              dplyr::rename(null_link_fit = fit, null_link_se = se.fit)
          }, silent = TRUE)
        } else {
          null_sumy_df <- try({
            as.data.frame(broom.mixed::tidy(null_mod)) # saves a few bytes by converting from tibble
          }, silent = TRUE)
          null_pred_df <- try({
            data.frame(stats::predict(null_mod, type = "link", se.fit = TRUE)[1:2]) %>%
              dplyr::rename(null_link_fit = fit, null_link_se = se.fit)
          }, silent = TRUE)
        }
        # generate marge model summary stat table & fitted values w/ standard errors (need to do so manually for GEE)
        if (is.gee) {
          robust_vcov_mat <- as.matrix(marge_mod$final_mod$var)
          marge_pred_df <- try({
            data.frame(marge_link_fit = predict(marge_mod$final_mod),
                       marge_link_se = sqrt(apply((tcrossprod(marge_mod$final_mod$X, robust_vcov_mat)) * marge_mod$final_mod$X, 1, sum)))
          }, silent = TRUE)
          marge_sumy_df <- try({
            marge_gee_summary <- summary(marge_mod$final_mod)
            data.frame(term = marge_gee_summary$coefnames,
                       estimate = unname(marge_gee_summary$beta),
                       `std.error` = unname(marge_gee_summary$se.robust),
                       statistic = unname(marge_gee_summary$wald.test),
                       `p.value` = unname(marge_gee_summary$p))
          }, silent = TRUE)
        } else if (is.glmm) {
          marge_pred_df <- try({
            data.frame(predict(marge_mod$final_mod, type = "link", se.fit = TRUE)[1:2]) %>%
              dplyr::rename(marge_link_fit = fit, marge_link_se = se.fit)
          }, silent = TRUE)
          marge_sumy_df <- try({
            marge_glmm_summary <- broom.mixed::tidy(marge_mod$final_mod, effects = "fixed")
            data.frame(term = marge_mod$coef_names,  # use hinge function names instead of X1, X2, ..., XP (makes interpretability easier)
                       estimate = marge_glmm_summary$estimate,
                       `std.error` = marge_glmm_summary$std.error,
                       statistic = marge_glmm_summary$statistic,
                       `p.value` = marge_glmm_summary$p.value)
          }, silent = TRUE)
        } else {
          marge_pred_df <- try({
            data.frame(stats::predict(marge_mod$final_mod, type = "link", se.fit = TRUE)[1:2]) %>%
              dplyr::rename(marge_link_fit = fit, marge_link_se = se.fit)
          }, silent = TRUE)
          marge_sumy_df <- try({
            as.data.frame(broom.mixed::tidy(marge_mod$final_mod)) %>%
              lapply(unname) %>%
              as.data.frame()
          }, silent = TRUE)
        }
        # create slope test input table
        marge_slope_df <- scLANE:::createSlopeTestData(marge.model = marge_mod,
                                                       pt = pt[lineage_cells, j, drop = FALSE],
                                                       is.gee = is.gee,
                                                       is.glmm = is.glmm) %>%
                          dplyr::mutate(Gene = genes[i], Lineage = LETTERS[j]) %>%
                          dplyr::relocate(Gene, Lineage)
        # compute LRT stat using asymptotic Chi-squared approximation
        if (is.gee) {
          test_res <- scLANE::waldTestGEE(mod.1 = marge_mod$final_mod, mod.0 = null_mod)
        } else {
          test_res <- scLANE::modelLRT(mod.1 = marge_mod$final_mod,
                                       mod.0 = null_mod,
                                       is.glmm = is.glmm)
        }
        # check positive-definiteness of Hessian for GLMM -- might have an effect on LRT stat / accompanying p-value
        if (is.glmm) {
          if (!marge_mod$final_mod$sdr$pdHess) {
            marge_fit_notes <- "Non-positive definite Hessian in GLMM, probably due to shallow log-likelihood. Be careful!"
          } else {
            marge_fit_notes <- NA_character_
          }
        } else {
          marge_fit_notes <- NA_character_
        }
        # create lineage result list
        lineage_list[[j]] <- list(Gene = genes[i],
                                  Lineage = LETTERS[j],
                                  Test_Stat = ifelse(is.gee, test_res$Wald_Stat, test_res$LRT_Stat),
                                  Test_Stat_Type = ifelse(is.gee, "Wald", "LRT"),
                                  Test_Stat_Note = test_res$Notes,
                                  P_Val = test_res$P_Val,
                                  LogLik_MARGE = ifelse(is.gee, NA_real_, test_res$Alt_Mod_LL),
                                  LogLik_Null = ifelse(is.gee, NA_real_, test_res$Null_Mod_LL),
                                  Dev_MARGE = ifelse((is.gee || is.glmm), NA_real_, stats::deviance(marge_mod$final_mod)),
                                  Dev_Null = ifelse((is.gee || is.glmm), NA_real_, stats::deviance(null_mod)),
                                  Model_Status = mod_status,
                                  MARGE_Fit_Notes = marge_fit_notes,
                                  MARGE_Summary = marge_sumy_df,
                                  Null_Summary = null_sumy_df,
                                  MARGE_Preds = marge_pred_df,
                                  Null_Preds = null_pred_df,
                                  MARGE_Slope_Data = marge_slope_df)
      } else {
        stop(paste0("Conditions for marge or null model fits not met for gene ",
                    genes[i],
                    " on lineage ",
                    j))
      }
    }
    names(lineage_list) <- paste0("Lineage_", LETTERS[1:n_lineages])
    lineage_list
  }
  # end parallelization & clean up
  withr::with_output_sink(tempfile(), {
    if (parallel.exec) {
      parallel::clusterEvalQ(cl, expr = {
        rm(list = ls(all.names = TRUE)); gc(verbose = FALSE, full = TRUE)
      })
      parallel::stopCluster(cl)
    }
    rm(cl)
    gc(verbose = FALSE, full = TRUE)
  })
  # prepare results
  names(test_stats) <- genes
  if (track.time) {
    end_time <- Sys.time()
    total_time <- end_time - start_time
    total_time_units <- attributes(total_time)$units
    total_time_numeric <- as.numeric(total_time)
    print(paste0("testDynamic evaluated ",
                 length(genes),
                 " genes with ",
                  n_lineages,
                 " ",
                 ifelse(n_lineages == 1, "lineage ", "lineages "),
                 "apiece in ",
                 round(total_time_numeric, 3),
                 " ",
                 total_time_units))
  }
  class(test_stats) <- "scLANE"
  return(test_stats)
}
