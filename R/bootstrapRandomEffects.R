#' Generate bootstrapped confidence intervals for random effects.
#'
#' @name bootstrapRandomEffects
#' @author Jack R. Leary
#' @description This function leverages the parametric bootstrap to generate empirical confidence intervals for the random effects terms of a fitted model.
#' @importFrom stats simulate update quantile
#' @importFrom dplyr mutate with_groups summarise if_else
#' @importFrom glmmTMB ranef
#' @importFrom parallel makeCluster clusterSetRNGStream clusterEvalQ stopCluster
#' @importFrom doSNOW registerDoSNOW
#' @importFrom utils txtProgressBar setTxtProgressBar
#' @importFrom foreach foreach %dopar%
#' @importFrom withr with_output_sink
#' @importFrom tidyr pivot_longer
#' @importFrom purrr reduce
#' @param glmm.mod The output from \code{\link{fitGLMM}}. Defaults to NULL.
#' @param id.vec A vector of subject IDs. Defaults to NULL.
#' @param Y.offset An offset to be included in the final model fit. Defaults to NULL.
#' @param n.boot (Optional) The number of bootstrap resamples to generate. Defaults to 500.
#' @param alpha (Optional) The desired confidence level. Defaults to good old 0.05.
#' @param n.cores (Optional) The number of threads to use in parallel processing of the bootstrap resampling procedure. Defaults to 4.
#' @param random.seed (Optional) The seed used to control stochasticity during bootstrap resampling. Defaults to 312.
#' @param verbose (Optional) A boolean indicating whether a progress bar should be printed to the console. Defaults to TRUE.
#' @return An object of class \code{data.frame} containing the upper and lower quantiles of the per-subject random effects.
#' @seealso \code{\link{fitGLMM}}
#' @seealso \code{\link[glmmTMB]{glmmTMB}}
#' @export
#' @examples
#' data(sim_counts)
#' data(sim_pseudotime)
#' cell_offset <- createCellOffset(sim_counts)
#' glmm_mod <- fitGLMM(X_pred = sim_pseudotime,
#'                     Y = BiocGenerics::counts(sim_counts)[4, ],
#'                     Y.offset = cell_offset,
#'                     id.vec = sim_counts$subject,
#'                     return.basis = TRUE)
#' ranef_sumy <- bootstrapRandomEffects(glmm_mod,
#'                                      id.vec = sim_counts$subject,
#'                                      Y.offset = cell_offset, 
#'                                      n.boot = 100L, 
#'                                      n.cores = 1L)

bootstrapRandomEffects <- function(glmm.mod = NULL,
                                   id.vec = NULL,
                                   Y.offset = NULL,
                                   n.boot = 500L,
                                   alpha = 0.05,
                                   n.cores = 4L,
                                   random.seed = 312, 
                                   verbose = TRUE) {
  # check inputs
  if (is.null(glmm.mod) || is.null(id.vec) || is.null(Y.offset)) { stop("You forgot some arguments to bootstrapRandomEffects().") }
  if (is.null(glmm.mod$basis_mtx)) { stop("fitGLMM() must be run with return.basis = TRUE.") }
  # define simulation function
  simFunction <- function(fitted.model = NULL,
                          basis.mtx = NULL,
                          id.vec = NULL,
                          Y.offset = NULL) {
    sim_data <- stats::simulate(fitted.model$final_mod)
    sim_data <- dplyr::mutate(basis.mtx,
                              Y = sim_data$sim_1,
                              Y.offset = Y.offset,
                              subject = id.vec)
    refit_model <- stats::update(fitted.model$final_mod, data = sim_data)
    sim_ranef <- glmmTMB::ranef(refit_model)$cond$subject %>%
                 dplyr::mutate(subject = rownames(.), .before = 1)
    return(sim_ranef)
  }
  # setup parallelism
  if (n.cores > 1L) {
    cl <- parallel::makeCluster(n.cores)
    doSNOW::registerDoSNOW(cl)
    parallel::clusterSetRNGStream(cl, iseed = random.seed)
  } else {
    cl <- foreach::registerDoSEQ()
  }
  # set up progress bar
  if (verbose) {
    withr::with_output_sink(tempfile(), {
      pb <- utils::txtProgressBar(0, length(n.boot), style = 3)
    })
    progress_fun <- function(n) utils::setTxtProgressBar(pb, n)
    snow_opts <- list(progress = progress_fun)
  } else {
    snow_opts <- list()
  }
  # bootstrap random effects
  ranef_boot <- foreach::foreach(b = seq(n.boot),
                                 .combine = "list",
                                 .multicombine = ifelse(n.boot > 1, TRUE, FALSE),
                                 .maxcombine = ifelse(n.boot > 1, n.boot, 2),
                                 .packages = c("glmmTMB", "dplyr", "stats", "tidyr"),
                                 .errorhandling = "pass",
                                 .inorder = TRUE,
                                 .verbose = FALSE, 
                                 .options.snow = snow_opts) %dopar% {
    ranef_iter_b <- simFunction(glmm.mod,
                                basis.mtx = glmm.mod$basis_mtx,
                                id.vec = id.vec,
                                Y.offset = Y.offset) %>%
                    tidyr::pivot_longer(cols = !subject,
                                        names_to = "term",
                                        values_to = "effect") %>%
                    dplyr::mutate(iter = b, .before = 1)
    return(ranef_iter_b)
  }
  # end parallelization & clean up each worker node
  withr::with_output_sink(tempfile(), {
    if (n.cores > 1L) {
      parallel::clusterEvalQ(cl, expr = {
        gc(verbose = FALSE, full = TRUE)
      })
      parallel::stopCluster(cl)
    }
  })
  # summarize bootstrap resample
  lower_bound <- alpha / 2
  upper_bound <- 1 - alpha / 2
  ranef_sumy <- purrr::reduce(ranef_boot, rbind) %>%
                dplyr::with_groups(c(subject, term),
                                   dplyr::summarise,
                                   QL = stats::quantile(effect, probs = lower_bound),
                                   QU = stats::quantile(effect, probs = upper_bound)) %>%
                dplyr::mutate(term = dplyr::if_else(term == "(Intercept)", "Intercept", term))
  return(ranef_sumy)
}
