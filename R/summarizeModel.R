#' Represent a \code{marge} model as a series of piecewise equations.
#'
#' @name summarizeModel
#' @author Jack Leary
#' @author Rhonda Bacher
#' @import magrittr
#' @importFrom stats vcov
#' @importFrom dplyr mutate select
#' @importFrom purrr map
#' @importFrom tidyr pivot_longer
#' @importFrom tidyselect everything
#' @description This function summarizes the model for each gene and allows for quantitative interpretation of fitted gene dynamics.
#' @param marge.model The fitted model output from \code{\link{marge2}} (this function is internally called by \code{\link{testDynamic}}). Defaults to NULL.
#' @param pt The predictor matrix of pseudotime. Defaults to NULL.
#' @param is.glmm Flag specifying whether or not data were generated using GLMM mode. Defaults to FALSE.
#' @return A data.frame of the model coefficients, cutpoint intervals, and formatted equations.
#' @seealso \code{\link{marge2}}
#' @export
#' @examples
#' data(sim_counts)
#' data(sim_pseudotime)
#' cell_offset <- createCellOffset(sim_counts)
#' marge_model <- marge2(sim_pseudotime,
#'                       Y = BiocGenerics::counts(sim_counts)[4, ],
#'                       Y.offset = cell_offset)
#' model_summary <- summarizeModel(marge.model = marge_model, pt = sim_pseudotime)

summarizeModel <- function(marge.model = NULL, 
                           pt = NULL, 
                           is.glmm = FALSE) {
  # check inputs
  if (is.null(marge.model)) { stop("Please provide a non-NULL input argument for marge.model.") }
  if (is.null(pt)) { stop("Please provide a non-NULL input argument for pt.") }
  if (inherits(marge.model, "try-error")) {
    mod_summ <- list(Breakpoint = NA_real_,
                     Slope.Segment = NA_real_,
                     Trend.Segment = NA_real_)
  } else {
    # extract model equation, slopes, & covariances
    if (marge.model$model_type == "GEE") {
      coef_vcov <- as.matrix(marge.model$final_mod$var)
    } else {
      coef_vcov <- stats::vcov(marge.model$final_mod)
    }
    if (length( stats::coef(marge.model$final_mod)) > 1) {
    if (is.glmm) {
      coef_df_full <- stats::coef(marge.model$final_mod)[[1]][[1]] %>% 
                      dplyr::mutate(subject = rownames(.), .before = 1)
      coef_df_list <- split(coef_df_full, coef_df_full$subject) %>% 
                      purrr::map(\(x) dplyr::select(x, -subject)) %>% 
                      purrr::map(\(x) { colnames(x)[colnames(x) == "(Intercept)"] <- "Intercept"; x }) %>% 
                      purrr::map(\(x) {
                        tidyr::pivot_longer(x, 
                                            cols = tidyselect::everything(), 
                                            names_to = "coef_name", 
                                            values_to = "coef_value") %>% 
                        dplyr::mutate(coef_var = NA_real_) %>% 
                        as.data.frame()
                      })
    } else {
      coef_df <- data.frame(coef_name = names(coef(marge.model$final_mod)),
                            coef_value = unname(coef(marge.model$final_mod)),
                            coef_var = unname(diag(coef_vcov)))
      coef_df_list <- list(coef_df)
    }
    coef_df_list <- purrr::map(coef_df_list, \(x) x[-which(x$coef_name == "Intercept"), ]) %>% 
                    purrr::map(\(x) cbind(x, extractBreakpoints(marge.model))) %>% 
                    purrr::map(\(x) x[order(x$Breakpoint, x$Direction), ])
    # summarize hinge coefficients
    MIN <- min(pt[, 1])
    MAX <- max(pt[, 1])
    mod_summ_list <- purrr::map(coef_df_list, \(z) {
      coef_ranges <- mapply(function(x, y) {
        if (x == "Right") {
          c(y, MAX)
        } else if (x == "Left") {
          c(MIN, y)
        }
      }, z$Direction, z$Breakpoint)
      num_segments <- length(unique(z$Breakpoint)) + 1
      mod_seg <- list()
      for (i in seq(num_segments)) {
        if (i == 1) {
          mod_seg[[i]] <- c(MIN, z$Breakpoint[i])
        } else if (i == num_segments) {
          mod_seg[[i]] <- c(z$Breakpoint[i - 1], MAX)
        } else {
          mod_seg[[i]] <- c(z$Breakpoint[i - 1], z$Breakpoint[i])
        }
      }
      coef_ranges <- t(coef_ranges)
      mod_seg_overlaps <- lapply(mod_seg, function(x) {
        overlap_ind <- c()
        for (i in seq(nrow(coef_ranges))) {
          if (x[1] < coef_ranges[i, 2] && x[2] > coef_ranges[i, 1]) {
            overlap_ind <- c(overlap_ind, i)
          }
        }
        return(overlap_ind)
      })
      
      seg_slopes <- lapply(mod_seg_overlaps, function(x) {
        if (length(x) == 1) {
          if (z$Direction[x] == "Right") {
            temp_slp <- z$coef_value[x]
          } else if (z$Direction[x] == "Left") {
            temp_slp <- -1 * z$coef_value[x]
          }
        } else {
          to_rev <- which(z$Direction[x] == "Left")
          z$coef_value[x][to_rev] <- -1 * z$coef_value[x][to_rev]
          temp_slp <- sum(z$coef_value[x])
        }
        return(temp_slp)
      })
      # combine results
      seg_slopes <- do.call(c, seg_slopes)
      # create discretized trend summary
      seg_trends <- ifelse(seg_slopes > 0, 1, -1)
      seg_trends[seg_slopes == 0] <- 0
      # prepare results
      mod_summ <- list(Breakpoint = unique(z$Breakpoint),
                       Slope.Segment = seg_slopes,
                       Trend.Segment = seg_trends)
      return(mod_summ)
    })
    if (!is.glmm) {
      mod_summ <- mod_summ_list[[1]]
    } else {
      mod_summ <- mod_summ_list
    }
    }  else {
    ## If intercept only model:
    mod_summ <- list(Breakpoint = NA,
                       Slope.Segment = 0,
                       Trend.Segment = 0)
    }
  }
  return(mod_summ)
}
