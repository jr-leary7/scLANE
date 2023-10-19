#' Represent a \code{marge} model as a series of piecewise equations.
#'
#' @name summarizeModel
#' @author Jack Leary
#' @author Rhonda Bacher
#' @import magrittr
#' @importFrom stats vcov
#' @description This function summarizes the model for each gene and allows for quantiative interpretation of fitted gene dynamics.
#' @param marge.model The fitted model output from \code{\link{marge2}} (this function is internally called by \code{\link{testDynamic}}). Defaults to NULL.
#' @param pt The predictor matrix of pseudotime. Defaults to NULL.
#' @return A data.frame of the model coefficients, cutpoint intervals, and formatted equations.
#' @seealso \code{\link{marge2}}
#' @examples
#' \dontrun{
#' summarizeModel(marge.model = marge_mod, pt = pt_df)
#' }

summarizeModel <- function(marge.model = NULL, pt = NULL) {
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
    coef_df <- data.frame(coef_name = names(coef(marge.model$final_mod)),
                          coef_value = unname(coef(marge.model$final_mod)),
                          coef_var = unname(diag(coef_vcov)))

    coef_df <- coef_df[-which(coef_df$coef_name == "Intercept"), ]

    coef_df <- cbind(coef_df, extractBreakpoints(marge.model))

    coef_df <- coef_df[order(coef_df$Breakpoint, coef_df$Direction), ]

    MIN <- min(pt[, 1])
    MAX <- max(pt[, 1])

    coef_ranges <- mapply(function(x, y) {
      if (x == "Right") {
        c(y, MAX)
      } else if (x == "Left") {
        c(MIN, y)
      }
    }, coef_df$Direction, coef_df$Breakpoint)

    num_segments <- length(unique(coef_df$Breakpoint)) + 1
    mod_seg <- list()
    for (i in seq(num_segments)) {
      if (i == 1) {
        mod_seg[[i]] <- c(MIN, coef_df$Breakpoint[i])
      } else if (i == num_segments) {
        mod_seg[[i]] <- c(coef_df$Breakpoint[i - 1], MAX)
      } else {
        mod_seg[[i]] <- c(coef_df$Breakpoint[i - 1], coef_df$Breakpoint[i])
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
        if (coef_df$Direction[x] == "Right") {
          temp_slp <- coef_df$coef_value[x]
        } else if (coef_df$Direction[x] == "Left") {
          temp_slp <- -1 * coef_df$coef_value[x]
        }
      } else {
        to_rev <- which(coef_df$Direction[x] == "Left")
        coef_df$coef_value[x][to_rev] <- -1 * coef_df$coef_value[x][to_rev]
        temp_slp <- sum(coef_df$coef_value[x])
      }
      return(temp_slp)
    })
    # combine results
    seg_slopes <- do.call(c, seg_slopes)
    # create discretized trend summary
    seg_trends <- ifelse(seg_slopes > 0, 1, -1)
    seg_trends[seg_slopes == 0] <- 0
    # prepare results
    mod_summ <- list(Breakpoint = unique(coef_df$Breakpoint),
                     Slope.Segment = seg_slopes,
                     Trend.Segment = seg_trends)

  }
  return(mod_summ)
}
