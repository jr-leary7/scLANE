#' Plot results of \code{marge} and other models using \code{ggplot2}.
#'
#' @name plotModels
#' @author Jack Leary
#' @import magrittr
#' @importFrom stats qnorm predict
#' @importFrom dplyr mutate rename bind_cols pull
#' @importFrom ggplot2 ggplot aes geom_point geom_line scale_y_continuous labs theme_classic theme element_text geom_ribbon scale_color_manual scale_fill_manual geom_vline
#' @importFrom scales comma_format
#' @description This function visualizes the fitted values of a \code{marge} model over the expression and pseudotime values of each cell.
#' @param gene.counts A vector of integer expression values for each cell. Defaults to NULL.
#' @param pt A data.frame of pseudotime values for each cell. Defaults to NULL.
#' @param gene The name of the gene that's being analyzed. Used as the title of the \code{ggplot} object. Defaults to NULL.
#' @param marge.preds (Optional) The fitted values for a \code{marge} model as returned by \code{\link{testDynamic}}. Defaults to NULL.
#' @param marge.ci A boolean indicating whether a (\eqn{1 - \alpha})\% CI ribbon should be drawn for the \code{marge} model. Defaults to TRUE.
#' @param plot.breakpoints A boolean indicating whether or not vertical lines for each changepoint in the \code{marge} model should be plotted. Defaults to FALSE.
#' @param null.preds (Optional) The fitted values for a null (intercept-only) model as returned by \code{\link{testDynamic}}. Defaults to NULL.
#' @param null.ci A boolean indicating whether a (\eqn{1 - \alpha})\% CI ribbon should be drawn for the null model. Defaults to FALSE.
#' @param marge.ci A boolean indicating whether a (\eqn{1 - \alpha})\% CI ribbon should be drawn for the optional intercept-only model. Defaults to FALSE
#' @param plot.gam (Optional) Should the fitted values from an NB GAM be plotted? Defaults to FALSE.
#' @param gam.ci A boolean indicating whether a (\eqn{1 - \alpha})\% CI ribbon should be drawn for the optional \code{gam} model. Defaults to TRUE.
#' @param plot.glm (Optional) Should the fitted values from an NB GLM be plotted? Defaults to FALSE.
#' @param glm.ci A boolean indicating whether a (\eqn{1 - \alpha})\% CI ribbon should be drawn for the optional \code{glm} model. Defaults to TRUE.
#' @param tradeseq.mod (Optional) An \code{tradeseq} object specifying a GAM fit using the \code{tradeSeq} package. NOTE: this functionality is not fully implemented yet. Defaults to NULL.
#' @param ci.alpha (Optional) The pre-specified Type I Error rate used in generating (\eqn{1 - \alpha})\% CIs. Defaults to good old 0.05.
#' @param plot.alpha (Optional) The opacity of the user-specified confidence interval bands. Defaults to 0.25.
#' @param line.size (Optional) The size of the lines showing the fitted values for each model. Defaults to 1.25.
#' @param plot.theme (Optional) A \code{ggplot2} theme to be added to the plot. Defaults to \code{theme_classic(base_size = 14)}.
#' @return A \code{ggplot} object.
#' @export
#' @examples
#' \dontrun{PlotMARGE(marge_mod, gene.counts = exp_vec, pt = pt_df, gene = "BRCA2", plot.gam = TRUE)}
#' \dontrun{PlotMARGE(marge_mod, gene.counts = exp_vec, pt = pt_df, gene = "BRCA2", marge.ci = TRUE, ci.alpha = 0.1)}

plotModels <- function(marge.preds = NULL,
                       gene.counts = NULL,
                       pt = NULL,
                       gene = NULL,
                       marge.ci = TRUE,
                       plot.breakpoints = FALSE,
                       null.preds = NULL,
                       null.ci = TRUE,
                       plot.gam = FALSE,
                       gam.ci = TRUE,
                       plot.glm = FALSE,
                       glm.ci = TRUE,
                       tradeseq.mod = NULL,
                       ci.alpha = 0.05,
                       plot.alpha = 0.25,
                       line.size = 1.25,
                       plot.theme = NULL) {
  # check inputs
  if (is.null(gene.counts) | is.null(pt) | is.null(gene) | is.null(ci.alpha)) { stop("You forgot one or more of the arguments to plotModels().") }
  # generate parameters for CIs
  Z <- abs(stats::qnorm(ci.alpha / 2))
  if (is.null(plot.theme)) {
    plot.theme <- ggplot2::theme_classic(base_size = 14)
  }
  col_values <- c()
  fill_values <- c()
  # create dataframe of gene / PT data
  plot_df <- data.frame(PT = pt[, 1], counts = gene.counts)
  # generate base plot
  p <- ggplot2::ggplot(plot_df, ggplot2::aes(x = PT, y = counts)) +
       ggplot2::geom_point(alpha = 0.75, size = 2.5) +
       ggplot2::scale_y_continuous(labels = scales::comma_format()) +
       ggplot2::labs(x = "Pseudotime", y = "Expression", title = gene) +
       ggplot2::theme_classic(base_size = 14) +
       ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0.5))

  # optionally add intercept-only model fitted values + CI to plot
  if (!is.null(null.preds)) {
    null_plot_df <- null.preds %>%
                    dplyr::mutate(null_resp = exp(null_link_fit),
                                  null_LL = exp(null_link_fit - Z * null_link_se),
                                  null_UL = exp(null_link_fit + Z * null_link_se),
                                  PT = pt[, 1],
                                  counts = gene.counts)
    p <- p + ggplot2::geom_line(null_plot_df, mapping = ggplot2::aes(x = PT, y = null_resp, color = "Null"), size = line.size)
    col_values <- c(col_values, "forestgreen")
    # optionally add CI for null model (lol)
    if (null.ci) {
      p <- p + ggplot2::geom_ribbon(null_plot_df, mapping = ggplot2::aes(x = PT, ymin = null_LL, ymax = null_UL, fill = "Null"), alpha = plot.alpha) +
               ggplot2::labs(color = "Model", fill = "Model")
      fill_values <- c(fill_values, "forestgreen")
    }
  }

  # optionally add NB GLM fitted values + CI to plot
  if (plot.glm) {
    glm_mod <- MASS::glm.nb(gene.counts ~ ., pt, x = FALSE, y = FALSE, method = "glm.fit2", init.theta = 1)
    glm_plot_df <- data.frame(stats::predict(glm_mod, type = "link", se.fit = TRUE)[1:2]) %>%
                   dplyr::rename(glm_link_fit = fit, glm_link_se = se.fit) %>%
                   dplyr::mutate(glm_resp = exp(glm_link_fit),
                                 glm_LL = exp(glm_link_fit - Z * glm_link_se),
                                 glm_UL = exp(glm_link_fit + Z * glm_link_se),
                                 PT = pt[, 1],
                                 counts = gene.counts)
    p <- p + ggplot2::geom_line(glm_plot_df, mapping = ggplot2::aes(x = PT, y = glm_resp, color = "NB GLM"), size = line.size) +
             ggplot2::labs(color = "Model", fill = "Model")
    col_values <- c(col_values, "goldenrod2")
    # optionally add CI for NB GLM
    if (glm.ci) {
      p <- p + ggplot2::geom_ribbon(glm_plot_df, mapping = ggplot2::aes(x = PT, ymin = glm_LL, ymax = glm_UL, fill = "NB GLM"), alpha = plot.alpha)
      fill_values <- c(fill_values, "goldenrod2")
    }
  }

  # optionally add NB GAM fitted values + CI to plot
  if (plot.gam) {
    gam_mod <- nbGAM(expr = gene.counts, pt = pt)
    gam_plot_df <- data.frame(stats::predict(gam_mod, type = "link", se.fit = TRUE)[1:2]) %>%
                   dplyr::rename(gam_link_fit = fit, gam_link_se = se.fit) %>%
                   dplyr::mutate(gam_resp = exp(gam_link_fit),
                                 gam_LL = exp(gam_link_fit - Z * gam_link_se),
                                 gam_UL = exp(gam_link_fit + Z * gam_link_se),
                                 PT = pt[, 1],
                                 counts = gene.counts)
    p <- p + ggplot2::geom_line(gam_plot_df, mapping = ggplot2::aes(x = PT, y = gam_resp, color = "NB GAM"), size = line.size) +
             ggplot2::labs(color = "Model", fill = "Model")
    col_values <- c(col_values, "firebrick")
    # optionally add CI for NB GAM
    if (gam.ci) {
      p <- p + ggplot2::geom_ribbon(gam_plot_df, mapping = ggplot2::aes(x = PT, ymin = gam_LL, ymax = gam_UL, fill = "NB GAM"), alpha = plot.alpha)
      fill_values <- c(fill_values, "firebrick")
    }
  }

  # optionally add tradeSeq GAM fitted values + CI to plot TODO: Jack
  # if (!is.null(tradeseq.mod)) {
    # tradeSeq_plot_df <- data.frame(stats::predict(tradeseq.mod, type = "link", se.fit = TRUE)[1:2]) %>%
    #                     dplyr::rename(trade_link_fit = fit, trade_link_se = se.fit) %>%
    #                     dplyr::mutate(trade_resp = exp(trade_link_fit),
    #                                   trade_LL = exp(trade_link_fit - Z * trade_link_se),
    #                                   trade_UL = exp(trade_link_fit + Z * trade_link_se),
    #                                   PT = pt[, 1],
    #                                   counts = gene.counts)
  # }

  #  optionally add MARGE model fitted values + CI to plot
  if (!is.null(marge.preds)) {
    marge_plot_df <- marge.preds %>%
                     dplyr::mutate(marge_resp = exp(marge_link_fit),
                                   marge_LL = exp(marge_link_fit - Z * marge_link_se),
                                   marge_UL = exp(marge_link_fit + Z * marge_link_se),
                                   PT = pt[, 1],
                                   counts = gene.counts)
    p <- p + ggplot2::geom_line(marge_plot_df, mapping = ggplot2::aes(x = PT, y = marge_resp, color = "MARGE"), size = line.size) +
             ggplot2::labs(color = "Model", fill = "Model")
    col_values <- c(col_values, "dodgerblue")
    # optionally add CI for MARGE model
    if (marge.ci) {
      p <- p + ggplot2::geom_ribbon(marge_plot_df, mapping = ggplot2::aes(ymin = marge_LL, ymax = marge_UL, fill = "MARGE"), alpha = plot.alpha)
      fill_values <- c(fill_values, "dodgerblue")
    }
  }

  # finalize colors
  p <- p + ggplot2::scale_color_manual(values = col_values) +
           ggplot2::scale_fill_manual(values = fill_values)
  # optionally add vertical lines for MARGE model breakpoints
  if (plot.breakpoints) {
    breakpoints <- extractBreakpoints(model = marge.mod) %>% pull(Breakpoint)
    p <- p + ggplot2::geom_vline(xintercept = breakpoints)
  }
  return(p)
}
