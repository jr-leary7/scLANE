#' Plot results of \code{marge} and other models using \code{ggplot2}.
#'
#' @name plotModels
#' @author Jack Leary
#' @import magrittr
#' @import glm2
#' @importFrom stats qnorm predict as.formula
#' @importFrom purrr map map2 reduce
#' @importFrom dplyr relocate mutate select contains case_when filter if_else rowwise ungroup
#' @importFrom geeM geem
#' @importFrom glmmTMB glmmTMB nbinom2
#' @importFrom MASS negative.binomial theta.mm
#' @importFrom tidyr pivot_longer
#' @importFrom scales label_comma label_number
#' @importFrom ggplot2 ggplot aes geom_point geom_line geom_ribbon facet_wrap scale_y_continuous labs theme element_text guides guide_legend
#' @description This function visualizes the fitted values of several types of models over the expression and pseudotime values of each cell.
#' @param test.dyn.res The output from \code{\link{testDynamic}}. Defaults to NULL.
#' @param gene The name of the gene that's being analyzed. Used as the title of the \code{ggplot} object & to subset the counts matrix. Defaults to NULL.
#' @param pt A data.frame of pseudotime values for each cell. Defaults to NULL.
#' @param expr.mat Either a \code{SingleCellExperiment} or \code{Seurat} object from which counts can be extracted, or a matrix of integer-valued counts. Defaults to NULL.
#' @param size.factor.offset (Optional) An offset to be included in the final model fit. Can be generated easily with \code{\link{createCellOffset}}. Defaults to NULL.
#' @param log1p.norm (Optional) Should log1p-normalized versions of expression & model predictions be returned instead of raw counts? Defaults to TRUE.
#' @param is.gee Should a GEE framework be used instead of the default GLM? Defaults to FALSE.
#' @param is.glmm Should a GLMM framework be used instead of the default GLM? Defaults to FALSE.
#' @param id.vec If the GEE or GLMM framework is being used, a vector of subject IDs to use as input to \code{\link[geeM]{geem}} or \code{\link[glmmTMB]{glmmTMB}}. Defaults to NULL.
#' @param cor.structure If the GEE framework is used, specifies the desired working correlation structure. Must be one of "ar1", "independence", or "exchangeable". Defaults to "ar1".
#' @param ci.alpha (Optional) The pre-specified Type I Error rate used in generating (\eqn{1 - \alpha})\% CIs. Defaults to good old 0.05.
#' @param plot.null (Optional) Should the fitted values from the intercept-only null model be plotted? Defaults to FALSE.
#' @param plot.glm (Optional) Should the fitted values from an NB GLM be plotted? If the data are multi-subject, the "GLM" model can be a GEE or GLMM depending on the desired framework. See Examples for more detail. Defaults to FALSE.
#' @param plot.gam (Optional) Should the fitted values from an NB GAM be plotted? Defaults to FALSE.
#' @param plot.scLANE (Optional) Should the fitted values from the \code{scLANE} model be plotted? Defaults to TRUE.
#' @param filter.lineage (Optional) A character vector of lineages to filter out before generating the final plot. Should be letters, i.e. lineage "A" or "B". Defaults to NULL.
#' @param gg.theme (Optional) A \code{ggplot2} theme to be added to the plot. Defaults to \code{\link{theme_scLANE}}.
#' @return A \code{ggplot} object.
#' @export
#' @examples
#' \dontrun{
#' plotModels(gene_stats,
#'            gene = "AURKA",
#'            pt = pt_df,
#'            expr.mat = count_mat,
#'            size.factor.offset = cell_offset)
#' plotModels(gene_stats,
#'            gene = "CD3E",
#'            pt = pt_df,
#'            expr.mat = seu_obj,
#'            size.factor.offset = cell_offset,
#'            ci.alpha = 0.1,
#'            filter.lineage = c("A", "C"))
#' plotModels(gene_stats,
#'            gene = "CD14",
#'            pt = pt_df,
#'            expr.mat = sce_obj,
#'            size.factor.offset = cell_offset,
#'            is.glmm = TRUE,
#'            id.vec = subject_ids,
#'            plot.glm = TRUE,  # plots an NB GLMM with random intercepts & slopes per-subject
#'            plot.gam = TRUE,  # plots an NB GAMM with random intercepts per-subject
#'            gg.theme = ggplot2::theme_minimal())
#' }

plotModels <- function(test.dyn.res = NULL,
                       gene = NULL,
                       pt = NULL,
                       expr.mat = NULL,
                       size.factor.offset = NULL,
                       log1p.norm = TRUE,
                       is.gee = FALSE,
                       is.glmm = FALSE,
                       id.vec = NULL,
                       cor.structure = "ar1",
                       ci.alpha = 0.05,
                       plot.null = FALSE,
                       plot.glm = FALSE,
                       plot.gam = FALSE,
                       plot.scLANE = TRUE,
                       filter.lineage = NULL,
                       gg.theme = theme_scLANE()) {
  # check inputs
  if (is.null(expr.mat) || is.null(pt) || is.null(gene) || is.null(test.dyn.res)) { stop("You forgot one or more of the arguments to plotModels().") }
  # get raw counts from SingleCellExperiment or Seurat object & transpose to cell x gene dense matrix
  if (inherits(expr.mat, "SingleCellExperiment")) {
    expr.mat <- BiocGenerics::counts(expr.mat)
  } else if (inherits(expr.mat, "Seurat")) {
    expr.mat <- Seurat::GetAssayData(expr.mat,
                                     slot = "counts",
                                     assay = Seurat::DefaultAssay(expr.mat))
  }
  if (!(inherits(expr.mat, "matrix") || inherits(expr.mat, "array") || inherits(expr.mat, "dgCMatrix"))) { stop("Input expr.mat must be coerceable to a matrix of integer counts.") }
  # generate parameters for CIs
  Z <- stats::qnorm(ci.alpha / 2, lower.tail = FALSE)
  # select sublist for gene of interest
  td_res <- test.dyn.res[[gene]]
  # make sure lineages are named nicely
  colnames(pt) <- paste0("Lineage_", LETTERS[seq_len(ncol(pt))])
  # create base list w/ elements being lineage-specific dataframes
  counts_df_list <- purrr::map2(pt,
                                LETTERS[seq_len(ncol(pt))],
                                \(x, y) {
                                  mod_df <- data.frame(CELL = rownames(pt)[!is.na(x)],
                                                       LINEAGE = y,
                                                       PT = x[!is.na(x)],
                                                       COUNT = as.numeric(expr.mat[gene, !is.na(x)]))
                                  if (is.gee || is.glmm) {
                                    mod_df <- dplyr::mutate(mod_df, ID = id.vec[!is.na(x)])
                                  } else {
                                    mod_df <- dplyr::mutate(mod_df, ID = NA_character_)
                                  }
                                  if (!is.null(size.factor.offset)) {
                                    mod_df <- dplyr::mutate(mod_df, CELL_OFFSET = unname(size.factor.offset[!is.na(x)]))
                                  }
                                  return(mod_df)
                                })
  # create list of dataframes w/ predicted values, standard errors, CIs for null, GLM, GAM, & MARGE models
  counts_df_list <- purrr::map(counts_df_list, \(x) dplyr::relocate(x, ID, CELL, LINEAGE, COUNT, PT)) %>%
                    purrr::map2(.y = td_res, \(.x, .y) dplyr::mutate(.x,
                                                                     RESP_MARGE = .y$MARGE_Preds$marge_link_fit,
                                                                     SE_MARGE = .y$MARGE_Preds$marge_link_se,
                                                                     RESP_NULL = .y$Null_Preds$null_link_fit,
                                                                     SE_NULL = .y$Null_Preds$null_link_se)) %>%
                    purrr::map(function(x) {
                      if (plot.glm) {
                        if (is.gee) {
                          theta_hat <- MASS::theta.mm(y = x$COUNT,
                                                      mu = mean(x$COUNT),
                                                      dfr = nrow(x) - 1)
                          gee_formula <- "COUNT ~ PT"
                          if (!is.null(size.factor.offset)) {
                            gee_formula <- paste0(gee_formula, " + offset(log(1 / CELL_OFFSET))")
                          }
                          gee_formula <- stats::as.formula(gee_formula)
                          glm_mod <- geeM::geem(gee_formula,
                                                id = ID,
                                                data = x,
                                                family = MASS::negative.binomial(theta_hat, link = log),
                                                corstr = cor.structure,
                                                scale.fix = FALSE,
                                                sandwich = TRUE)
                          robust_vcov_mat <- as.matrix(glm_mod$var)
                          glm_preds <- data.frame(fit = predict(glm_mod),
                                                  se.fit = sqrt(apply((tcrossprod(glm_mod$X, robust_vcov_mat)) * glm_mod$X, 1, sum)))
                        } else if (is.glmm) {
                          if (is.null(size.factor.offset)) {
                            glm_mod <- glmmTMB::glmmTMB(COUNT ~ PT + (1 + PT | ID),
                                                        data = x,
                                                        family = glmmTMB::nbinom2(link = "log"),
                                                        se = TRUE,
                                                        REML = FALSE)
                          } else {
                            glm_mod <- glmmTMB::glmmTMB(COUNT ~ PT + (1 + PT | ID),
                                                        data = x,
                                                        family = glmmTMB::nbinom2(link = "log"),
                                                        offset = log(1 / x$CELL_OFFSET),
                                                        se = TRUE,
                                                        REML = FALSE)
                          }
                          glm_preds <- data.frame(predict(glm_mod, type = "link", se.fit = TRUE)[seq(2)])
                        } else {
                          glm_formula <- "COUNT ~ PT"
                          if (!is.null(size.factor.offset)) {
                            glm_formula <- paste0(glm_formula, " + offset(log(1 / CELL_OFFSET))")
                          }
                          glm_formula <- stats::as.formula(glm_formula)
                          glm_mod <- MASS::glm.nb(glm_formula,
                                                  data = x,
                                                  x = FALSE,
                                                  y = FALSE,
                                                  method = "glm.fit2",
                                                  link = log,
                                                  init.theta = 1)
                          glm_preds <- data.frame(stats::predict(glm_mod, type = "link", se.fit = TRUE)[seq(2)])
                        }
                        x <- dplyr::mutate(x,
                                           RESP_GLM = glm_preds$fit,
                                           SE_GLM = glm_preds$se.fit)
                      } else {
                        x <- dplyr::mutate(x,
                                           RESP_GLM = NA_real_,
                                           SE_GLM = NA_real_)
                      }
                      return(x)
                    }) %>%
                    purrr::map(function(x) {
                      if (plot.gam) {
                        if (is.null(size.factor.offset)) {
                          if (is.glmm) {
                            gam_mod <- nbGAM(expr = x$COUNT,
                                             pt = x$PT,
                                             id.vec = x$ID)
                          } else {
                            gam_mod <- nbGAM(expr = x$COUNT, pt = x$PT)
                          }
                        } else {
                          if (is.glmm) {
                            gam_mod <- nbGAM(expr = x$COUNT,
                                             Y.offset = x$CELL_OFFSET,
                                             pt = x$PT,
                                             id.vec = x$ID)
                          } else {
                            gam_mod <- nbGAM(expr = x$COUNT,
                                             Y.offset = x$CELL_OFFSET,
                                             pt = x$PT)
                          }
                        }
                        gam_preds <- data.frame(predict(gam_mod, type = "link", se.fit = TRUE)[seq(2)])
                        x <- dplyr::mutate(x,
                                           RESP_GAM = gam_preds$fit,
                                           SE_GAM = gam_preds$se.fit)
                      } else {
                        x <- dplyr::mutate(x,
                                           RESP_GAM = NA_real_,
                                           SE_GAM = NA_real_)
                      }

                      return(x)
                    }) %>%
                    purrr::map(function(x) {
                      dplyr::select(x, CELL, CELL_OFFSET, ID, LINEAGE, COUNT, PT) %>%
                      dplyr::left_join((dplyr::select(x, CELL, ID, LINEAGE, dplyr::contains("RESP")) %>%
                                        tidyr::pivot_longer(cols = dplyr::contains("RESP"),
                                                            values_to = "RESPONSE",
                                                            names_to = "MODEL") %>%
                                        dplyr::mutate(MODEL = gsub("RESP_", "", MODEL))),
                                       by = c("CELL" = "CELL",
                                              "ID" = "ID",
                                              "LINEAGE" = "LINEAGE")) %>%
                      dplyr::left_join((dplyr::select(x, CELL, ID, LINEAGE, dplyr::contains("SE")) %>%
                                        tidyr::pivot_longer(cols = dplyr::contains("SE"),
                                                            values_to = "SE",
                                                            names_to = "MODEL") %>%
                                        dplyr::mutate(MODEL = gsub("SE_", "", MODEL))),
                                       by = c("CELL" = "CELL",
                                              "ID" = "ID",
                                              "LINEAGE" = "LINEAGE",
                                              "MODEL" = "MODEL"))
    })
  counts_df <- purrr::reduce(counts_df_list, rbind) %>%
               dplyr::mutate(MODEL = dplyr::case_when(MODEL == "NULL" ~ "Intercept-only",
                                                      MODEL == "MARGE" ~ "scLANE",
                                                      TRUE ~ MODEL),
                             MODEL = factor(MODEL, levels = c("Intercept-only", "GLM", "GAM", "scLANE")),
                             PRED = exp(RESPONSE),
                             CI_LL = exp(RESPONSE - Z * SE),
                             CI_UL = exp(RESPONSE + Z * SE))
  if (!is.null(size.factor.offset)) {
    if (is.gee) {
      counts_df <- dplyr::rowwise(counts_df) %>%
                   dplyr::mutate(COUNT = COUNT * CELL_OFFSET,
                                 dplyr::across(c(PRED, CI_LL, CI_UL), \(x) {
                                   if (as.character(MODEL) == "GAM") {
                                     x * CELL_OFFSET
                                   } else {
                                     x
                                   }
                                 })) %>%
                   dplyr::ungroup()
    } else {
      counts_df <- dplyr::mutate(counts_df, dplyr::across(c(COUNT, PRED, CI_LL, CI_UL), \(x) x * CELL_OFFSET))
    }
  }
  # add conditional filters here
  if (!plot.null) {
    counts_df <- dplyr::filter(counts_df, MODEL != "Intercept-only")
  }
  if (!plot.glm) {
    counts_df <- dplyr::filter(counts_df, MODEL != "GLM")
  }
  if (!plot.gam) {
    counts_df <- dplyr::filter(counts_df, MODEL != "GAM")
  }
  if (!plot.scLANE) {
    counts_df <- dplyr::filter(counts_df, MODEL != "scLANE")
  }
  if (!is.null(filter.lineage)) {
    counts_df <- dplyr::filter(counts_df, !LINEAGE %in% filter.lineage)
  }
  # change model labels as necessary
  if (is.gee) {
    counts_df <- dplyr::mutate(counts_df,
                               MODEL = dplyr::if_else(as.character(MODEL) == "GLM", "GEE", as.character(MODEL)),
                               MODEL = factor(MODEL, levels = c("Intercept-only", "GEE", "GAM", "scLANE")))
  }
  if (is.glmm) {
    counts_df <- dplyr::mutate(counts_df,
                               MODEL = dplyr::if_else(as.character(MODEL) == "GLM", "GLMM", as.character(MODEL)),
                               MODEL = factor(MODEL, levels = c("Intercept-only", "GLMM", "GAM", "scLANE")))
  }
  if (log1p.norm) {
    counts_df <- dplyr::mutate(counts_df, dplyr::across(c(COUNT, PRED, CI_LL, CI_UL), log1p))
  }
  # generate plot
  if (is.glmm) {
    p <- ggplot2::ggplot(counts_df, ggplot2::aes(x = PT, y = COUNT, group = ID)) +
         ggplot2::geom_point(mapping = ggplot2::aes(color = ID),
                             alpha = 0.5,
                             size = 0.5)
    if (requireNamespace("ggh4x", quietly = TRUE)) {
      p <- p + ggh4x::facet_nested_wrap(~paste0("Lineage ", LINEAGE) + MODEL + ID,
                                        nrow = length(levels(counts_df$MODEL)),
                                        strip = ggh4x::strip_nested(clip = "off",
                                                                    background_x = list(ggplot2::element_rect(linewidth = gg.theme$line$linewidth))))
    } else {
      p <- p + ggplot2::facet_wrap(~paste0("Lineage ", LINEAGE) + MODEL + ID)
    }
    p <- p +
         ggplot2::geom_line(mapping = ggplot2::aes(y = PRED),
                            linewidth = 0.75,
                            color = "black",
                            show.legend = FALSE) +
         ggplot2::geom_ribbon(mapping = ggplot2::aes(ymin = CI_LL, ymax = CI_UL),
                              fill = "grey50",
                              alpha = 0.4,
                              linewidth = 0,
                              show.legend = FALSE) +
         ggplot2::scale_y_continuous(labels = scales::label_comma()) +
         ggplot2::scale_x_continuous(labels = scales::label_number(accuracy = 0.1)) +
         ggplot2::labs(x = "Pseudotime",
                       y = "Expression",
                       color = "Subject",
                       fill = "Subject",
                       title = gene) +
         gg.theme +
         ggplot2::theme(plot.title = ggplot2::element_text(face = "italic"),
                        strip.clip = "off",
                        strip.background = ggplot2::element_rect(linewidth = gg.theme$line$linewidth)) +
         ggplot2::guides(color = ggplot2::guide_legend(override.aes = list(size = 4, alpha = 1)))
  } else {
    p <- ggplot2::ggplot(counts_df, ggplot2::aes(x = PT, y = COUNT, color = LINEAGE)) +
         ggplot2::geom_point(alpha = 0.5,
                             size = 0.5,
                             show.legend = ifelse(ncol(pt) > 1, TRUE, FALSE))
    if (requireNamespace("ggh4x", quietly = TRUE)) {
      p <- p + ggh4x::facet_nested_wrap(~paste0("Lineage ", LINEAGE) + MODEL,
                                        strip = ggh4x::strip_nested(clip = "off",
                                                                    background_x = list(ggplot2::element_rect(linewidth = gg.theme$line$linewidth))))
    } else {
      p <- p + ggplot2::facet_wrap(~paste0("Lineage ", LINEAGE) + MODEL)
    }
    p <- p +
         ggplot2::geom_line(mapping = ggplot2::aes(y = PRED),
                            linewidth = 0.75,
                            color = "black") +
         ggplot2::geom_ribbon(mapping = ggplot2::aes(ymin = CI_LL, ymax = CI_UL),
                              alpha = 0.4,
                              linewidth = 0,
                              color = "grey") +
         ggplot2::scale_y_continuous(labels = scales::label_comma()) +
         ggplot2::scale_x_continuous(labels = scales::label_number(accuracy = 0.1)) +
         ggplot2::labs(x = "Pseudotime",
                       y = "Expression",
                       color = "Lineage",
                       fill = "Lineage",
                       title = gene) +
         gg.theme +
         ggplot2::theme(plot.title = ggplot2::element_text(face = "italic"),
                        strip.clip = "off",
                        strip.background = ggplot2::element_rect(linewidth = gg.theme$line$linewidth)) +
         ggplot2::guides(color = ggplot2::guide_legend(override.aes = list(size = 4, alpha = 1)))
  }
  return(p)
}
