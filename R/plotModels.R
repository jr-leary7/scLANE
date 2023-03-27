#' Plot results of \code{marge} and other models using \code{ggplot2}.
#'
#' @name plotModels
#' @author Jack Leary
#' @import magrittr
#' @import glm2
#' @importFrom stats qnorm predict
#' @importFrom purrr map map2 reduce
#' @importFrom dplyr relocate mutate select contains case_when filter if_else
#' @importFrom geeM geem
#' @importFrom glmmTMB glmmTMB nbinom2
#' @importFrom MASS negative.binomial theta.mm
#' @importFrom tidyr pivot_longer
#' @importFrom scales label_comma label_number
#' @importFrom ggplot2 theme_classic ggplot aes geom_point geom_line geom_ribbon facet_wrap scale_y_continuous labs theme element_text guides guide_legend
#' @description This function visualizes the fitted values of several types of models over the expression and pseudotime values of each cell.
#' @param test.dyn.res The output from \code{\link{testDynamic}}. Defaults to NULL.
#' @param gene The name of the gene that's being analyzed. Used as the title of the \code{ggplot} object & to subset the counts matrix. Defaults to NULL.
#' @param pt A data.frame of pseudotime values for each cell. Defaults to NULL.
#' @param gene.counts A matrix of integer expression values for each cell & gene. Must have genes as columns & cells as rows, with column names being gene names. Defaults to NULL.
#' @param is.gee Should a GEE framework be used instead of the default GLM? Defaults to FALSE.
#' @param is.glmm Should a GLMM framework be used instead of the default GLM? Defaults to FALSE.
#' @param id.vec If the GEE or GLMM framework is being used, a vector of subject IDs to use as input to \code{\link[geeM]{geem}} or \code{\link[glmmTMB]{glmmTMB}}. Defaults to NULL.
#' @param cor.structure If the GEE framework is used, specifies the desired working correlation structure. Must be one of "ar1", "independence", "unstructured", or "exchangeable". Defaults to "exchangeable".
#' @param ci.alpha (Optional) The pre-specified Type I Error rate used in generating (\eqn{1 - \alpha})\% CIs. Defaults to good old 0.05.
#' @param plot.null (Optional) Should the fitted values from the intercept-only null model be plotted? Defaults to TRUE.
#' @param plot.glm (Optional) Should the fitted values from an NB GLM be plotted? If the data are multi-subject, the "GLM" model can be a GEE or GLMM depending on the desired framework. See Examples for more detail. Defaults to TRUE.
#' @param plot.gam (Optional) Should the fitted values from an NB GAM be plotted? Defaults to TRUE.
#' @param plot.marge (Optional) Should the fitted values from the \code{marge} model be plotted? Defaults to TRUE.
#' @param filter.lineage (Optional) A character vector of lineages to filter out before generating the final plot. Should be letters, i.e. lineage "A" or "B". Defaults to NULL.
#' @param gg.theme (Optional) A \code{ggplot2} theme to be added to the plot. Defaults to \code{theme_classic(base_size = 14)}.
#' @return A \code{ggplot} object.
#' @export
#' @examples
#' \dontrun{
#' plotModels(gene_stats,
#'            gene = "AURKA",
#'            pt = pt_df,
#'            gene.counts = count_mat)
#' plotModels(gene_stats,
#'            gene = "BRCA2",
#'            pt = pt_df,
#'            gene.counts = count_mat,
#'            ci.alpha = 0.1,
#'            filter.lineage = c("A", "C"))
#' plotModels(gene_stats,
#'            gene = "CD14",
#'            pt = pt_df,
#'            is.glmm = TRUE,
#'            id.vec = subject_ids,
#'            plot.glm = TRUE,  # plots an NB GLMM with random intercepts & slopes per-subject
#'            plot.gam = TRUE,  # plots an NB GAMM with random intercepts per-subject
#'            gg.theme = ggplot2::theme_minimal())
#' }

plotModels <- function(test.dyn.res = NULL,
                       gene = NULL,
                       pt = NULL,
                       gene.counts = NULL,
                       is.gee = FALSE,
                       is.glmm = FALSE,
                       id.vec = NULL,
                       cor.structure = "exchangeable",
                       ci.alpha = 0.05,
                       plot.null = TRUE,
                       plot.glm = TRUE,
                       plot.gam = TRUE,
                       plot.marge = TRUE,
                       filter.lineage = NULL,
                       gg.theme = ggplot2::theme_classic(base_size = 14)) {
  # check inputs
  if (is.null(gene.counts) || is.null(pt) || is.null(gene) || is.null(test.dyn.res)) { stop("You forgot one or more of the arguments to plotModels().") }
  # generate parameters for CIs
  Z <- stats::qnorm(ci.alpha / 2, lower.tail = FALSE)
  # select sublist for gene of interest
  td_res <- test.dyn.res[[gene]]
  # make sure lineages are named nicely
  colnames(pt) <- paste0("Lineage_", LETTERS[1:ncol(pt)])
  # create base list w/ elements being lineage-specific dataframes
  counts_df_list <- purrr::map2(pt,
                                LETTERS[1:ncol(pt)],
                                \(x, y) {
                                  mod_df <- data.frame(CELL = rownames(pt)[!is.na(x)],
                                                       LINEAGE = y,
                                                       PT = x[!is.na(x)],
                                                       COUNT = gene.counts[!is.na(x), gene])
                                  if (is.gee || is.glmm) {
                                    mod_df$ID <- id.vec[!is.na(x)]
                                  } else {
                                    mod_df$ID <- NA_character_
                                  }
                                  return(mod_df)
                                })
  # create list of dataframes w/ predicted values, standard errors, CIs for null, GLM, GAM, & MARGE models
  counts_df_list <- purrr::map(counts_df_list, \(x) dplyr::relocate(x, ID, CELL, LINEAGE, COUNT, PT)) %>%
                    purrr::map2(.y = td_res, \(.x, .y) dplyr::mutate(.x,
                                                                     RESP_MARGE = .y$MARGE_Preds$marge_link_fit,
                                                                     SE_MARGE = .y$MARGE_Preds$marge_link_se,
                                                                     PRED_MARGE = exp(RESP_MARGE),
                                                                     CI_LL_MARGE = exp(RESP_MARGE - Z * SE_MARGE),
                                                                     CI_UL_MARGE = exp(RESP_MARGE + Z * SE_MARGE),
                                                                     RESP_NULL = .y$Null_Preds$null_link_fit,
                                                                     SE_NULL = .y$Null_Preds$null_link_se,
                                                                     PRED_NULL = exp(RESP_NULL),
                                                                     CI_LL_NULL = exp(RESP_NULL - Z * SE_NULL),
                                                                     CI_UL_NULL = exp(RESP_NULL + Z * SE_NULL))) %>%
                    purrr::map(function(x) {
                      if (is.gee) {
                        theta_hat <- MASS::theta.mm(y = x$COUNT,
                                                    mu = mean(x$COUNT),
                                                    dfr = nrow(x) - 1)
                        glm_mod <- geeM::geem(x$COUNT ~ x$PT,
                                              id = x$ID,
                                              family = MASS::negative.binomial(theta_hat),
                                              corstr = cor.structure,
                                              sandwich = TRUE)
                        robust_vcov_mat <- as.matrix(glm_mod$var)
                        glm_preds <- data.frame(fit = predict(glm_mod),
                                                se.fit = sqrt(apply((tcrossprod(glm_mod$X, robust_vcov_mat)) * glm_mod$X, 1, sum)))
                      } else if (is.glmm) {
                        mod_df <- data.frame(COUNT = x$COUNT,
                                             PT = x$PT,
                                             ID = x$ID)
                        glm_mod <- glmmTMB::glmmTMB(COUNT ~ PT + (1 | ID) + (1 + PT | ID),
                                                    data = mod_df,
                                                    family = glmmTMB::nbinom2(link = "log"),
                                                    se = TRUE,
                                                    REML = FALSE)
                        glm_preds <- data.frame(predict(glm_mod, type = "link", se.fit = TRUE)[1:2])
                      } else {
                        glm_mod <- MASS::glm.nb(x$COUNT ~ x$PT,
                                                x = FALSE,
                                                y = FALSE,
                                                method = "glm.fit2",
                                                init.theta = 1)
                        glm_preds <- data.frame(stats::predict(glm_mod, type = "link", se.fit = TRUE)[1:2])
                      }
                      x %<>% dplyr::mutate(RESP_GLM = glm_preds$fit,
                                           SE_GLM = glm_preds$se.fit,
                                           PRED_GLM = exp(RESP_GLM),
                                           CI_LL_GLM = exp(RESP_GLM - Z * SE_GLM),
                                           CI_UL_GLM = exp(RESP_GLM + Z * SE_GLM))
                      return(x)
                    }) %>%
                    purrr::map(function(x) {
                      if (is.glmm) {
                        gam_mod <- nbGAM(expr = x$COUNT,
                                         pt = x$PT,
                                         id.vec = x$ID)
                      } else {
                        gam_mod <- nbGAM(expr = x$COUNT, pt = x$PT)
                      }
                      gam_preds <- data.frame(predict(gam_mod, type = "link", se.fit = TRUE)[1:2])
                      x %<>% dplyr::mutate(RESP_GAM = gam_preds$fit,
                                           SE_GAM = gam_preds$se.fit,
                                           PRED_GAM = exp(RESP_GAM),
                                           CI_LL_GAM = exp(RESP_GAM - Z * SE_GAM),
                                           CI_UL_GAM = exp(RESP_GAM + Z * SE_GAM))
                      return(x)
                    }) %>%
                    purrr::map(function(x) {
                      dplyr::select(x, CELL, ID, LINEAGE, COUNT, PT) %>%
                      dplyr::left_join((dplyr::select(x, CELL, ID, LINEAGE, dplyr::contains("RESP")) %>%
                                        tidyr::pivot_longer(cols = dplyr::contains("RESP"), values_to = "RESPONSE", names_to = "MODEL") %>%
                                        dplyr::mutate(MODEL = gsub("RESP_", "", MODEL))),
                                       by = c("CELL" = "CELL", "ID" = "ID", "LINEAGE" = "LINEAGE")) %>%
                      dplyr::left_join((dplyr::select(x, CELL, ID, LINEAGE, dplyr::contains("SE")) %>%
                                        tidyr::pivot_longer(cols = dplyr::contains("SE"), values_to = "SE", names_to = "MODEL") %>%
                                        dplyr::mutate(MODEL = gsub("SE_", "", MODEL))),
                                       by = c("CELL" = "CELL", "ID" = "ID", "LINEAGE" = "LINEAGE", "MODEL" = "MODEL")) %>%
                      dplyr::left_join((dplyr::select(x, CELL, ID, LINEAGE, dplyr::contains("PRED")) %>%
                                        tidyr::pivot_longer(cols = dplyr::contains("PRED"), values_to = "PRED", names_to = "MODEL") %>%
                                        dplyr::mutate(MODEL = gsub("PRED_", "", MODEL))),
                                       by = c("CELL" = "CELL", "ID" = "ID", "LINEAGE" = "LINEAGE", "MODEL" = "MODEL")) %>%
                      dplyr::left_join((dplyr::select(x, CELL, ID, LINEAGE, dplyr::contains("CI_LL")) %>%
                                        tidyr::pivot_longer(cols = dplyr::contains("CI_LL"), values_to = "CI_LL", names_to = "MODEL") %>%
                                        dplyr::mutate(MODEL = gsub("CI_LL_", "", MODEL))),
                                       by = c("CELL" = "CELL", "ID" = "ID", "LINEAGE" = "LINEAGE", "MODEL" = "MODEL")) %>%
                      dplyr::left_join((dplyr::select(x, CELL, ID, LINEAGE, dplyr::contains("CI_UL")) %>%
                                        tidyr::pivot_longer(cols = dplyr::contains("CI_UL"), values_to = "CI_UL", names_to = "MODEL") %>%
                                        dplyr::mutate(MODEL = gsub("CI_UL_", "", MODEL))),
                                       by = c("CELL" = "CELL", "ID" = "ID", "LINEAGE" = "LINEAGE", "MODEL" = "MODEL"))
    })
  counts_df <- purrr::reduce(counts_df_list, rbind) %>%
               dplyr::mutate(MODEL = dplyr::case_when(MODEL == "NULL" ~ "Intercept-only",
                                                      MODEL == "MARGE" ~ "scLANE",
                                                      TRUE ~ MODEL))
  # add conditional filters here
  if (!plot.null) { counts_df %<>% dplyr::filter(MODEL != "Intercept-only") }
  if (!plot.glm) { counts_df %<>% dplyr::filter(MODEL != "GLM") }
  if (!plot.gam) { counts_df %<>% dplyr::filter(MODEL != "GAM") }
  if (!plot.marge) { counts_df %<>% dplyr::filter(MODEL != "scLANE") }
  if (!is.null(filter.lineage)) { counts_df %<>% dplyr::filter(!LINEAGE %in% filter.lineage) }
  # change model labels as necessary
  if (is.gee) {
    counts_df %<>% dplyr::mutate(MODEL = dplyr::if_else(MODEL == "GLM", "GEE", MODEL),
                                 MODEL = factor(MODEL, levels = c("Intercept-only", "GEE", "GAM", "scLANE")))
  }
  if (is.glmm) {
    counts_df %<>% dplyr::mutate(MODEL = dplyr::if_else(MODEL == "GLM", "GLMM", MODEL),
                                 MODEL = factor(MODEL, levels = c("Intercept-only", "GLMM", "GAM", "scLANE")))
  }
  # generate plot
  if (is.glmm) {
    p <- ggplot2::ggplot(counts_df, ggplot2::aes(x = PT, y = COUNT, color = ID, fill = ID)) +
         ggplot2::geom_point(alpha = 0.4, size = 0.5, show.legend = FALSE) +
         ggplot2::facet_wrap(~LINEAGE + MODEL) +
         ggplot2::geom_line(mapping = ggplot2::aes(x = PT, y = PRED, group = ID), linewidth = 1) +
         ggplot2::geom_ribbon(mapping = ggplot2::aes(x = PT, ymin = CI_LL, ymax = CI_UL),
                              alpha = 0.4,
                              size = 0,
                              show.legend = FALSE) +
         ggplot2::scale_y_continuous(labels = scales::label_comma()) +
         ggplot2::scale_x_continuous(labels = scales::label_number(accuracy = 0.1)) +
         ggplot2::labs(x = "Pseudotime",
                       y = "Expression",
                       color = "Subject",
                       fill = "Subject",
                       title = gene) +
         gg.theme +
         ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0.5)) +
         ggplot2::guides(color = ggplot2::guide_legend(override.aes = list(size = 2, alpha = 1)))
  } else {
    p <- ggplot2::ggplot(counts_df, ggplot2::aes(x = PT, y = COUNT, color = LINEAGE)) +
         ggplot2::geom_point(alpha = 0.4, size = 0.5) +
         ggplot2::facet_wrap(~LINEAGE + MODEL) +
         ggplot2::geom_line(mapping = ggplot2::aes(x = PT, y = PRED), linewidth = 1, color = "black") +
         ggplot2::geom_ribbon(mapping = ggplot2::aes(x = PT, ymin = CI_LL, ymax = CI_UL),
                              alpha = 0.4,
                              size = 0,
                              color = "grey") +
         ggplot2::scale_y_continuous(labels = scales::label_comma()) +
         ggplot2::scale_x_continuous(labels = scales::label_number(accuracy = 0.1)) +
         ggplot2::labs(x = "Pseudotime",
                       y = "Expression",
                       color = "Lineage",
                       fill = "Lineage",
                       title = gene) +
         gg.theme +
         ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0.5)) +
         ggplot2::guides(color = ggplot2::guide_legend(override.aes = list(size = 2, alpha = 1)))
  }
  return(p)
}
