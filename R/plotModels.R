#' Plot results of \code{marge} and other models using \code{ggplot2}.
#'
#' @name plotModels
#' @author Jack Leary
#' @import magrittr
#' @importFrom ggplot2 theme_classic ggplot aes geom_point geom_line geom_ribbon facet_wrap scale_y_continuous labs theme element_text
#' @importFrom stats qnorm predict
#' @importFrom purrr map map2 reduce
#' @importFrom dplyr relocate mutate select contains case_when filter
#' @importFrom tidyr pivot_longer
#' @importFrom scales comma_format
#' @description This function visualizes the fitted values of several types of models over the expression and pseudotime values of each cell.
#' @param test.dyn.res The output from \code{\link{testDynamic}}. Defaults to NULL.
#' @param gene The name of the gene that's being analyzed. Used as the title of the \code{ggplot} object & to subset the counts matrix. Defaults to NULL.
#' @param pt A data.frame of pseudotime values for each cell. Defaults to NULL.
#' @param gene.counts A matrix of integer expression values for each cell & gene. Must have genes as columns & cells as rows, with column names being gene names. Defaults to NULL.
#' @param ci.alpha (Optional) The pre-specified Type I Error rate used in generating (\eqn{1 - \alpha})\% CIs. Defaults to good old 0.05.
#' @param plot.null (Optional) Should the fitted values from the intercept-only null model be plotted? Defaults to TRUE.
#' @param plot.glm (Optional) Should the fitted values from an NB GLM be plotted? Defaults to TRUE.
#' @param plot.gam (Optional) Should the fitted values from an NB GAM be plotted? Defaults to TRUE.
#' @param plot.marge (Optional) Should the fitted values from the \code{marge} model be plotted? Defaults to TRUE.
#' @param filter.lineage (Optional) A character vector of lineages to filter out before generating the final plot. Should be letters, i.e. lineage "A" or "B". Defaults to NULL.
#' @param gg.theme (Optional) A \code{ggplot2} theme to be added to the plot. Defaults to \code{theme_classic(base_size = 14)}.
#' @return A \code{ggplot} object.
#' @export
#' @examples
#' \dontrun{PlotMARGE(gene_stats, gene = "AURKA", pt = pt_df, gene.counts = count_mat)}
#' \dontrun{PlotMARGE(gene_stats, gene = "BRCA2", pt = pt_df, gene.counts = count_mat, ci.alpha = 0.1, filter.lineage = c("A", "C"))}

plotModels <- function(test.dyn.res = NULL,
                       gene = NULL,
                       pt = NULL,
                       gene.counts = NULL,
                       ci.alpha = 0.05,
                       plot.null = TRUE,
                       plot.glm = TRUE,
                       plot.gam = TRUE,
                       plot.marge = TRUE,
                       filter.lineage = NULL,
                       gg.theme = ggplot2::theme_classic(base_size = 14)) {
  # check inputs
  if (is.null(gene.counts) | is.null(pt) | is.null(gene) | is.null(test.dyn.res)) { stop("You forgot one or more of the arguments to plotModels().") }
  # generate parameters for CIs
  Z <- abs(stats::qnorm(ci.alpha / 2))
  # select sublist for gene of interest
  td_res <- test.dyn.res[[gene]]
  # make sure lineages are named nicely
  colnames(pt) <- paste0("Lineage_", LETTERS[1:ncol(pt)])
  # create base list w/ elements being lineage-specific dataframes
  counts_df_list <- purrr::map(pt, function(x) data.frame(CELL = rownames(pt)[!is.na(x)],
                                                          PT = x[!is.na(x)],
                                                          COUNT = gene.counts[!is.na(x), gene]))
  for (i in seq_along(counts_df_list)) { counts_df_list[[i]]$LINEAGE <- LETTERS[i] }
  # create list of dataframes w/ predicted values, standard errors, CIs for null, GLM, GAM, & MARGE models
  counts_df_list <- purrr::map(counts_df_list, function(x) x %>% dplyr::relocate(CELL, LINEAGE, COUNT, PT)) %>%
                    purrr::map2(.y = td_res, function(.x, .y) .x %>% dplyr::mutate(RESP_MARGE = .y$MARGE_Preds$marge_link_fit,
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
                      glm_mod <- MASS::glm.nb(x$COUNT ~ x$PT, x = FALSE, y = FALSE, method = "glm.fit2", init.theta = 1)
                      glm_preds <- data.frame(stats::predict(glm_mod, type = "link", se.fit = TRUE)[1:2])
                      x %<>% dplyr::mutate(RESP_GLM = glm_preds$fit,
                                           SE_GLM = glm_preds$se.fit,
                                           PRED_GLM = exp(RESP_GLM),
                                           CI_LL_GLM = exp(RESP_GLM - Z * SE_GLM),
                                           CI_UL_GLM = exp(RESP_GLM + Z * SE_GLM))
                      return(x)
                    }) %>%
                    purrr::map(function(x) {
                      gam_mod <- nbGAM(expr = x$COUNT, pt = x$PT, theta.init = 1)
                      gam_preds <- data.frame(stats::predict(gam_mod, type = "link", se.fit = TRUE)[1:2])
                      x %<>% dplyr::mutate(RESP_GAM = gam_preds$fit,
                                           SE_GAM = gam_preds$se.fit,
                                           PRED_GAM = exp(RESP_GAM),
                                           CI_LL_GAM = exp(RESP_GAM - Z * SE_GAM),
                                           CI_UL_GAM = exp(RESP_GAM + Z * SE_GAM))
                      return(x)
                    }) %>%
                    purrr::map(function(x) {
                      dplyr::select(x, CELL, LINEAGE, COUNT, PT) %>%
                      dplyr::left_join((x %>%
                                        dplyr::select(CELL, LINEAGE, dplyr::contains("RESP")) %>%
                                        tidyr::pivot_longer(cols = dplyr::contains("RESP"), values_to = "RESPONSE", names_to = "MODEL") %>%
                                        dplyr::mutate(MODEL = gsub("RESP_", "", MODEL))),
                                       by = c("CELL" = "CELL", "LINEAGE" = "LINEAGE")) %>%
                      dplyr::left_join((x %>%
                                        dplyr::select(CELL, LINEAGE, dplyr::contains("SE")) %>%
                                        tidyr::pivot_longer(cols = dplyr::contains("SE"), values_to = "SE", names_to = "MODEL") %>%
                                        dplyr::mutate(MODEL = gsub("SE_", "", MODEL))),
                                       by = c("CELL" = "CELL", "LINEAGE" = "LINEAGE", "MODEL" = "MODEL")) %>%
                      dplyr::left_join((x %>%
                                        dplyr::select(CELL, LINEAGE, dplyr::contains("PRED")) %>%
                                        tidyr::pivot_longer(cols = dplyr::contains("PRED"), values_to = "PRED", names_to = "MODEL") %>%
                                        dplyr::mutate(MODEL = gsub("PRED_", "", MODEL))),
                                       by = c("CELL" = "CELL", "LINEAGE" = "LINEAGE", "MODEL" = "MODEL")) %>%
                      dplyr::left_join((x %>%
                                        dplyr::select(CELL, LINEAGE, dplyr::contains("CI_LL")) %>%
                                        tidyr::pivot_longer(cols = dplyr::contains("CI_LL"), values_to = "CI_LL", names_to = "MODEL") %>%
                                        dplyr::mutate(MODEL = gsub("CI_LL_", "", MODEL))),
                                       by = c("CELL" = "CELL", "LINEAGE" = "LINEAGE", "MODEL" = "MODEL")) %>%
                      dplyr::left_join((x %>%
                                        dplyr::select(CELL, LINEAGE, dplyr::contains("CI_UL")) %>%
                                        tidyr::pivot_longer(cols = dplyr::contains("CI_UL"), values_to = "CI_UL", names_to = "MODEL") %>%
                                        dplyr::mutate(MODEL = gsub("CI_UL_", "", MODEL))),
                                       by = c("CELL" = "CELL", "LINEAGE" = "LINEAGE", "MODEL" = "MODEL"))
    })
  counts_df <- counts_df_list %>%
               purrr::reduce(rbind) %>%
               dplyr::mutate(MODEL = factor(dplyr::case_when(MODEL == "NULL" ~ "Intercept-only", TRUE ~ MODEL),
                                           levels = c("Intercept-only", "GLM", "GAM", "MARGE")))
  # add conditional filters here
  if (!plot.null) { counts_df %<>% dplyr::filter(MODEL != "Intercept-only") }
  if (!plot.glm) { counts_df %<>% dplyr::filter(MODEL != "GLM") }
  if (!plot.gam) { counts_df %<>% dplyr::filter(MODEL != "GAM") }
  if (!plot.marge) { counts_df %<>% dplyr::filter(MODEL != "MARGE") }
  if (!is.null(filter.lineage)) { counts_df %<>% dplyr::filter(!LINEAGE %in% filter.lineage) }
  # generate plot
  p <- ggplot2::ggplot(counts_df, ggplot2::aes(x = PT, y = COUNT, color = LINEAGE)) +
       ggplot2::geom_point(alpha = 0.5, size = 0.5) +
       ggplot2::geom_line(ggplot2::aes(x = PT, y = PRED, color = LINEAGE), size = 1) +
       ggplot2::geom_ribbon(mapping = ggplot2::aes(x = PT, ymin = CI_LL, ymax = CI_UL, fill = LINEAGE), alpha = 0.5, size = 0) +
       ggplot2::facet_wrap(~MODEL) +
       ggplot2::scale_y_continuous(labels = scales::comma_format()) +
       ggplot2::labs(x = "Pseudotime", y = "Expression", color = "Lineage", fill = "Lineage", title = gene) +
       gg.theme +
       ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0.5))
  return(p)
}
