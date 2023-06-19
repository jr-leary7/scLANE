#' Generate a table of fitted values and celltype metadata for genes of interest.
#'
#' @name getFittedValues
#' @author Jack Leary
#' @import magrittr
#' @importFrom stats qnorm
#' @importFrom purrr map map2 reduce
#' @importFrom dplyr mutate bind_cols relocate select filter across
#' @description Generate a table of expression counts, model fitted values, celltype metadata, etc. in order to create custom plots of gene dynamics.
#' @param test.dyn.res The output from \code{\link{testDynamic}}. Defaults to NULL.
#' @param genes A character vector of genes to generate fitted values for. Defaults to NULL.
#' @param pt A data.frame of pseudotime values for each cell. Defaults to NULL.
#' @param expr.mat A matrix of integer expression values for each cell & gene. Must have genes as columns & cells as rows, with column names being gene names. Defaults to NULL.
#' @param size.factor.offset (Optional) An offset to be used to rescale the fitted values. Can be generated easily with \code{\link{createCellOffset}}. No need to provide if the GEE backend was used. Defaults to NULL.
#' @param cell.meta.data (Optional) A data.frame of metadata values for each cell (celltype label, subject characteristics, tissue type, etc.) that will be included in the result table. Defaults to NULL.
#' @param id.vec (Optional) A vector of subject IDs used in fitting GEE or GLMM models. Defaults to NULL.
#' @param ci.alpha (Optional) The pre-specified Type I Error rate used in generating (\eqn{1 - \alpha})\% CIs. Defaults to good old 0.05.
#' @param filter.lineage (Optional) A character vector of lineages to filter out before generating the final plot. Should be letters, i.e. lineage "A" or "B". Defaults to NULL.
#' @return A data.frame containing expression, fitted values, and metadata.
#' @export
#' @examples
#' \dontrun{
#' getFittedValues(gene_stats,
#'                 genes = c("Neurog3", "Epcam", "Krt19"),
#'                 pt = pt_df,
#'                 expr.mat = gene_counts,
#'                 cell.meta.data = seurat_object@meta.data,
#'                 ci.alpha = 0.05)
#' }

getFittedValues <- function(test.dyn.res = NULL,
                            genes = NULL,
                            pt = NULL,
                            expr.mat = NULL,
                            size.factor.offset = NULL,
                            cell.meta.data = NULL,
                            id.vec = NULL,
                            ci.alpha = 0.05,
                            filter.lineage = NULL) {
  # check inputs
  if (is.null(expr.mat) || is.null(pt) || is.null(genes) || is.null(test.dyn.res)) { stop("You forgot one or more of the arguments to getFittedValues().") }
  # generate parameters for CIs
  Z <- stats::qnorm(ci.alpha / 2, lower.tail = FALSE)
  # select sublist for gene of interest
  test.dyn.res <- test.dyn.res[genes]
  # make sure lineages are named nicely
  colnames(pt) <- paste0("Lineage_", LETTERS[1:ncol(pt)])
  lineages <- LETTERS[1:ncol(pt)]
  # create list of lineage-specific dataframes containing expression of each gene
  mod_df_list <- purrr::map2(pt, lineages, \(x, y) {
    mod_df <- purrr::map(genes,
                         \(g) {
                           gene_df <- data.frame(cell = rownames(pt)[!is.na(x)],
                                                 lineage = y,
                                                 pt = x[!is.na(x)],
                                                 gene = g,
                                                 expression = expr.mat[!is.na(x), g])
                           pred_df <- try({
                             data.frame(scLANE_fitted_link = test.dyn.res[[g]][[paste0("Lineage_", y)]]$MARGE_Preds$marge_link_fit,
                                        scLANE_se_link = test.dyn.res[[g]][[paste0("Lineage_", y)]]$MARGE_Preds$marge_link_se)
                           }, silent = TRUE)
                           if (inherits(pred_df, "try-error")) {
                             gene_df <- dplyr::mutate(gene_df,
                                                      scLANE_fitted_link = NA_real_,
                                                      scLANE_se_link = NA_real_)
                           } else {
                             gene_df <- dplyr::mutate(gene_df,
                                                      scLANE_fitted_link = pred_df$scLANE_fitted_link,
                                                      scLANE_se_link = pred_df$scLANE_se_link)
                           }
                           gene_df <- dplyr::mutate(gene_df,
                                                    scLANE_fitted = exp(scLANE_fitted_link),
                                                    scLANE_ci_ll = exp(scLANE_fitted_link - Z * scLANE_se_link),
                                                    scLANE_ci_ul = exp(scLANE_fitted_link + Z * scLANE_se_link))
                           if (!is.null(size.factor.offset)) {
                             gene_df <- dplyr::mutate(gene_df,
                                                      dplyr::across(c(scLANE_fitted, scLANE_ci_ll, scLANE_ci_ul), \(m) m * unname(size.factor.offset)[!is.na(x)]))
                           }
                           gene_df <- dplyr::bind_cols(gene_df, cell.meta.data[!is.na(x), , drop = FALSE])
                           return(gene_df)
                         }) %>%
              purrr::reduce(rbind)
    if (!is.null(id.vec)) {
      mod_df$id <- id.vec[!is.na(x)]
    } else {
      mod_df$id <- NA_character_
    }
    return(mod_df)
  })
  final_df <- purrr::reduce(mod_df_list, rbind) %>%
              dplyr::relocate(cell, id, lineage)
  if (all(is.na(final_df$id))) {
    final_df <- dplyr::select(final_df, -id)
  }
  if (!is.null(filter.lineage)) {
    final_df <- dplyr::filter(final_df, !lineage %in% filter.lineage)
  }
  if (length(genes) > 1) {
    final_df <- dplyr::mutate(final_df, gene = factor(gene, levels = genes))
  }
  return(final_df)
}
