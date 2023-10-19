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
#' @param expr.mat Either a \code{SingleCellExperiment} or \code{Seurat} object from which counts can be extracted, or a matrix of integer-valued counts with genes as rows & cells as columns. Defaults to NULL.
#' @param size.factor.offset (Optional) An offset to be used to rescale the fitted values. Can be generated easily with \code{\link{createCellOffset}}. No need to provide if the GEE backend was used. Defaults to NULL.
#' @param log1p.norm (Optional) Should log1p-normalized versions of expression & model predictions be returned as well? Defaults to TRUE.
#' @param cell.meta.data (Optional) A data.frame of metadata values for each cell (celltype label, subject characteristics, tissue type, etc.) that will be included in the result table. Defaults to NULL.
#' @param id.vec (Optional) A vector of subject IDs used in fitting GEE or GLMM models. Defaults to NULL.
#' @param ci.alpha (Optional) The pre-specified Type I Error rate used in generating (\eqn{1 - \alpha})\% CIs. Defaults to good old 0.05.
#' @param filter.lineage (Optional) A character vector of lineages to filter out before generating the final plot. Should be letters, i.e. lineage "A" or "B". Defaults to NULL.
#' @return A data.frame containing depth- and log1p-normalized expression, model predictions, and cell-level metadata.
#' @export
#' @examples
#' \dontrun{
#' getFittedValues(gene_stats,
#'                 genes = c("Neurog3", "Epcam", "Krt19"),
#'                 pt = pt_df,
#'                 expr.mat = seu_obj,
#'                 cell.meta.data = seu_obj@meta.data,
#'                 ci.alpha = 0.05)
#' }

getFittedValues <- function(test.dyn.res = NULL,
                            genes = NULL,
                            pt = NULL,
                            expr.mat = NULL,
                            size.factor.offset = NULL,
                            log1p.norm = TRUE,
                            cell.meta.data = NULL,
                            id.vec = NULL,
                            ci.alpha = 0.05,
                            filter.lineage = NULL) {
  # check inputs
  if (is.null(expr.mat) || is.null(pt) || is.null(genes) || is.null(test.dyn.res)) { stop("You forgot one or more of the arguments to getFittedValues().") }
  # get raw counts from SingleCellExperiment or Seurat object & transpose to cell x gene dense matrix
  if (inherits(expr.mat, "SingleCellExperiment")) {
    expr.mat <- BiocGenerics::counts(expr.mat)[genes, , drop = FALSE]
    expr.mat <- as.matrix(expr.mat)
  } else if (inherits(expr.mat, "Seurat")) {
    expr.mat <- Seurat::GetAssayData(expr.mat,
                                     slot = "counts",
                                     assay = Seurat::DefaultAssay(expr.mat))
    expr.mat <- as.matrix(expr.mat[genes, , drop = FALSE])
  } else if (inherits(expr.mat, "dgCMatrix") || inherits(expr.mat, "dgTMatrix")) {
    expr.mat <- as.matrix(expr.mat[genes, , drop = FALSE])
  }
  if (!(inherits(expr.mat, "matrix") || inherits(expr.mat, "array"))) { stop("Input expr.mat must be coerceable to a matrix of integer counts.") }
  expr.mat <- t(expr.mat)  # transpose to cell x gene matrix
  if (is.null(genes)) {
    genes <- colnames(expr.mat)
  }
  # generate parameters for CIs
  Z <- stats::qnorm(ci.alpha / 2, lower.tail = FALSE)
  # select sublist for gene of interest
  test.dyn.res <- test.dyn.res[genes]
  # make sure lineages are named nicely
  colnames(pt) <- paste0("Lineage_", LETTERS[seq_len(ncol(pt))])
  lineages <- LETTERS[seq_len(ncol(pt))]
  # create list of lineage-specific dataframes containing expression of each gene
  mod_df_list <- purrr::map2(pt, lineages, \(x, y) {
    mod_df <- purrr::map(genes,
                         \(g) {
                           gene_df <- data.frame(cell = rownames(expr.mat)[!is.na(x)],
                                                 lineage = y,
                                                 pt = x[!is.na(x)],
                                                 gene = g,
                                                 rna = expr.mat[!is.na(x), g])
                           if (!is.null(size.factor.offset)) {
                             gene_df <- dplyr::mutate(gene_df,
                                                      size_factor = unname(size.factor.offset)[!is.na(x)],
                                                      model_offset = log(1 / unname(size.factor.offset)[!is.na(x)]))
                           }
                           pred_df <- try({
                             data.frame(scLANE_pred_link = test.dyn.res[[g]][[paste0("Lineage_", y)]]$MARGE_Preds$marge_link_fit,
                                        scLANE_se_link = test.dyn.res[[g]][[paste0("Lineage_", y)]]$MARGE_Preds$marge_link_se)
                           }, silent = TRUE)
                           if (inherits(pred_df, "try-error") || is.null(pred_df) || all(is.na(pred_df))) {
                             gene_df <- dplyr::mutate(gene_df,
                                                      scLANE_pred_link = NA_real_,
                                                      scLANE_se_link = NA_real_,
                                                      scLANE_ci_ll_link = NA_real_,
                                                      scLANE_ci_ul_link = NA_real_)
                           } else {
                             gene_df <- dplyr::mutate(gene_df,
                                                      scLANE_pred_link = pred_df$scLANE_pred_link,
                                                      scLANE_se_link = pred_df$scLANE_se_link,
                                                      scLANE_ci_ll_link = scLANE_pred_link - Z * scLANE_se_link,
                                                      scLANE_ci_ul_link = scLANE_pred_link + Z * scLANE_se_link)
                           }
                           gene_df <- dplyr::mutate(gene_df,
                                                    scLANE_pred = exp(scLANE_pred_link),
                                                    scLANE_ci_ll = exp(scLANE_pred_link - Z * scLANE_se_link),
                                                    scLANE_ci_ul = exp(scLANE_pred_link + Z * scLANE_se_link))
                           if (!is.null(size.factor.offset)) {
                             gene_df <- dplyr::mutate(gene_df,
                                                      dplyr::across(c(rna, scLANE_pred, scLANE_ci_ll, scLANE_ci_ul), \(m) m * size_factor))
                           }
                           if (log1p.norm) {
                             gene_df <- dplyr::mutate(gene_df,
                                                      rna_log1p = log1p(rna),
                                                      scLANE_pred_log1p = log1p(scLANE_pred),
                                                      scLANE_ci_ll_log1p = log1p(scLANE_ci_ll),
                                                      scLANE_ci_ul_log1p = log1p(scLANE_ci_ul))
                           }
                           gene_df <- dplyr::bind_cols(gene_df, cell.meta.data[!is.na(x), , drop = FALSE])
                           return(gene_df)
                         }) %>%
              purrr::reduce(rbind)
    if (!is.null(id.vec)) {
      mod_df <- dplyr::mutate(mod_df, subj_id = id.vec[!is.na(x)])
    } else {
      mod_df <- dplyr::mutate(mod_df, subj_id = NA_character_)
    }
    return(mod_df)
  })
  final_df <- purrr::reduce(mod_df_list, rbind) %>%
              dplyr::relocate(cell, subj_id, lineage)
  if (is.null(id.vec)) {
    final_df <- dplyr::select(final_df, -subj_id)
  }
  if (!is.null(filter.lineage)) {
    final_df <- dplyr::filter(final_df, !lineage %in% filter.lineage)
  }
  if (length(genes) > 1) {
    final_df <- dplyr::mutate(final_df,
                              gene = as.character(gene),
                              gene = factor(gene, levels = genes))
  }
  return(final_df)
}
