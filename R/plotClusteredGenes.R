#' Generate tidy results from \code{\link{clusterGenes}} to use in plotting.
#'
#' @name plotClusteredGenes
#' @author Jack Leary
#' @description Generate a table of per-lineage, per-cluster fitted values from \code{scLANE} to be used in visualizations.
#' @import magrittr
#' @importFrom future plan multisession sequential
#' @importFrom furrr future_imap
#' @importFrom purrr reduce
#' @importFrom dplyr inner_join rename mutate
#' @param test.dyn.res The list returned by \code{\link{testDynamic}} - no extra processing required. Defaults to NULL.
#' @param gene.clusters The data.frame returned by \code{\link{clusterGenes}}. Defaults to NULL.
#' @param pt A data.frame containing the pseudotime or latent time estimates for each cell. Defaults to NULL.
#' @param size.factor.offset (Optional) An offset to be used to rescale the fitted values. Can be generated easily with \code{\link{createCellOffset}}. No need to provide if the GEE backend was used. Defaults to NULL.
#' @param n.cores If parallel execution is desired, how many cores should be utilized? Defaults to 2.
#' @details
#' \itemize{
#' \item Due to some peculiarities of how the fitted values (on the link scale) are generated for \code{geeM} models, it's not necessary to multiply them by the offset as this is done internally. For GLM & GEE models, the opposite is true, and \code{size.factor.offset} must be provided in order to rescale the fitted values correctly.
#' }
#' @return A \code{data.frame} with ready-to-plot tidy data. Includes columns for gene name, pseudotime lineage, cell name, fitted values on link & response scale, pseudotime, & gene cluster.
#' @seealso \code{\link{clusterGenes}}
#' @export
#' @examples
#' \dontrun{
#' plotClusteredGenes(test.dyn.res = gene_stats,
#'                    gene.clusters = gene_clusters,
#'                    pt = pt_df) %>%
#'   ggplot(aes(x = PT, y = FITTED, color = CLUSTER, group = GENE)) +
#'   facet_wrap(~LINEAGE + CLUSTER) +
#'   geom_line() +
#'   theme_classic()
#' }

plotClusteredGenes <- function(test.dyn.res = NULL,
                               gene.clusters = NULL,
                               pt = NULL,
                               size.factor.offset = NULL,
                               n.cores = 2L) {
  # check inputs
  if (is.null(test.dyn.res) || is.null(gene.clusters) || is.null(pt)) { stop("Arguments to plotClusteredGenes() are missing.") }
  colnames(pt) <- paste0("Lineage_", LETTERS[seq_len(ncol(pt))])
  if (n.cores > 1L) {
    future::plan(future::multisession, workers = n.cores)
  } else {
    future::plan(future::sequential)
  }
  furrr::future_imap(test.dyn.res, function(x, y) {
    df_list <- vector("list", ncol(pt))
    for (l in seq_len(ncol(pt))) {
      lineage_name <- colnames(pt)[l]
      if (grepl("MARGE model error", x[[lineage_name]]$Model_Status)) {
        fitted_vals_mat <- data.frame(GENE = character(),
                                      LINEAGE = character(),
                                      CELL = character(),
                                      FITTED_LINK = numeric(),
                                      FITTED = numeric(),
                                      PT = numeric())
      } else {
        fitted_vals_mat <- data.frame(GENE = y,
                                      LINEAGE = LETTERS[l],
                                      CELL = rownames(pt)[!is.na(pt[, l])],
                                      FITTED_LINK = x[[lineage_name]]$MARGE_Preds$marge_link_fit,
                                      FITTED = exp(x[[lineage_name]]$MARGE_Preds$marge_link_fit),
                                      PT = pt[!is.na(pt[, l]), l])
        if (!is.null(size.factor.offset)) {
          fitted_vals_mat <- dplyr::mutate(fitted_vals_mat, FITTED = FITTED * unname(size.factor.offset)[!is.na(pt[, l])])
        }
      }
      df_list[[l]] <- fitted_vals_mat
    }
    df_temp <- purrr::reduce(df_list, rbind)
    return(df_temp)
  }) -> all_genes
  gene_res <- purrr::reduce(all_genes, rbind) %>%
              dplyr::inner_join(gene.clusters,
                                by = c("GENE" = "Gene", "LINEAGE" = "Lineage")) %>%
              dplyr::rename(CLUSTER = Cluster) %>%
              dplyr::mutate(CLUSTER = as.factor(CLUSTER))
  return(gene_res)
}
