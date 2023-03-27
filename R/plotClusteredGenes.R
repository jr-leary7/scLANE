#' Plot the results from \code{\link{clusterGenes}}.
#'
#' @name plotClusteredGenes
#' @author Jack Leary
#' @description Plot per-lineage, per-cluster fitted values from \code{scLANE}.
#' @import magrittr
#' @importFrom purrr imap reduce
#' @importFrom dplyr inner_join rename mutate
#' @param test.dyn.results The list returned by \code{\link{testDynamic}} - no extra processing required. Defaults to NULL.
#' @param gene.clusters The data.frame returned by \code{\link{clusterGenes}}. Defaults to NULL.
#' @param pt A data.frame containing the pseudotime or latent time estimates for each cell. Defaults to NULL.
#' @param parallel.exec Should \code{furrr} be used to speed up execution at the cost of increased memory usage? Defaults to TRUE.
#' @param n.cores If parallel execution is desired, how many cores should be utilized? Defaults to 2.
#' @return A \code{data.frame} object with ready-to-plot tidy data.
#' @seealso \code{\link{clusterGenes}}
#' @export
#' @examples
#' \dontrun{
#' plotClusteredGenes(test.dyn.results = gene_stats,
#'                    gene.clusters = gene_clusters,
#'                    pt = pt_df) %>%
#'   ggplot(aes(x = PT, y = FITTED, color = CLUSTER, group = GENE) +
#'   geom_line() +
#'   facet_wrap(~LINEAGE))
#' }

plotClusteredGenes <- function(test.dyn.results = NULL,
                               gene.clusters = NULL,
                               pt = NULL,
                               parallel.exec = TRUE,
                               n.cores = 2) {
  # check inputs
  if (is.null(test.dyn.results) || is.null(gene.clusters) || is.null(pt)) { stop("Arguments to plotClusteredGenes() are missing.") }
  colnames(pt) <- paste0("Lineage_", LETTERS[1:ncol(pt)])
  if (parallel.exec) {
    future::plan(future::multisession, workers = n.cores)
    furrr::future_imap(test.dyn.results, function(x, y) {
      df_list <- vector("list", ncol(pt))
      for (l in seq(ncol(pt))) {
        lineage_name <- colnames(pt)[l]
        if (grepl("MARGE model error", x[[lineage_name]]$Model_Status)) {
          fitted_vals_mat <- data.frame(GENE = character(),
                                        LINEAGE = character(),
                                        CELL = character(),
                                        FITTED_LINK = numeric(),
                                        FITTED = numeric(),
                                        # EXP = numeric(),
                                        PT = numeric())
        } else {
          fitted_vals_mat <- data.frame(GENE = y,
                                        LINEAGE = LETTERS[l],
                                        CELL = rownames(pt)[!is.na(pt[, l])],
                                        FITTED_LINK = x[[lineage_name]]$MARGE_Preds$marge_link_fit,
                                        FITTED = exp(x[[lineage_name]]$MARGE_Preds$marge_link_fit),
                                        # EXP = gene.counts[!is.na(pt[, l]), y],  -- maybe bring this back in later, but not needed right now
                                        PT = pt[!is.na(pt[, l]), l])
        }
        df_list[[l]] <- fitted_vals_mat
      }
      df_temp <- purrr::reduce(df_list, rbind)
      return(df_temp)
    }) -> all_genes
  } else {
    purrr::imap(test.dyn.results, function(x, y) {
      df_list <- vector("list", ncol(pt))
      for (l in seq(ncol(pt))) {
        lineage_name <- colnames(pt)[l]
        if (grepl("MARGE model error", x[[lineage_name]]$Model_Status)) {
          fitted_vals_mat <- data.frame(GENE = character(),
                                        LINEAGE = character(),
                                        CELL = character(),
                                        FITTED_LINK = numeric(),
                                        FITTED = numeric(),
                                        # EXP = numeric(),
                                        PT = numeric())
        } else {
          fitted_vals_mat <- data.frame(GENE = y,
                                        LINEAGE = LETTERS[l],
                                        CELL = rownames(pt)[!is.na(pt[, l])],
                                        FITTED_LINK = x[[lineage_name]]$MARGE_Preds$marge_link_fit,
                                        FITTED = exp(x[[lineage_name]]$MARGE_Preds$marge_link_fit),
                                        # EXP = gene.counts[!is.na(pt[, l]), y],  -- maybe bring this back in later, but not needed right now
                                        PT = pt[!is.na(pt[, l]), l])
        }
        df_list[[l]] <- fitted_vals_mat
      }
      df_temp <- purrr::reduce(df_list, rbind)
      return(df_temp)
    }) -> all_genes
  }
  gene_res <- purrr::reduce(all_genes, rbind) %>%
              dplyr::inner_join(gene.clusters, by = c("GENE" = "Gene", "LINEAGE" = "Lineage")) %>%
              dplyr::rename(CLUSTER = Cluster) %>%
              dplyr::mutate(CLUSTER = as.factor(CLUSTER))
  return(gene_res)
}
