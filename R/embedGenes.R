#' Generate a table of fitted values and celltype metadata for genes of interest.
#'
#' @name embedGenes
#' @author Jack Leary
#' @import magrittr
#' @importFrom purrr map map_dbl
#' @importFrom dplyr bind_cols
#' @importFrom stats as.dist
#' @description Embed genes in dimension-reduced space given a smoothed counts matrix.
#' @param smoothed.counts The output from \code{\link{smoothedCountsMatrix}}. Defaults to NULL.
#' @param genes A character vector of genes to embed. If not specified, all genes in \code{test.dyn.res} are used. Defaults to NULL.
#' @param pc.embed (Optional) How many PCs should be used to cluster the genes and run UMAP? Defaults to 30.
#' @param pcs.return (Optional) How many principal components should be included in the output? Defaults to 2.
#' @param cluster.genes (Optional) Should genes be clustered in PCA space using the Leiden algorithm? Defaults to TRUE.
#' @param gene.meta.data (Optional) A data.frame of metadata values for each gene (HVG status, Ensembl ID, gene biotype, etc.) that will be included in the result table. Defaults to NULL.
#' @param k.param (Optional) The value of nearest-neighbors used in creating the SNN graph prior to clustering & in running UMAP. Defaults to 20.
#' @param random.seed (Optional) The random seed used to control stochasticity in the clustering algorithm. Defaults to 312.
#' @return A data.frame containing embedding coordinates, cluster IDs, and metadata for each gene.
#' @export
#' @examples
#' \dontrun{
#' embedGenes(smoothed_counts$Lineage_A,
#'            pcs.return = 3,
#'            cluster.genes = TRUE)
#' }

embedGenes <- function(smoothed.counts = NULL,
                       genes = NULL,
                       pc.embed = 30,
                       pcs.return = 2,
                       cluster.genes = TRUE,
                       gene.meta.data = NULL,
                       k.param = 20,
                       random.seed = 312) {
  # check inputs
  if (is.null(smoothed.counts)) { stop("You forgot to provide a smoothed counts matrix to embedGenes().") }
  # embeddings
  set.seed(random.seed)
  smoothed_counts_pca <- irlba::prcomp_irlba(t(smoothed.counts),
                                             pc.embed = pc.embed,
                                             center = TRUE,
                                             scale. = TRUE)
  smoothed_counts_umap <- uwot::umap(smoothed_counts_pca$x,
                                     n_components = 2,
                                     metric = "cosine",
                                     n_neighbors = k.param,
                                     init = "spectral")
  # clustering w/ silhouette score parameter tuning
  smoothed_counts_snn <- bluster::makeSNNGraph(smoothed_counts_pca$x,
                                               k = k.param,
                                               type = "jaccard",
                                               BNPARAM = BiocNeighbors::AnnoyParam(distance = "Cosine"))
  dist_matrix <- stats::as.dist(1 - coop::tcosine(x = smoothed_counts_pca$x))
  clust_runs <- purrr::map(c(0.1, 0.2, 0.3, 0.4, 0.5), \(r) {
    smoothed_counts_clust <- igraph::cluster_leiden(smoothed_counts_snn,
                                                    objective_function = "modularity",
                                                    resolution_parameter = r)
    if (smoothed_counts_clust$nb_clusters == 1) {
      sil_val <- 0
    } else {
      sil_val <- mean(cluster::silhouette(as.integer(smoothed_counts_clust$membership - 1L), dist_matrix)[, 3])
    }
    clust_res <- list(clusters = as.factor(smoothed_counts_clust$membership - 1L),
                      resolution = r,
                      silhouette = sil_val)
    return(clust_res)
  })
  best_clustering <- which.max(purrr::map_dbl(clust_runs, \(x) x$silhouette))
  # prepare results
  pca_df <- as.data.frame(smoothed_counts_pca$x[, seq(pcs.return)])
  colnames(pca_df) <- paste0("pc", seq(pcs.return))
  gene_clust_df <- data.frame(gene = colnames(smoothed.counts),
                              leiden = clust_runs[[best_clustering]]$clusters,
                              umap1 = smoothed_counts_umap[, 1],
                              umap2 = smoothed_counts_umap[, 2])
  gene_clust_df <- dplyr::bind_cols(gene_clust_df, pca_df)
  return(gene_clust_df)
}
