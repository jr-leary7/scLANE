#' Generate PCA & UMAP embeddings of fitted gene dynamics.
#'
#' @name embedGenes
#' @author Jack R. Leary
#' @import magrittr
#' @importFrom purrr map map_dbl
#' @importFrom dplyr bind_cols
#' @importFrom stats as.dist
#' @description Embed genes in dimension-reduced space given a smoothed counts matrix.
#' @param smoothed.counts The output from \code{\link{smoothedCountsMatrix}}. Defaults to NULL.
#' @param genes A character vector of genes to embed. If not specified, all genes in \code{smoothed.counts} are used. Defaults to NULL.
#' @param pca.init A boolean specifying whether or not the embedded PCs should be used as initialization for clustering and UMAP. The default is to cluster / embed the raw dynamics i.e., defaults to FALSE.
#' @param pc.embed (Optional) How many PCs should be used to cluster the genes and run UMAP? Defaults to 30.
#' @param pc.return (Optional) How many principal components should be included in the output? Defaults to 2.
#' @param cluster.genes (Optional) Should genes be clustered in PCA space using the Leiden algorithm? Defaults to TRUE.
#' @param gene.meta.data (Optional) A data.frame of metadata values for each gene (HVG status, Ensembl ID, gene biotype, etc.) that will be included in the result table. Defaults to NULL.
#' @param k.param (Optional) The value of nearest-neighbors used in creating the SNN graph prior to clustering & in running UMAP. Defaults to 20.
#' @param resolution.param (Optional) The value of the resolution parameter for the Leiden algorithm. If unspecified, silhouette scoring is used to select an optimal value. Defaults to NULL.
#' @param random.seed (Optional) The random seed used to control stochasticity in the clustering algorithm. Defaults to 312.
#' @param n.cores (Optional) Integer specifying the number of threads used by \code{\link[uwot]{umap}} and in \code{\link[bluster]{makeSNNGraph}}. Defaults to 2.
#' @return A data.frame containing embedding coordinates, cluster IDs, and metadata for each gene.
#' @export
#' @examples
#' data(sim_pseudotime)
#' data(scLANE_models)
#' smoothed_dynamics <- smoothedCountsMatrix(scLANE_models,
#'                                           pt = sim_pseudotime,
#'                                           n.cores = 1L)
#' gene_embed <- embedGenes(smoothed_dynamics$Lineage_A, n.cores = 1L)

embedGenes <- function(smoothed.counts = NULL,
                       genes = NULL,
                       pca.init = FALSE,
                       pc.embed = 30,
                       pc.return = 2,
                       cluster.genes = TRUE,
                       gene.meta.data = NULL,
                       k.param = 20,
                       resolution.param = NULL,
                       random.seed = 312,
                       n.cores = 2L) {
  # check inputs
  if (is.null(smoothed.counts)) { stop("You forgot to provide a smoothed counts matrix to embedGenes().") }
  genes <- colnames(smoothed.counts)
  smoothed.counts <- t(smoothed.counts)
  # embeddings
  smoothed_counts_pca <- irlba::prcomp_irlba(smoothed.counts,
                                             n = pc.embed,
                                             center = TRUE,
                                             scale. = TRUE)
  if (pca.init) {
    smoothed_counts_umap <- uwot::umap(smoothed_counts_pca$x,
                                       n_components = 2,
                                       metric = "cosine",
                                       n_neighbors = k.param,
                                       init = "spectral",
                                       nn_method = "annoy",
                                       seed = random.seed,
                                       n_threads = n.cores,
                                       verbose = FALSE)
  } else {
    smoothed_counts_umap <- uwot::umap(smoothed.counts,
                                       n_components = 2,
                                       metric = "cosine",
                                       n_neighbors = k.param,
                                       init = "spectral",
                                       nn_method = "annoy",
                                       seed = random.seed,
                                       n_threads = n.cores,
                                       verbose = FALSE)
  }
  # clustering w/ silhouette score parameter tuning
  if (cluster.genes) {
    if (pca.init) {
      smoothed_counts_snn <- bluster::makeSNNGraph(smoothed_counts_pca$x,
                                                   k = k.param,
                                                   type = "jaccard",
                                                   BNPARAM = BiocNeighbors::AnnoyParam(distance = "Cosine"),
                                                   BPPARAM = BiocParallel::SnowParam(workers = n.cores, RNGseed = random.seed))
    } else {
      smoothed_counts_snn <- bluster::makeSNNGraph(smoothed.counts,
                                                   k = k.param,
                                                   type = "jaccard",
                                                   BNPARAM = BiocNeighbors::AnnoyParam(distance = "Cosine"),
                                                   BPPARAM = BiocParallel::SnowParam(workers = n.cores, RNGseed = random.seed))
    }
    if (is.null(resolution.param)) {
      if (pca.init) {
        dist_matrix <- stats::as.dist(1 - coop::tcosine(x = smoothed_counts_pca$x))
      } else {
        dist_matrix <- stats::as.dist(1 - coop::tcosine(x = smoothed.counts))
      }
      clust_runs <- purrr::map(seq(0.1, 0.7, by = 0.1), \(r) {
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
      cluster_vec <- clust_runs[[best_clustering]]$clusters
    } else {
      smoothed_counts_clust <- igraph::cluster_leiden(smoothed_counts_snn,
                                                      objective_function = "modularity",
                                                      resolution_parameter = resolution.param)
      cluster_vec <- as.factor(smoothed_counts_clust$membership - 1L)
    }
  } else {
    cluster_vec <- NA_integer_
  }
  # prepare results
  pca_df <- as.data.frame(smoothed_counts_pca$x[, seq(pc.return)])
  colnames(pca_df) <- paste0("pc", seq(pc.return))
  gene_df <- data.frame(gene = genes,
                        leiden = cluster_vec,
                        umap1 = smoothed_counts_umap[, 1],
                        umap2 = smoothed_counts_umap[, 2])
  gene_df <- dplyr::bind_cols(gene_df, pca_df)
  if (!cluster.genes) {
    gene_df <- dplyr::select(gene_df, -leiden)
  }
  return(gene_df)
}
