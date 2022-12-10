#' Cluster the fitted values from a set of \code{scLANE} models.
#'
#' @name clusterGenes
#' @author Jack Leary
#' @description This function takes as input the output from \code{\link{testDynamic}} and clusters the fitted values from the model for each gene using one of several user-chosen algorithms.
#' @importFrom purrr map discard map2 reduce
#' @importFrom stats hclust cutree kmeans
#' @importFrom rlang is_na
#' @importFrom magrittr set_rownames
#' @importFrom irlba prcomp_irlba
#' @importFrom cluster silhouette
#' @importFrom bluster makeSNNGraph
#' @importFrom BiocNeighbors AnnoyParam
#' @importFrom igraph cluster_leiden
#' @param test.dyn.results The list returned by \code{\link{testDynamic}} - no extra processing required. Defaults to NULL.
#' @param clust.algo The clustering method to use. Can be one of "hclust", "kmeans", "leiden". Defaults to "leiden".
#' @param use.pca Should PCA be performed prior to clustering? Defaults to FALSE.
#' @param n.PC The number of principal components to use when performing dimension reduction prior to clustering. Defaults to 15.
#' @param lineages Should one or more lineages be isolated? If so, specify which one(s). Otherwise, all lineages will be clustered independently. Defaults to NULL.
#' @return A data.frame of with three columns: \code{Gene}, \code{Lineage}, and \code{Cluster}.
#' @seealso \code{\link{testDynamic}}
#' @seealso \code{\link{plotClusteredGenes}}
#' @export
#' @examples
#' \dontrun{
#' clusterGenes(test.dyn.results = gene_stats, clust.algo = "leiden")
#' clusterGenes(test.dyn.results = gene_stats,
#'              use.PCA = TRUE,
#'              n.PC = 10,
#'              lineages = "B")
#' clusterGenes(test.dyn.results = gene_stats, lineages = c("A", "C"))
#' }

clusterGenes <- function(test.dyn.results = NULL,
                         clust.algo = "leiden",
                         use.pca = FALSE,
                         n.PC = 15,
                         lineages = NULL) {
  # check inputs
  if (is.null(test.dyn.results)) { stop("test.dyn.results must be supplied to clusterGenes().") }
  clust.algo <- tolower(clust.algo)
  if (!clust.algo %in% c("hclust", "kmeans", "leiden")) { stop("clust.algo must be one of 'hclust', 'kmeans', or 'leiden'.") }
  if ((use.pca & is.null(n.PC)) | (use.pca & n.PC <= 0)) { stop("n.PC must be a non-zero integer when clustering on principal components.") }
  if (is.null(lineages)) {
    lineages <- LETTERS[1:length(test.dyn.results[[1]])]
  }
  gene_cluster_list <- vector("list", length = length(lineages))
  for (l in seq_along(lineages)) {
    # coerce fitted values to a gene x cell matrix
    lineage_name <- paste0("Lineage_", lineages[l])
    fitted_vals_mat <- purrr::map(test.dyn.results, function(x) { x[[lineage_name]]$MARGE_Preds}) %>%
                       setNames(names(test.dyn.results)) %>%
                       purrr::discard(rlang::is_na) %>%
                       purrr::map2(.y = names(.), function(x, y) {
                         t(as.data.frame(exp(x$marge_link_fit))) %>%
                           magrittr::set_rownames(y)
                       }) %>%
                       purrr::reduce(rbind)
    if (use.pca) {
      fitted_vals_pca <- irlba::prcomp_irlba(fitted_vals_mat,
                                             n = n.PC,
                                             center = TRUE,
                                             scale. = TRUE)
    }
    # hierarchical clustering routine w/ Ward's linkage
    if (clust.algo  == "hclust") {
      if (use.pca) {
        hclust_tree <- stats::hclust(dist(fitted_vals_pca$x), method = "ward.D2")
      } else {
        hclust_tree <- stats::hclust(dist(fitted_vals_mat), method = "ward.D2")
      }
      k_vals <- c(2:10)
      sil_vals <- vector("numeric", 9L)
      for (k in seq_along(k_vals)) {
        clust_res <- stats::cutree(hclust_tree, k = k_vals[k])
        if (use.pca) {
          sil_res <- cluster::silhouette(clust_res, dist(fitted_vals_pca$x))
        } else {
          sil_res <- cluster::silhouette(clust_res, dist(fitted_vals_mat))
        }
        sil_vals[k] <- mean(sil_res[, 3])  # silhouette widths stored in third column
      }
      k_to_use <- k_vals[which.max(sil_vals)]
      clust_res <- stats::cutree(hclust_tree, k = k_to_use)
      gene_clusters <- data.frame(Gene = rownames(fitted_vals_mat),
                                  Lineage = lineages[l],
                                  Cluster = clust_res)
    # k-means clustering routine w/ Hartigan-Wong algorithm
    } else if (clust.algo == "kmeans") {
      k_vals <- c(2:10)
      sil_vals <- vector("numeric", 9L)
      for (k in seq_along(k_vals)) {
        if (use.pca) {
          clust_res <- stats::kmeans(fitted_vals_pca$x,
                                     centers = k_vals[k],
                                     nstart = 5,
                                     algorithm = "Hartigan-Wong")
          sil_res <- cluster::silhouette(clust_res$cluster, dist(fitted_vals_pca$x))
        } else {
          clust_res <- stats::kmeans(fitted_vals_mat,
                                     centers = k_vals[k],
                                     nstart = 5,
                                     algorithm = "Hartigan-Wong")
          sil_res <- cluster::silhouette(clust_res$cluster, dist(fitted_vals_mat))
        }
        sil_vals[k] <- mean(sil_res[, 3])  # silhouette widths stored in third column
      }
      k_to_use <- k_vals[which.max(sil_vals)]
      if (use.pca) {
        clust_res <- stats::kmeans(fitted_vals_pca$x,
                                   centers = k_to_use,
                                   nstart = 5,
                                   algorithm = "Hartigan-Wong")
      } else {
        clust_res <- stats::kmeans(fitted_vals_mat,
                                   centers = k_to_use,
                                   nstart = 5,
                                   algorithm = "Hartigan-Wong")
      }
      gene_clusters <- data.frame(Gene = rownames(fitted_vals_mat),
                                  Lineage = lineages[l],
                                  Cluster = clust_res$cluster)
    # Leiden clustering routine
    } else if (clust.algo == "leiden") {
      if (use.pca) {
        fitted_vals_graph <- bluster::makeSNNGraph(x = fitted_vals_pca$x,
                                                   k = 30,
                                                   type = "jaccard",
                                                   BNPARAM = BiocNeighbors::AnnoyParam(distance = "Euclidean"))
      } else {
        fitted_vals_graph <- bluster::makeSNNGraph(x = fitted_vals_mat,
                                                   k = 30,
                                                   type = "jaccard",
                                                   BNPARAM = BiocNeighbors::AnnoyParam(distance = "Euclidean"))
      }
      res_vals <- seq(0.1, 1, by = 0.1)
      quality_vals <- vector("numeric", 10L)
      for (r in seq_along(res_vals)) {
        clust_res <- igraph::cluster_leiden(graph = fitted_vals_graph,
                                            objective_function = "modularity",
                                            resolution_parameter = res_vals[r])
        quality_vals[r] <- clust_res$quality  # not sure what metric this uses lol
      }
      res_to_use <- res_vals[which.max(quality_vals)]
      clust_res <- igraph::cluster_leiden(graph = fitted_vals_graph,
                                          objective_function = "modularity",
                                          resolution_parameter = res_to_use)
      gene_clusters <- data.frame(Gene = rownames(fitted_vals_mat),
                                  Lineage = lineages[l],
                                  Cluster = clust_res$membership)
    }
    gene_cluster_list[[l]] <- gene_clusters
  }
  res <- purrr::reduce(gene_cluster_list, rbind)
  return(res)
}
