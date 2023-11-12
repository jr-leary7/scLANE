#' Cluster the fitted values from a set of \code{scLANE} models.
#'
#' @name clusterGenes
#' @author Jack Leary
#' @description This function takes as input the output from \code{\link{testDynamic}} and clusters the fitted values from the model for each gene using one of several user-chosen algorithms. An approximately optimal clustering is determined by iterating over reasonable hyperparameter values & choosing the value with the highest mean silhouette score based on the cosine distance.
#' @import magrittr
#' @importFrom purrr map discard map2 reduce
#' @importFrom stats setNames hclust cutree kmeans as.dist
#' @param test.dyn.res The list returned by \code{\link{testDynamic}} - no extra processing required. Defaults to NULL.
#' @param pt A data.frame containing the pseudotime or latent time estimates for each cell. Defaults to NULL.
#' @param size.factor.offset (Optional) An offset to be used to rescale the fitted values. Can be generated easily with \code{\link{createCellOffset}}. No need to provide if the GEE backend was used. Defaults to NULL.
#' @param clust.algo The clustering method to use. Can be one of "hclust", "kmeans", "leiden". Defaults to "leiden".
#' @param use.pca Should PCA be performed prior to clustering? Defaults to FALSE.
#' @param n.PC The number of principal components to use when performing dimension reduction prior to clustering. Defaults to 15.
#' @param lineages Should one or more lineages be isolated? If so, specify which one(s). Otherwise, all lineages will be clustered independently. Defaults to NULL.
#' @details
#' \itemize{
#' \item Due to some peculiarities of how the fitted values (on the link scale) are generated for \code{geeM} models, it's not necessary to multiply them by the offset as this is done internally. For GLM & GEE models, the opposite is true, and \code{size.factor.offset} must be provided in order to rescale the fitted values correctly.
#' }
#' @return A data.frame of with three columns: \code{Gene}, \code{Lineage}, and \code{Cluster}.
#' @seealso \code{\link{testDynamic}}
#' @seealso \code{\link{embedGenes}}
#' @seealso \code{\link{plotClusteredGenes}}
#' @export
#' @examples
#' data(sim_pseudotime)
#' data(scLANE_models)
#' cell_offset <- createCellOffset(sim_counts)
#' gene_clusters <- clusterGenes(scLANE_models,
#'                               pt = sim_pseudotime,
#'                               size.factor.offset = cell_offset)

clusterGenes <- function(test.dyn.res = NULL,
                         pt = NULL,
                         size.factor.offset = NULL,
                         clust.algo = "leiden",
                         use.pca = FALSE,
                         n.PC = 15L,
                         lineages = NULL) {
  # check inputs
  if (is.null(test.dyn.res) || is.null(pt)) { stop("test.dyn.res & pt must be supplied to clusterGenes().") }
  clust.algo <- tolower(clust.algo)
  if (!clust.algo %in% c("hclust", "kmeans", "leiden")) { stop("clust.algo must be one of 'hclust', 'kmeans', or 'leiden'.") }
  if ((use.pca & is.null(n.PC)) || (use.pca & n.PC <= 0)) { stop("n.PC must be a non-zero integer when clustering on principal components.") }
  if (is.null(lineages)) {
    lineages <- LETTERS[seq_along(test.dyn.res[[1]])]
  }
  gene_cluster_list <- vector("list", length = length(lineages))
  for (l in seq_along(lineages)) {
    # coerce fitted values to a gene x cell matrix, dropping genes w/ model errors
    lineage_name <- paste0("Lineage_", lineages[l])
    fitted_vals_mat <- purrr::map(test.dyn.res, \(x) x[[lineage_name]]$MARGE_Preds) %>%
                       stats::setNames(names(test.dyn.res)) %>%
                       purrr::discard(rlang::is_na) %>%
                       purrr::discard(rlang::is_null) %>%
                       purrr::discard(\(p) rlang::inherits_only(p, "try-error")) %>%
                       purrr::map2(.y = names(.), function(x, y) {
                         if (is.null(size.factor.offset)) {
                           t(as.data.frame(exp(x$marge_link_fit))) %>%
                             magrittr::set_rownames(y)
                         } else {
                           t(as.data.frame(exp(x$marge_link_fit)) * unname(size.factor.offset)[!is.na(pt[, l])]) %>%
                             magrittr::set_rownames(y)
                         }
                       }) %>%
                       purrr::reduce(rbind)
    if (use.pca) {
      fitted_vals_pca <- irlba::prcomp_irlba(fitted_vals_mat,
                                             n = n.PC,
                                             center = TRUE,
                                             scale. = TRUE)
      dist_matrix <- stats::as.dist(1 - coop::tcosine(x = fitted_vals_pca$x))
    } else {
      dist_matrix <- stats::as.dist(1 - coop::tcosine(x = fitted_vals_mat))
    }
    # hierarchical clustering routine w/ Ward's linkage
    if (clust.algo  == "hclust") {
      hclust_tree <- stats::hclust(dist_matrix, method = "ward.D2")
      k_vals <- seq(2, 10)
      sil_vals <- vector("numeric", 9L)
      clust_list <- vector("list", 9L)
      for (k in seq_along(k_vals)) {
        clust_res <- stats::cutree(hclust_tree, k = k_vals[k])
        sil_res <- cluster::silhouette(clust_res, dist_matrix)
        sil_vals[k] <- mean(sil_res[, 3])
        clust_list[[k]] <- clust_res
      }
      gene_clusters <- data.frame(Gene = rownames(fitted_vals_mat),
                                  Lineage = lineages[l],
                                  Cluster = clust_list[[which.max(sil_vals)]])
    # k-means clustering routine w/ Hartigan-Wong algorithm
    } else if (clust.algo == "kmeans") {
      k_vals <- seq(2, 10)
      sil_vals <- vector("numeric", 9L)
      clust_list <- vector("list", 9L)
      for (k in seq_along(k_vals)) {
        if (use.pca) {
          clust_res <- stats::kmeans(fitted_vals_pca$x,
                                     centers = k_vals[k],
                                     nstart = 5,
                                     algorithm = "Hartigan-Wong")
        } else {
          clust_res <- stats::kmeans(fitted_vals_mat,
                                     centers = k_vals[k],
                                     nstart = 5,
                                     algorithm = "Hartigan-Wong")
        }
        clust_list[[k]] <- clust_res
        sil_res <- cluster::silhouette(clust_res$cluster, dist_matrix)
        sil_vals[k] <- mean(sil_res[, 3])
      }
      gene_clusters <- data.frame(Gene = rownames(fitted_vals_mat),
                                  Lineage = lineages[l],
                                  Cluster = clust_list[[which.max(sil_vals)]]$cluster)
    # Leiden clustering routine
    } else if (clust.algo == "leiden") {
      if (use.pca) {
        fitted_vals_graph <- bluster::makeSNNGraph(x = fitted_vals_pca$x,
                                                   k = 20,
                                                   type = "jaccard",
                                                   BNPARAM = BiocNeighbors::AnnoyParam(distance = "Cosine"))
      } else {
        fitted_vals_graph <- bluster::makeSNNGraph(x = fitted_vals_mat,
                                                   k = 20,
                                                   type = "jaccard",
                                                   BNPARAM = BiocNeighbors::AnnoyParam(distance = "Cosine"))
      }
      res_vals <- seq(0.1, 1, by = 0.1)
      sil_vals <- vector("numeric", 10L)
      clust_list <- vector("list", 10L)
      for (r in seq_along(res_vals)) {
        clust_res <- igraph::cluster_leiden(graph = fitted_vals_graph,
                                            objective_function = "modularity",
                                            resolution_parameter = res_vals[r])
        if (clust_res$nb_clusters == 1) {
          sil_vals[r] <- 0
        } else {
          sil_res <- cluster::silhouette(clust_res$membership, dist_matrix)
          sil_vals[r] <- mean(sil_res[, 3])
        }
        clust_list[[r]] <- clust_res$membership
      }
      gene_clusters <- data.frame(Gene = rownames(fitted_vals_mat),
                                  Lineage = lineages[l],
                                  Cluster = clust_list[[which.max(sil_vals)]])
    }
    gene_cluster_list[[l]] <- gene_clusters
  }
  res <- purrr::reduce(gene_cluster_list, rbind)
  return(res)
}
