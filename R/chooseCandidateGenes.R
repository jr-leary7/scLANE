#' Choose candidate genes for trajectory DE analysis.
#'
#' @name chooseCandidateGenes
#' @author Jack R. Leary
#' @description This function identifies good gene candidates for trajectory differential expression modeling by ranking genes based on their mean expression, SD of expression, and sparsity across cells.
#' @importFrom purrr map reduce
#' @importFrom Matrix Matrix rowMeans
#' @importFrom dplyr with_groups summarise arrange desc mutate row_number rowwise c_across ungroup slice_head pull
#' @importFrom tidyselect starts_with
#' @param obj An object of class \code{\link[SingleCellExperiment]{SingleCellExperiment}}, \code{\link[Seurat]{Seurat}}, or \code{CellDataSet}, or a gene-by-cell matrix. Defaults to NULL.
#' @param group.by.subject Boolean specifying whether or not the summary statistics should be computed per-subject and then mean-aggregated. Defaults to TRUE.
#' @param id.vec A vector of subject IDs. Defaults to NULL.
#' @param n.desired.genes An integer specifying the number of candidate genes to return. Defaults to 2000.
#' @return A vector of candidate gene names.
#' @export
#' @examples
#' data(sim_counts)
#' candidate_genes <- chooseCandidateGenes(sim_counts, id.vec = sim_counts$subject)

chooseCandidateGenes <- function(obj = NULL,
                                 group.by.subject = TRUE,
                                 id.vec = NULL,
                                 n.desired.genes = 2000L) {
  # check inputs
  if (is.null(obj) || !(inherits(obj, "SingleCellExperiment") || inherits(obj, "Seurat"))) { stop("Please provide a SingleCellExperiment or Seurat object.") }
  if (group.by.subject && is.null(id.vec)) { stop("Grouping by subject requires a vector of subject IDs.") }
  # extract counts matrix from SingleCellExperiment, Seurat, or CellDataSet object 
  if (inherits(obj, "SingleCellExperiment")) {
    counts_matrix <- BiocGenerics::counts(obj)
  } else if (inherits(obj, "Seurat")) {
    counts_matrix <- Seurat::GetAssayData(obj,
                                          slot = "counts",
                                          assay = Seurat::DefaultAssay(obj))
  } else if (inherits(expr.mat, "cell_data_set")) {
    counts_matrix <- BiocGenerics::counts(slot)
  } else if (inherits(obj, "dgCMatrix") || inherits(obj, "dgRMatrix") || inherits(obj, "matrix")) {
    counts_matrix <- obj
  }
  if (!inherits(counts_matrix, "dgCMatrix")) {
    counts_matrix <- Matrix::Matrix(counts_matrix, sparse = TRUE)
  }
  # compute gene summary statistics
  if (group.by.subject) {
    grouped_stats <- purrr::map(seq(unique(id.vec)), \(i) {
      sub_matrix <- counts_matrix[, which(id.vec == unique(id.vec)[i])]
      gene_means <- Matrix::rowMeans(sub_matrix)
      gene_sds <- sqrt(Matrix::rowMeans((sub_matrix - gene_means)^2))
      gene_sparsity <- Matrix::rowMeans(sub_matrix == 0)
      res <- data.frame(subject = unique(id.vec)[i],
                        gene = rownames(sub_matrix),
                        mu = unname(gene_means),
                        sigma = unname(gene_sds),
                        lambda = unname(gene_sparsity))
      return(res)
    })
    gene_df <- purrr::reduce(grouped_stats, rbind) %>%
               dplyr::with_groups(gene,
                                  dplyr::summarise,
                                  mu = mean(mu),
                                  sigma = mean(sigma),
                                  lambda = mean(lambda))

  } else {
    gene_means <- Matrix::rowMeans(counts_matrix)
    gene_sds <- sqrt(Matrix::rowMeans((counts_matrix - gene_means)^2))
    gene_sparsity <- Matrix::rowMeans(counts_matrix == 0)
    gene_df <- data.frame(gene = rownames(counts_matrix),
                          mu = unname(gene_means),
                          sigma = unname(gene_sds),
                          lambda = unname(gene_sparsity))
  }
  # rank genes
  gene_df <- dplyr::arrange(gene_df, dplyr::desc(mu)) %>%
             dplyr::mutate(rank_mu = dplyr::row_number()) %>%
             dplyr::arrange(dplyr::desc(sigma)) %>%
             dplyr::mutate(rank_sigma = dplyr::row_number()) %>%
             dplyr::arrange(lambda) %>%
             dplyr::mutate(rank_lambda = dplyr::row_number()) %>%
             dplyr::rowwise() %>%
             dplyr::mutate(mean_rank = mean(dplyr::c_across(tidyselect::starts_with("rank")))) %>%
             dplyr::ungroup() %>%
             dplyr::arrange(mean_rank)
  # select candidates
  selected_genes <- dplyr::slice_head(gene_df, n = n.desired.genes) %>%
                    dplyr::pull(gene)
  return(selected_genes)
}
