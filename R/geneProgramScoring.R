#' Add per-cell module scores for gene programs.
#'
#' @name geneProgramScoring
#' @author Jack Leary
#' @importFrom Matrix Matrix
#' @description This function uses \code{\link[UCell]{ScoreSignatures_UCell}} to create a per-cell module score for each of the provided gene clusters. If the input matrix is a \code{Seurat} or \code{SingleCellExperiment} object, then the resulting scores will be added to the \code{meta.data} or the \code{colData} slot, respectively. Otherwise, a data.frame of the per-program scores is returned.
#' @param expr.mat Either a \code{SingleCellExperiment} or \code{Seurat} object from which counts can be extracted, or a matrix of integer-valued counts with genes as rows & cells as columns. Defaults to NULL.
#' @param genes A character vector of gene IDs. Defaults to NULL.
#' @param gene.clusters A factor containing the cluster assignment of each gene in \code{genes}. Defaults to NULL.
#' @param program.labels (Optional) A character vector specifying a label for each gene cluster. Defaults to NULL.
#' @param minmax.norm (Optional) Should each program's score be min-max normalized to be on (0, 1)? Defaults to TRUE.
#' @param minmax.epsilon (Optional) The tolerance used to ensure that program scores equal to 0 or 1 do not occur. Defaults to 0.01.
#' @param n.cores (Optional) The number of cores used under the hood in \code{\link[UCell]{ScoreSignatures_UCell}}. Defaults to 2.
#' @return Either a \code{Seurat} or \code{SingleCellExperiment} object if \code{expr.mat} is in either form, or a data.frame containing per-cell program scores if \code{expr.mat} is a matrix.
#' @seealso \code{\link[UCell]{ScoreSignatures_UCell}}
#' @seealso \code{\link{geneProgramDrivers}}
#' @export
#' @examples
#' data(sim_counts)
#' data(scLANE_models)
#' data(sim_pseudotime)
#' smoothed_dynamics <- smoothedCountsMatrix(scLANE_models,
#'                                           pt = sim_pseudotime,
#'                                           n.cores = 1L)
#' gene_embed <- embedGenes(smoothed_dynamics$Lineage_A, n.cores = 1L)
#' sim_counts <- geneProgramScoring(sim_counts,
#'                                  genes = gene_embed$gene,
#'                                  gene.clusters = gene_embed$leiden,
#'                                  n.cores = 1L)

geneProgramScoring <- function(expr.mat = NULL,
                               genes = NULL,
                               gene.clusters = NULL,
                               program.labels = NULL,
                               minmax.norm = TRUE,
                               minmax.epsilon = 1e-2,
                               n.cores = 2L) {
  # check inputs
  if (is.null(expr.mat) || is.null(genes) || is.null(gene.clusters)) { stop("Arguments to geneProgramScoring() are missing.") }
  if (!is.factor(gene.clusters)) {
    gene.clusters <- as.factor(gene.clusters)
  }
  # set program labels
  cluster.labels <- unique(gene.clusters)
  if (is.null(program.labels)) {
    program.labels <- paste0("cluster_", cluster.labels)
  } else {
    program.labels <- gsub(" ", "_", program.labels)
    if (length(program.labels) != length(levels(gene.clusters))) {
      stop("Each cluster must have a label.")
    }
  }
  # set up query matrix
  if (inherits(expr.mat, "SingleCellExperiment")) {
    counts_matrix <- BiocGenerics::counts(expr.mat)
  } else if (inherits(expr.mat, "Seurat")) {
    counts_matrix <- Seurat::GetAssayData(expr.mat,
                                          slot = "counts",
                                          assay = Seurat::DefaultAssay(expr.mat))
  } else if (inherits(expr.mat, "matrix") || inherits(expr.mat, "array")) {
    counts_matrix <- Matrix::Matrix(expr.mat, sparse = TRUE)
  }
  # set up feature list
  program_list <- split(genes, gene.clusters)
  names(program_list) <- program.labels
  # run UCell
  program_scores <- UCell::ScoreSignatures_UCell(counts_matrix,
                                                 features = program_list,
                                                 ncores = n.cores)
  # min-max normalize with tol to (0, 1) if desired
  if (minmax.norm) {
    program_scores <- purrr::map(seq(ncol(program_scores)), \(i) {
      scores <- program_scores[, i]
      normed_scores <- minmax.epsilon + (((scores - min(scores)) * (1 - 2 * minmax.epsilon)) / (max(scores) - min(scores)))
      return(normed_scores)
    })
    program_scores <- purrr::reduce(program_scores, cbind)
  }
  # reformat program scores depending on input format
  if (inherits(expr.mat, "matrix") || inherits(expr.mat, "array") || inherits(expr.mat, "dgCMatrix")) {
    colnames(program_scores) <- program.labels
  } else {
    for (g in seq(ncol(program_scores))) {
      if (inherits(expr.mat, "SingleCellExperiment")) {
        SummarizedExperiment::colData(expr.mat)[, program.labels[g]] <- program_scores[, g]
      } else if (inherits(expr.mat, "Seurat")) {
        expr.mat <- Seurat::AddMetaData(expr.mat,
                                        metadata = program_scores[, g],
                                        program.labels[g])
      }
    }
  }
  # return results
  if (inherits(expr.mat, "matrix") || inherits(expr.mat, "array") || inherits(expr.mat, "dgCMatrix")) {
    return(program_scores)
  } else {
    return(expr.mat)
  }
}
