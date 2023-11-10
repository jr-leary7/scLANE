#' Add per-cell module scores for gene programs.
#'
#' @name geneProgramScoring
#' @author Jack Leary
#' @import magrittr
#' @importFrom Matrix Matrix
#' @importFrom ggplot2 ggplot aes geom_point geom_vline geom_ribbon geom_line scale_x_continuous labs
#' @importFrom scales label_number
#' @description This function uses \code{\link[UCell]{ScoreSignatures_UCell}} to create a per-cell module score for each of the provided gene clusters. If the
#' @param expr.mat Either a \code{SingleCellExperiment} or \code{Seurat} object from which counts can be extracted, or a matrix of integer-valued counts with genes as rows & cells as columns. Defaults to NULL.
#' @param genes A character vector of gene IDs. Defaults to NULL.
#' @param gene.clusters A factor containing the cluster assignment of each gene in \code{genes}. Defaults to NULL.
#' @param program.labels (Optional) A character vector specifying a label for each gene cluster. Defaults to NULL.
#' @param n.cores (Optional) The number of cores used under the hood in \code{\link[UCell]{ScoreSignatures_UCell}}. Defaults to 2.
#' @return Either a \code{Seurat} or \code{SingleCellExperiment} object if \code{expr.mat} is in either form, or a data.frame containing per-cell program scores if \code{expr.mat} is a matrix.
#' @export
#' @examples
#' \dontrun{
#' geneProgramScoring(seu_obj,
#'                    genes = gene_embed$gene,
#'                    gene.clusters = gene_embed$leiden,
#'                    program.labels = c("cell cycle", "organogenesis"))
#' }

geneProgramScoring <- function(expr.mat = NULL,
                               genes = NULL,
                               gene.clusters = NULL,
                               program.labels = NULL,
                               n.cores = 2) {
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
