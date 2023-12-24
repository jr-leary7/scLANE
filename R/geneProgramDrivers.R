#' Identify driver genes for a given gene program.
#'
#' @name geneProgramDrivers
#' @author Jack Leary
#' @importFrom Matrix Matrix
#' @importFrom purrr map reduce
#' @importFrom foreach foreach registerDoSEQ
#' @importFrom parallel makeCluster clusterEvalQ stopCluster
#' @importFrom withr with_output_sink
#' @importFrom utils txtProgressBar setTxtProgressBar
#' @importFrom doSNOW registerDoSNOW
#' @importFrom stats cor.test p.adjust
#' @importFrom dplyr arrange desc mutate filter
#' @description This function computes the correlation
#' @param expr.mat Either a \code{SingleCellExperiment} or \code{Seurat} object from which counts can be extracted, or a matrix of normalized counts with genes as rows & cells as columns. Defaults to NULL.
#' @param genes A character vector of genes to test. Defaults to NULL.
#' @param gene.program A vector of program scores as returned by \code{\link{geneProgramScoring}}. Defaults to NULL.
#' @param cor.method (Optional) The correlation method to be used. Defaults to "spearman".
#' @param fdr.cutoff (Optional) The FDR threshold for determining statistical significance. Defaults to 0.01.
#' @param p.adj.method (Optional) The method used to adjust \emph{p}-values for multiple hypothesis testing. Defaults to "holm".
#' @param n.cores (Optional) The number of cores used when iterating over genes to perform testing. Defaults to 2.
#' @param verbose (Optional) Should a progress bar be printed to the console during processing? Defaults to TRUE.
#' @return Either a \code{Seurat} or \code{SingleCellExperiment} object if \code{expr.mat} is in either form, or a data.frame containing per-cell program scores if \code{expr.mat} is a matrix.
#' @seealso \code{\link{geneProgramScoring}}
#' @seealso \code{\link[stats]{cor.test}}
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
#' program_drivers <- geneProgramDrivers(sim_counts,
#'                                       genes = gene_embed$gene,
#'                                       gene.program = sim_counts$cluster_0,
#'                                       fdr.cutoff = 0.05)

geneProgramDrivers <- function(expr.mat = NULL,
                               genes = NULL,
                               gene.program = NULL,
                               cor.method = "spearman",
                               fdr.cutoff = 0.01,
                               p.adj.method = "holm",
                               n.cores = 2L,
                               verbose = TRUE) {
  # check inputs
  if (is.null(expr.mat) || is.null(genes) || is.null(gene.program)) { stop("Arguments to geneProgramDrivers() are missing.") }
  cor.method <- tolower(cor.method)
  if (!cor.method %in% c("pearson", "spearman", "kendall")) { stop("Please specify a valid correlation metric.") }
  # set up counts matrix
  if (inherits(expr.mat, "SingleCellExperiment")) {
    counts_matrix <- SingleCellExperiment::logcounts(expr.mat)
  } else if (inherits(expr.mat, "Seurat")) {
    counts_matrix <- Seurat::GetAssayData(expr.mat,
                                          slot = "data",
                                          assay = Seurat::DefaultAssay(expr.mat))
  } else if (inherits(expr.mat, "dgCMatrix")) {
    counts_matrix <- Matrix::Matrix(expr.mat, sparse = FALSE)
  }
  # set up parallel processing
  if (n.cores > 1L) {
    cl <- parallel::makeCluster(n.cores)
    doSNOW::registerDoSNOW(cl)
  } else {
    cl <- foreach::registerDoSEQ()
  }
  # set up progress bar
  if (verbose) {
    withr::with_output_sink(tempfile(), {
      pb <- utils::txtProgressBar(0, length(genes), style = 3)
    })
    progress_fun <- function(n) utils::setTxtProgressBar(pb, n)
    snow_opts <- list(progress = progress_fun)
  } else {
    snow_opts <- list()
  }
  # iteratively compute correlations & p-values
  cor_tests <- foreach::foreach(g = seq(genes),
                                .combine = "list",
                                .multicombine = ifelse(length(genes) > 1, TRUE, FALSE),
                                .maxcombine = ifelse(length(genes) > 1, length(genes), 2),
                                .packages = c("stats"),
                                .errorhandling = "pass",
                                .inorder = TRUE,
                                .verbose = FALSE,
                                .options.snow = snow_opts) %dopar% {
    cor_res <- stats::cor.test(counts_matrix[genes[g], ],
                               gene.program,
                               method = cor.method,
                               exact = FALSE)
    cor_df <- data.frame(gene = genes[g],
                         corr = unname(cor_res$estimate),
                         pvalue = cor_res$p.value)
    cor_df
  }
  # end parallelization & clean up each worker node
  withr::with_output_sink(tempfile(), {
    if (n.cores > 1L) {
      parallel::clusterEvalQ(cl, expr = {
        gc(verbose = FALSE, full = TRUE)
      })
      parallel::stopCluster(cl)
    }
  })
  cor_tests <- purrr::reduce(cor_tests, rbind)
  rownames(cor_tests) <- genes
  cor_tests <- dplyr::arrange(cor_tests,
                              pvalue,
                              dplyr::desc(abs(corr))) %>%
               dplyr::mutate(pvalue_adj = stats::p.adjust(pvalue, method = p.adj.method)) %>%
               dplyr::filter(pvalue_adj < fdr.cutoff)
  return(cor_tests)
}
