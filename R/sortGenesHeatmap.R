#' Sort genes by where their peak expression occurs across pseudotime.
#'
#' @name sortGenesHeatmap
#' @author Jack Leary
#' @import magrittr
#' @importFrom purrr map reduce
#' @importFrom dplyr filter distinct select with_groups summarise arrange pull
#' @description Sort genes such that genes with peak expression occurring earlier in pseudotime are first, and vice versa for genes with late peak expression. Useful for ordering genes in order to create heatmaps of expression cascades.
#' @param heatmap.mat A matrix of raw or smoothed expression values with genes as columns and cells as rows. Defaults to NULL.
#' @param pt.vec A numeric vector of pseudotime values for each cell i.e., for each row in the heatmap matrix. Defaults to NULL.
#' @return A character vector of genes sorted by their peak expression values over pseudotime.
#' @seealso \code{\link{smoothedCountsMatrix}}
#' @export
#' @examples
#' \dontrun{
#' smoothed_counts <- smoothedCountsMatrix(gene_stats, pt = pt_df)
#' sortGenesheatmap(heatmap.mat = smoothed_counts$Lineage_A,
#'                  pt.vec = pt_df[!is.na(pt_df$Lineage_A), ]$Lineage_A)
#' }

sortGenesHeatmap <- function(heatmap.mat = NULL, pt.vec = NULL) {
  # check inputs
  if (!inherits(heatmap.mat, "matrix")) {
    heatmap.mat <- try(as.matrix(heatmap.mat), silent = TRUE)
    if (inherits(heatmap.mat, "try-error")) {
      stop("heatmap.mat must be coerceable to a matrix.")
    }
  }
  if (!is.numeric(pt.vec) || any(is.na(pt.vec))) { stop("pt.vec must be a numeric vector with no NA values.") }

  # identify point at which peak expression occurs for each gene across pseudotime
  gene_peak_order <- purrr::map(seq_len(ncol(heatmap.mat)), \(x) {
    data.frame(gene = colnames(heatmap.mat)[x],
               pseudotime = pt.vec,
               mRNA = heatmap.mat[, x]) %>%
      dplyr::filter(mRNA == max(mRNA)) %>%
      dplyr::distinct() %>%
      dplyr::select(gene,
                    pseudotime,
                    mRNA)
  }) %>%
    purrr::reduce(rbind) %>%
    dplyr::with_groups(gene,
                       dplyr::summarise,
                       mu = mean(pseudotime)) %>%
    dplyr::arrange(mu) %>%
    dplyr::pull(gene)
  return(gene_peak_order)
}
