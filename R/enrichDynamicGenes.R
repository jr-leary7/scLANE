#' Perform GSEA on dynamic genes identified by \code{scLANE}.
#'
#' @name enrichDynamicGenes
#' @author Jack Leary
#' @description This function uses the \code{msigdbr} and \code{clusterProfiler} packages to perform GSEA on a set of genes from one or more lineages that were determined to be dynamic with \code{\link{testDynamic}}.
#' @import magrittr
#' @importFrom dplyr filter distinct pull
#' @param scLANE.de.res The output from \code{\link{getResultsDE}}. Defaults to NULL.
#' @param lineage A character vector specifying lineages to isolate. Defaults to NULL.
#' @param species One of "hs" or "mm", specifying whether the gene sets are human or murine in origin.
#' @param gene.set.cat Corresponds to the \code{category} parameter of \code{\link[msigdbr]{msigdbr}}. Defaults to NULL.
#' @param gene.set.subcat Corresponds to the \code{subcategory} parameter of \code{\link[msigdbr]{msigdbr}}. Defaults to NULL.
#' @return The output from \code{\link[clusterProfiler]{enricher}}.
#' @seealso \code{\link[msigdbr]{msigdbr}}
#' @seealso \code{\link[msigdbr]{msigdbr_collections}}
#' @seealso \code{\link[clusterProfiler]{enricher}}
#' @export
#' @examples
#' \dontrun{
#' enrichDynamicGenes(scLANE.de.res = de_stats,
#'                    lineage = "A",
#'                    gene.set.cat = "C1",
#'                    gene.set.subcat = "CGP")
#' enrichDynamicGenes(scLANE.de.res = de_stats,
#'                    species = "mm",
#'                    gene.set.cat = "C2",
#'                    gene.set.subcat = "CP:REACTOME")
#' }

enrichDynamicGenes <- function(scLANE.de.res = NULL,
                               lineage = NULL,
                               species = "hs",
                               gene.set.cat = NULL,
                               gene.set.subcat = NULL) {
  # check inputs
  if (is.null(scLANE.de.res) || is.null(gene.set.cat)) { stop("Arguments to enrichDynamicGenes() are missing.") }
  species <- tolower(species)
  if (!species %in% c("hs", "mm")) { stop("species must be one of 'hs' or 'mm' at this time.") }
  if (!is.null(lineage)) {
    scLANE.de.res <- dplyr::filter(scLANE.de.res,
                                   Lineage %in% lineage,
                                   Gene_Dynamic_Lineage == 1)
  } else {
    scLANE.de.res <- dplyr::filter(scLANE.de.res, Gene_Dynamic_Overall == 1)
  }
  genes <- scLANE.de.res %>%
           dplyr::distinct(Gene) %>%
           dplyr::pull(Gene)
  if (species == "hs") {
    gene_sets <- msigdbr::msigdbr(species = "human",
                                  category = gene.set.cat,
                                  subcategory = gene.set.subcat)
  } else if (species == "mm") {
    gene_sets <- msigdbr::msigdbr(species = "mouse",
                                  category = gene.set.cat,
                                  subcategory = gene.set.subcat)
  }
  term_to_gene <- gene_sets %>% dplyr::distinct(gs_name, gene_symbol)
  gsea_res <- clusterProfiler::enricher(gene = genes,
                                        TERM2GENE = term_to_gene,
                                        universe = unique(term_to_gene$gene_symbol),
                                        qvalueCutoff = 0.05)
  return(gsea_res)
}
