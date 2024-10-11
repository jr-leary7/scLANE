#' Perform GSEA on dynamic genes identified by \code{scLANE}.
#'
#' @name enrichDynamicGenes
#' @author Jack R. Leary
#' @description This function uses the \code{gprofiler2} package to perform pathway analysis on a set of genes from one or more lineages that were determined to be dynamic with \code{\link{testDynamic}}.
#' @import magrittr
#' @importFrom dplyr filter arrange desc distinct pull
#' @param scLANE.de.res The output from \code{\link{getResultsDE}}. Defaults to NULL.
#' @param lineage A character vector specifying lineages to isolate. Defaults to NULL.
#' @param species The species against which to run enrichment analysis. Defaults to "hsapiens".
#' @return The output from \code{\link[gprofiler2]{gost}}.
#' @seealso \code{\link[gprofiler2]{gost}}
#' @export
#' @examples
#' data(scLANE_models)
#' scLANE_de_res <- getResultsDE(scLANE_models)
#' enr_res <- enrichDynamicGenes(scLANE_de_res)

enrichDynamicGenes <- function(scLANE.de.res = NULL,
                               lineage = NULL,
                               species = "hsapiens") {
  # check inputs
  if (is.null(scLANE.de.res)) { stop("Arguments to enrichDynamicGenes() are missing.") }
  species <- tolower(species)
  if (!is.null(lineage)) {
    scLANE.de.res <- dplyr::filter(scLANE.de.res,
                                   Lineage %in% lineage,
                                   Gene_Dynamic_Lineage == 1)
  } else {
    scLANE.de.res <- dplyr::filter(scLANE.de.res, Gene_Dynamic_Overall == 1)
  }
  genes <- dplyr::arrange(scLANE.de.res, dplyr::desc(Test_Stat)) %>%
           dplyr::distinct(Gene) %>%
           dplyr::pull(Gene)
  pathway_enr_res <- gprofiler2::gost(query = genes,
                                      organism = species,
                                      ordered_query = TRUE,
                                      multi_query = FALSE,
                                      significant = FALSE)
  return(pathway_enr_res)
}
