#' Pull the set of knots for dynamic genes across each lineage.
#'
#' @name getKnotDist
#' @author Jack Leary
#' @import magrittr
#' @importFrom purrr imap reduce
#' @description Pulls knot locations for dynamic genes across each lineage, allowing comparisons of where transcriptional switches occur between lineages.
#' @param test.dyn.res The output from \code{\link{testDynamic}}. Defaults to NULL.
#' @param dyn.genes The set of genes to pull knots for. If unspecified, pulls knots for all modeled genes. Defaults to NULL.
#' @return A data.frame containing gene name, lineage ID, and knot location in pseudotime.
#' @export
#' @examples
#' \dontrun{
#' getKnotDist(gene_stats)
#' }

getKnotDist <- function(test.dyn.res = NULL, dyn.genes = NULL) {
  # check inputs
  if (is.null(test.dyn.res)) { stop("You forgot one of the arguments to getKnotDist().") }
  if (is.null(dyn.genes)) {
    dyn.genes <- names(test.dyn.res)
  }
  # pull knots per-lineage
  knot_df <- purrr::imap(test.dyn.res[dyn.genes], \(x, y) {
    purrr::imap(x, \(z, w) {
      data.frame(gene = y,
                 lineage = w,
                 knot = z$MARGE_Slope_Data$Breakpoint)
    }) %>%
      purrr::reduce(rbind)
  }) %>%
    purrr::reduce(rbind)
  return(knot_df)
}