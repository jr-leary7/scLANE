#' Pull the set of knots for dynamic genes across each lineage.
#'
#' @name getKnotDist
#' @author Jack R. Leary
#' @import magrittr
#' @importFrom purrr imap reduce
#' @importFrom dplyr filter
#' @description Pulls knot locations for dynamic genes across each lineage, 
#' allowing comparisons of where transcriptional switches occur 
#' between lineages.
#' @param test.dyn.res The output from \code{\link{testDynamic}}. Defaults to 
#' NULL.
#' @param dyn.genes The set of genes to pull knots for. If unspecified, pulls 
#' knots for all modeled genes. Defaults to NULL.
#' @return A data.frame containing gene name, lineage ID, and knot location in 
#' pseudotime.
#' @export
#' @examples
#' data(scLANE_models)
#' knot_dist <- getKnotDist(scLANE_models)
getKnotDist <- function(test.dyn.res = NULL, dyn.genes = NULL) {
    # check inputs
    if (is.null(test.dyn.res)) {
        stop("You forgot one of the arguments to getKnotDist().")
    }
    if (is.null(dyn.genes)) {
        dyn.genes <- names(test.dyn.res)
    }
    # pull knots per-lineage
    knot_df <- purrr::imap(test.dyn.res[dyn.genes], \(x, y) {
        purrr::imap(x, \(z, w) {
            if (is.null(z$MARGE_Slope_Data$Breakpoint)) {
                knots <- NA_real_
            } else {
                knots <- z$MARGE_Slope_Data$Breakpoint
            }
            data.frame(
                gene = y,
                lineage = w,
                knot = knots
            )
        }) %>%
            purrr::reduce(rbind)
    }) %>%
        purrr::reduce(rbind) %>%
        dplyr::filter(!is.na(knot))
    return(knot_df)
}
