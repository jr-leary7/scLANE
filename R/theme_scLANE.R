#' A \code{ggplot2} theme for \code{scLANE}.
#'
#' @name theme_scLANE
#' @author Jack Leary
#' @importFrom ggplot2 theme_classic theme element_rect
#' @description A publication-ready theme for creating gene dynamics plots, embedding plots, etc.
#' @param base.size The base font size. Defaults to 12.
#' @param base.lwd The base linewidth. Defaults to 0.75.
#' @param base.family The font family to be used throughout. Defaults to "sans".
#' @return A \code{ggplot2} theme.
#' @export
#' @examples
#' \dontrun{
#' plotModels(gene_stats,
#'            gene = "CD14",
#'            pt = pt_df,
#'            expr.mat = count_mat
#' ) + theme_scLANE()
#' }

theme_scLANE <- function(base.size = 12,
                         base.lwd = 0.75,
                         base.family = "sans") {
  ggplot2::theme_classic(base_size = base.size,
                         base_family = base.family,
                         base_line_size = base.lwd,
                         base_rect_size = base.lwd) +
    ggplot2::theme(strip.clip = "off",
                   strip.background = ggplot2::element_rect(linewidth = base.lwd))
}
