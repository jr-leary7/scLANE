#' Plot gene dynamics with estimated coefficients.
#'
#' @name plotModelCoefs
#' @author Jack R. Leary
#' @import magrittr
#' @importFrom dplyr select mutate lag lead
#' @importFrom tidyr pivot_longer
#' @importFrom ggplot2 ggplot aes geom_point geom_vline geom_ribbon geom_line 
#' @importFrom ggplot2 scale_x_continuous labs
#' @importFrom scales label_number
#' @description Generate a plot of gene dynamics over a single pseudotime 
#' lineage, along with a table of coefficients across pseudotime intervals.
#' @param test.dyn.res The output from \code{\link{testDynamic}}. 
#' Defaults to NULL.
#' @param gene A character specifying which gene's dynamics should be plotted. 
#' Defaults to NULL.
#' @param pt A data.frame of pseudotime values for each cell. Defaults to NULL.
#' @param expr.mat Either a \code{SingleCellExperiment} or \code{Seurat} object 
#' from which counts can be extracted, or a matrix of integer-valued counts 
#' with genes as rows & cells as columns. Defaults to NULL.
#' @param size.factor.offset (Optional) An offset to be used to rescale the 
#' fitted values. Can be generated easily with \code{\link{createCellOffset}}. 
#' No need to provide if the GEE backend was used. Defaults to NULL.
#' @param log1p.norm (Optional) Should log1p-normalized versions of expression 
#' & model predictions be returned as well? Defaults to TRUE.
#' @param lineage A character vector specifying which lineage should be 
#' plotted. Should be letters, i.e. lineage "A" or "B". Defaults to "A".
#' @return A \code{ggplot2} object displaying a gene dynamics plot & a table 
#' of coefficients across pseudotime intervals.
#' @export
#' @examples
#' data(sim_counts)
#' data(sim_pseudotime)
#' data(scLANE_models)
#' cell_offset <- createCellOffset(sim_counts)
#' scLANE_de_res <- getResultsDE(scLANE_models)
#' plotModelCoefs(scLANE_models,
#'     gene = "ACLY",
#'     pt = sim_pseudotime,
#'     expr.mat = sim_counts,
#'     size.factor.offset = cell_offset
#' )
plotModelCoefs <- function(test.dyn.res = NULL,
    gene = NULL,
    pt = NULL,
    expr.mat = NULL,
    size.factor.offset = NULL,
    lineage = "A",
    log1p.norm = TRUE) {
    # check inputs
    if (is.null(test.dyn.res) || is.null(gene) ||
        is.null(pt) || is.null(expr.mat)) {
        stop("Arguments to plotModelCoefs() are missing.")
    }
    # pull fitted values
    all_lineages <- gsub("Lineage_", "", names(test.dyn.res[[1]]))
    if (length(all_lineages) == 1) {
        gfv_filter <- NULL
    } else {
        gfv_filter <- all_lineages[all_lineages != lineage]
    }
    fitted_vals <- getFittedValues(test.dyn.res,
        genes = gene,
        pt = pt,
        expr.mat = expr.mat,
        size.factor.offset = size.factor.offset,
        log1p.norm = log1p.norm,
        filter.lineage = gfv_filter
    )
    if (log1p.norm) {
        fitted_vals <- dplyr::select(fitted_vals,
            cell,
            lineage,
            pt,
            gene,
            rna = rna_log1p,
            scLANE_pred = scLANE_pred_log1p,
            scLANE_ci_ll = scLANE_ci_ll_log1p,
            scLANE_ci_ul = scLANE_ci_ul_log1p
        )
    } else {
        fitted_vals <- dplyr::select(
            fitted_vals,
            cell,
            lineage,
            pt,
            gene,
            rna,
            scLANE_pred,
            scLANE_ci_ll,
            scLANE_ci_ul
        )
    }
    # generate dynamics plot
    dyn_plot <- ggplot2::ggplot(fitted_vals, ggplot2::aes(x = pt, y = rna)) +
        ggplot2::geom_point(
            size = 1.5,
            alpha = 0.6,
            stroke = 0,
            color = "grey30"
        ) +
        ggplot2::geom_vline(
            data = data.frame(
                gene = gene,
                knot = unique(test.dyn.res[[gene]][[paste0(
                    "Lineage_",
                    lineage
                )]]$MARGE_Slope_Data$Breakpoint)
            ),
            mapping = ggplot2::aes(xintercept = knot),
            linetype = "dashed",
            color = "black",
            linewidth = 0.75
        ) +
        ggplot2::geom_ribbon(
            ggplot2::aes(
                ymin = scLANE_ci_ll,
                ymax = scLANE_ci_ul
            ),
            linewidth = 0,
            fill = "darkgreen",
            alpha = 0.35
        ) +
        ggplot2::geom_line(ggplot2::aes(y = scLANE_pred),
            color = "darkgreen",
            linewidth = 0.75
        ) +
        ggplot2::scale_x_continuous(
            labels =
                scales::label_number(accuracy = 0.01)
        ) +
        ggplot2::labs(
            x = "Pseudotime",
            y = ifelse(log1p.norm, "Normalized Expression", "Expression")
        ) +
        theme_scLANE()
    # generate coefficient summary
    min_pt <- min(pt[, which(LETTERS == lineage)], na.rm = TRUE)
    max_pt <- max(pt[, which(LETTERS == lineage)], na.rm = TRUE)
    n_breaks <- sum(grepl(
        "Breakpoint",
        colnames(test.dyn.res[[gene]][[paste0(
            "Lineage_",
            lineage
        )]]$Gene_Dynamics)
    ))
    if (n_breaks == 1) {
        coef_sumy <- dplyr::select(
            test.dyn.res[[gene]][[paste0(
                "Lineage_",
                lineage
            )]]$Gene_Dynamics,
            -dplyr::starts_with("Trend")
        ) %>%
            tidyr::pivot_longer(dplyr::starts_with("Slope"),
                names_to = "Segment",
                values_to = "Coef"
            ) %>%
            dplyr::mutate(
                Breakpoint_Lag = dplyr::lag(Breakpoint),
                Breakpoint_Lead = dplyr::lead(Breakpoint),
                Interval = NA_character_,
                .before = 4
            ) %>%
            dplyr::mutate(
                Breakpoint_Lag = dplyr::if_else(is.na(Breakpoint_Lag), min_pt,
                    Breakpoint_Lag
                ),
                Breakpoint_Lead = dplyr::if_else(is.na(Breakpoint_Lead),
                    max_pt, Breakpoint_Lead
                ),
                Interval = paste0(
                    "(", round(Breakpoint_Lag, 3), ", ",
                    round(Breakpoint_Lead, 3), ")"
                )
            ) %>%
            dplyr::select(Interval, Coef) %>%
            dplyr::mutate(Coef = round(Coef, 3))
    } else {
        stop("Number of knots > 1 not currently supported.")
    }
    # convert coefficient summary to grob
    coef_sumy_grob <- gridExtra::tableGrob(coef_sumy,
        cols = c("Interval", "Slope"),
        theme = gridExtra::ttheme_minimal(
            base_size = 11,
            core = list(fg_params = list(hjust = 0, x = 0.05)),
            colhead = list(fg_params = list(hjust = 0, x = 0.05))
        ),
        rows = NULL
    ) %>%
        gtable::gtable_add_grob(
            grobs = grid::rectGrob(gp = grid::gpar(fill = NA, lwd = 3)),
            t = 1,
            b = nrow(.),
            l = 1,
            r = ncol(.)
        ) %>%
        gtable::gtable_add_grob(
            grobs = grid::rectGrob(gp = grid::gpar(fill = NA, lwd = 3)),
            t = 1,
            l = 1,
            r = ncol(.)
        )

    # combine objects
    dyn_plot_anno <- ggpubr::ggarrange(dyn_plot,
        coef_sumy_grob,
        ncol = 2,
        nrow = 1,
        widths = c(2, 1)
    )
    return(dyn_plot_anno)
}
