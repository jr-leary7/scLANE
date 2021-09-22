#' Plot results of \code{marge} models using \code{ggplot2}.
#'
#' @name PlotMARGE
#' @description This function visualizes the fitted values of a \code{marge} model over the expression and pseudotime values of each cell.
#' @param model The \code{marge} model to extract fitted values from. Defaults to NULL.
#' @param gene.counts A vector of integer expression values for each cell. Defaults to NULL.
#' @param pt A data.frame of pseudotime values for each cell. Defaults to NULL.
#' @param gene The name of the gene that's being analyzed. Used as the title of the \code{ggplot} object. Defaults to NULL.
#' @param null.mod (Optional) An \code{lm} or \code{glm} object specifying the intercept-only null model for comparison. Defaults to NULL.
#' @return A \code{ggplot} object.
#' @export
#' @examples
#' PlotMARGE(model = marge_mod, gene.counts = exp_vec, pt = pt_df, gene = "BRCA2", null.mod = lm(exp_vec ~ 1))

PlotMARGE <- function(model = NULL,
                      gene.counts = NULL,
                      pt = NULL,
                      gene = NULL,
                      null.mod = NULL) {
  # check inputs
  if (any(sapply(c(model, gene.counts, pt, gene), is.null))) stop("You forgot one or more of the arguments.")
  # create dataframe
  plot_df <- data.frame(PT = pt[, 1],
                        Predicted = predict(model, type = "response"),
                        Counts = gene.counts)
  # generate base plot
  p <- ggplot2::ggplot(plot_df, aes(x = PT, y = Counts)) +
       ggplot2::geom_point(alpha = 0.75, size = 2.5) +
       ggplot2::geom_line(aes(x = PT, y = Predicted, color = "MARGE"), size = 2) +
       ggplot2::scale_y_continuous(labels = scales::comma_format()) +
       ggplot2::labs(x = "Pseudotime", y = "Expression", title = gene, color = "Model") +
       ggplot2::theme_classic() +
       ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0.5))
  # optionally add line of null model fitted values
  if (!is.null(null.mod)) {
    p <- p +
        ggplot2::geom_line(aes(x = PT, y = fitted.values(null.mod), color = "Null"), size = 2) +
        ggplot2::scale_color_manual(values = c("dodgerblue", "forestgreen"))
  } else {
    p <- p + ggplot2::scale_color_manual(values = "dodgerblue")
  }
  return(p)
}
