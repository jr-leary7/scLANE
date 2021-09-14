#' Plot results of \code{marge} models using \code{ggplot2}.
#'
#' @name PlotMARGE

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
  # plot
  p <- ggplot2::ggplot(plot_df, aes(x = PT, y = Counts)) +
       ggplot2::geom_point(alpha = 0.75, size = 2.5) +
       ggplot2::geom_line(aes(x = PT, y = Predicted, color = "MARGE"), size = 2) +
       ggplot2::scale_y_continuous(labels = scales::comma_format()) +
       ggplot2::labs(x = "Pseudotime", y = "Expression", title = gene, color = "Model") +
       ggplot2::theme_classic() +
       ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0.5))
  if (!is.null(null.mod)) {
    p <- p +
        ggplot2::geom_line(aes(x = PT, y = predict(null.mod, type = "response"), color = "Null"), size = 2) +
        ggplot2::scale_color_manual(values = c("dodgerblue", "forestgreen"))
  } else {
    p <- p + ggplot2::scale_color_manual(values = "dodgerblue")
  }
  return(p)
}
