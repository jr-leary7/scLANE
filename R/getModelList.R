#' Retrieve models from \code{testDynamic()} output.
#'
#' @name getModelList
#' @description This function retrieves a list of either \code{marge} or null models from the output of \code{testDynamic()}.
#' @param test.results The nested-list output of \code{testDynamic()}. Defaults to NULL.
#' @param which.model Whether the \code{marge} or null models should be returned. Defaults to "marge".
#' @return A list of \code{marge} or null NB GLM objects.
#' @export
#' @examples
#' getModelList(test.results = gene_test_list)
#' getModelList(test.results = gene_test_list, which.model = "null")

getModelList <- function(test.results = NULL, which.model = "marge") {
  # check inputs
  if (is.null(test.results)) { stop("You forgot to provide a list of models to getModelList().") }
  which.model <- tolower(which.model)
  if (!which.model %in% c("marge", "null")) { stop("which.model must be either 'marge' or 'null'.") }
  # get models
  if (which.model == "marge") {
    res <- lapply(gene_stats, function(x) x[[9]]$final_mod)
  } else {
    res <- lapply(gene_stats, function(x) x[[10]])
  }
  return(res)
}
