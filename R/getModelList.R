#' Retrieve models from \code{testDynamic()} output.
#'
#' @name getModelList
#' @author Jack Leary
#' @description This function retrieves a list of either \code{marge} or null models from the output of \code{testDynamic()}.
#' @param test.results The nested-list output of \code{testDynamic()}. Defaults to NULL.
#' @param which.model Whether the \code{marge} or null models should be returned. Defaults to "marge".
#' @return A list of \code{marge} or null model NB GLM objects.
#' @export
#' @examples
#' \dontrun{getModelList(test.results = gene_test_list)}
#' \dontrun{getModelList(test.results = gene_test_list, which.model = "null")}

getModelList <- function(test.results = NULL, which.model = "marge") {
  # check inputs
  if (is.null(test.results)) { stop("You forgot to provide a list of models to getModelList().") }
  which.model <- tolower(which.model)
  if (!which.model %in% c("marge", "null")) { stop("which.model must be either 'marge' or 'null'.") }
  # get models
  if (which.model == "marge") {
    res <- lapply(test.results, function(x) x[9][[1]][[1]])
  } else {
    res <- lapply(test.results, function(x) x[10][[1]])
  }
  return(res)
}
