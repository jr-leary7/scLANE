#' Represent a \code{marge} model as a series of piecewise equations.
#'
#' @name summarizeModel
#' @author Jack Leary
#' @import magrittr
#' @importFrom dplyr mutate case_when arrange lead lag rowwise ungroup select
#' @description This function converts a \code{marge} model from using truncated power basis functions to a more interpretable series of piecewise linear equations over a set of empirically chosen cutpoint intervals. Note that the function assumes that pseudotime has been normalized to exist over \eqn{[0, 1]}.
#' @param marge.model The fitted model output from \code{\link{marge2}} (this function is internally called by \code{\link{testDynamic}}). Defaults to NULL.
#' @return A data.frame of the model coefficients, cutpoint intervals, and formatted equations.
#' @examples
#' \dontrun{
#' summarizeModel(marge.model = marge_mod)
#' }

summarizeModel <- function(marge.model = NULL) {
  # check inputs
  if (is.null(marge.model)) { stop("Please provide a non-NULL input to summarizeModel().") }
  if (class(marge.model) != "marge") { stop("Input to summarizeModel() must be of class marge.") }

  # extract model equation & slopes
  coef_df <- data.frame(coef_name = names(coef(marge.model$final_mod)),
                        coef_value = unname(coef(marge.model$final_mod))) %>%
             dplyr::mutate(coef_value = round(coef_value, 2),
                           breakpoint = c(0, extractBreakpoints(marge.model$final_mod)$Breakpoint),
                           breakpoint_direction = c(NA_character_, extractBreakpoints(marge.model$final_mod)$Direction),
                           breakpoint_direction_int = dplyr::case_when(is.na(breakpoint_direction) ~ 0L,
                                                                       breakpoint_direction == "Left" ~ 1L,
                                                                       breakpoint_direction == "Right" ~ 2L,
                                                                       TRUE ~ NA_integer_)) %>%
             dplyr::arrange(breakpoint,
                            breakpoint_direction_int) %>%
             dplyr::mutate(next_breakpoint = dplyr::lead(breakpoint, n = 1),
                           prev_breakpoint = dplyr::lag(breakpoint, n = 1),
                           interval_min = dplyr::case_when(coef_name == "B_finalIntercept" ~ 0,
                                                           TRUE ~ prev_breakpoint),
                           interval_max = dplyr::case_when(coef_name == "B_finalIntercept" ~ 1,
                                                           is.na(next_breakpoint) ~ 1,
                                                           TRUE ~ next_breakpoint),
                           overall_intercept = coef_value[1]) %>%
             dplyr::rowwise() %>%
             dplyr::mutate(eqn = dplyr::case_when(coef_name == "B_finalIntercept" ~ paste0(coef_value, " + 0 * x"),
                                                  breakpoint_direction == "Right" ~ paste0(coef_value, " * max(0, x - ", breakpoint, ")"),
                                                  breakpoint_direction == "Left" ~ paste0(coef_value, " * max(0, ", breakpoint, " - x)")),
                           eqn_intercept = dplyr::case_when(coef_name == "B_finalIntercept" ~ as.numeric(strsplit(eqn, " \\+ ")[[1]][1]),
                                                            TRUE ~ NA_real_),
                           eqn_multiplier = dplyr::case_when(breakpoint_direction == "Left" ~ as.numeric(strsplit(eqn, " \\* ")[[1]][1]) * as.numeric(strsplit(gsub("max\\(0, ", "", strsplit(eqn, " \\* ")[[1]][2]), " - ")[[1]][1]),
                                                             breakpoint_direction == "Right" ~ as.numeric(strsplit(eqn, " \\* ")[[1]][1]) * as.numeric(paste0("-", gsub("\\)", "", strsplit(gsub("max\\(0, ", "", strsplit(eqn, " \\* ")[[1]][2]), " - ")[[1]][2]))),
                                                             TRUE ~ NA_real_),
                           eqn_full = dplyr::case_when(coef_name == "B_finalIntercept" ~ as.character(eqn_intercept),
                                                       breakpoint_direction == "Right" ~ paste0((overall_intercept + eqn_multiplier),
                                                                                                " + ",
                                                                                                strsplit(eqn, " \\* ")[[1]][1],
                                                                                                "x"),
                                                       breakpoint_direction == "Left" ~ paste0((overall_intercept + eqn_multiplier),
                                                                                               " + ",
                                                                                               "-",
                                                                                               strsplit(eqn, " \\* ")[[1]][1],
                                                                                               " * x"),
                                                       TRUE ~ NA_character_),
                           eqn_full = dplyr::case_when(breakpoint_direction == "Left" & as.numeric(strsplit(eqn, " \\* ")[[1]][1]) < 0 ~ paste0((overall_intercept + eqn_multiplier),
                                                                                                                                                " + ",
                                                                                                                                                gsub("-", "", strsplit(eqn, " \\* ")[[1]][1]),
                                                                                                                                                "x"),
                                                       TRUE ~ eqn_full),
                           eqn_slope = dplyr::case_when(breakpoint_direction %in% c("Left", "Right") ~ as.numeric(gsub("x", "", strsplit(eqn_full, " \\+ ")[[1]][2])),
                                                        TRUE ~ NA_real_),
                           interval_over = paste0("(", interval_min, ", ", interval_max, ")")) %>%
             dplyr::ungroup() %>%
             dplyr::select(coef_name,
                           coef_value,
                           interval_over,
                           eqn,
                           eqn_full,
                           eqn_slope)
  return(coef_df)
}
