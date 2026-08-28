
#' inhibition_data constructor
#'
#' @param data A data frame with the columns object, ki and source.
#' @param precipitant The precipitant name as character.
#'
#' @returns An inhibition_data object.
#' @export
inhibition_data <- function(data = NULL, precipitant = "") {
  if (is.null(data)) {
    data <- data.frame(
      object = character(),
      ki = numeric(),
      source = character()
    )
  }
  out <- data
  class(out) <- c("inhibition_data", "data.frame")
  attr(out, "precipitant") <- precipitant

  validate_inhibition_data(out)
  out
}


#' Reconstruct inhibition_data after dplyr verbs
#'
#' @param data The data as data frame.
#' @param template The object template.
#'
#' @returns A inhibition_data object.
#' @export
#' @noRd
dplyr_reconstruct.inhibition_data <- function(data, template) {
  inhibition_data(data, attr(template, "precipitant"))
}


#' Print method for inhibition_data objects
#'
#' @param x The inhibition
#' @param ... Further arguments.
#'
#' @returns Nothing.
#' @import dplyr
#' @import knitr
#' @export
print.inhibition_data <- function(x, ...) {
  validate_inhibition_data(x)

  if (isTRUE(getOption("knitr.in.progress"))) {
   x |>
      mutate(source = case_when(
        !is.na(source) & !source == "" ~ source,
        .default = ""
      )) |>
      knitr::kable(caption = paste(
        "In vitro inhibition data for precipitant",
        attr(x, "precipitant")
      ))
  } else {
    cat(paste0(hline(), " DDI inhibition data ", hline(), "\n"))

    cat(paste(
      "In vitro inhibition data for precipitant", attr(x, "precipitant"), "\n\n"))
    out <- x |>
      mutate(source = case_when(
        !is.na(source) & !source == "" ~ source,
        .default = ""
      ))

    cat(df_to_string(out, colnames = TRUE))
  }
}
