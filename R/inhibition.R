#' inhibition_data constructor
#'
#' @param data A data frame with the columns object, ki and source.
#' @param precipitant The precipitant name as character.
#'
#' @returns An inhibition_data object.
#' @noRd
new_inhibition_data <- function(data, precipitant = "") {
  structure(
    data,
    class = unique(c("inhibition_data", class(data))),
    precipitant = precipitant
  )
}

#' inhibition_data constructor alias
#'
#' @param data A data frame with the columns object, ki and source.
#' @param precipitant The precipitant name as character.
#'
#' @returns An inhibition_data object.
#' @export
inhibition_data <- function(data = NULL, precipitant = "") {
  if (is.null(data)) {
    data <- data.frame(
      object = character(), ki = numeric(), source = character()
    )
  }
  out <- new_inhibition_data(as.data.frame(data), precipitant)
  validate_inhibition_data(out)
  out
}


#' Reconstruct inhibition_data after dplyr verbs
#'
#' @param data The data as data frame.
#' @param template The object template.
#'
#' @returns A inhibition_data object.
#' @exportS3Method dplyr::dplyr_reconstruct
dplyr_reconstruct.inhibition_data <- function(data, template) {
  new_inhibition_data(data, attr(template, "precipitant"))
}


#' @importFrom vctrs vec_restore
#' @exportS3Method vctrs::vec_restore
vec_restore.inhibition_data <- function(x, to, ...) {
  new_inhibition_data(x, attr(to, "precipitant"))
}


#' Print method for inhibition_data objects
#'
#' @param x The inhibition
#' @param ... Further arguments.
#'
#' @returns Nothing.
#' @import dplyr
#' @import knitr
#' @exportS3Method base::print
print.inhibition_data <- function(x, ...) {
  # validate_inhibition_data(x)

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
