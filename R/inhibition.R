#' inhibition_data constructor
#'
#' @param data A data frame with the columns object, ki and source.
#' @param precipitant The precipitant name as character.
#'
#' @returns An inhibition_data object.
#' @noRd
new_inhibition_data <- function(data, precipitant = "") {
  structure(
    relocate(data, any_of(c("object", "ic50", "ki", "kinact", "source"))),
    class = unique(c("inhibition_data", "tbl_df", "tbl", "data.frame")),
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
    data <- tibble(
      object = character(), ki = numeric(), source = character()
    )
  }

  out <- new_inhibition_data(as_tibble(data), precipitant)
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
#' @noRd
dplyr_reconstruct.inhibition_data <- function(data, template) {
  new_inhibition_data(data, attr(template, "precipitant"))
}


#' @importFrom vctrs vec_restore
#' @exportS3Method vctrs::vec_restore
#' @noRd
vec_restore.inhibition_data <- function(x, to, ...) {
  new_inhibition_data(x, attr(to, "precipitant"))
}


#' Reconstruct ki from ic50
#'
#' @param data A data frame.
#'
#' @returns A data frame.
#' @noRd
ensure_ki <- function(data) {
  if (!"ki" %in% names(data)) {
    if ("ic50" %in% names(data)) {
      data <- data |>
        mutate(ki = .data$ic50 / 2)

      warning("ki derived from ic50/2 assuming substrate concentration is close to KM")
    } else {
      stop("Input must contain either a ki or an ic50 column!")
    }
  }
  data
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
      "In vitro inhibition data for precipitant", attr(x, "precipitant"),
      "\n\n"))

    out <- x |>
      as.data.frame() |>
      mutate(source = case_when(
        !is.na(.data$source) & !.data$source == "" ~ .data$source,
        .default = ""
      ))

    cat(df_to_string(out, colnames = TRUE))
  }
}
