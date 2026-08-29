#' induction_data constructor
#'
#' @param data A data frame with the columns object, ki and source.
#' @param precipitant The precipitant name as character.
#'
#' @returns An inhibition_data object.
#' @noRd
new_induction_data <- function(data = NULL, precipitant = "") {
  structure(
    data,
    class = unique(c("induction_data", class(data))),
    precipitant = precipitant
  )
}


#' induction_data constructor alias
#'
#' @param data A data frame with the columns object, ki and source.
#' @param precipitant The precipitant name as character.
#'
#' @returns An inhibition_data object.
#' @export
induction_data <- function(data = NULL, precipitant = "") {
  if (is.null(data)) {
    data <- data.frame(
      object = character(),
      emax = numeric(),
      ec50 = numeric(),
      max_c = numeric(),
      source = character()
    )
  }

  out <- new_induction_data(as.data.frame(data), precipitant)

  validate_induction_data(out)
  out
}


#' Reconstruct inhibition_data after dplyr verbs
#'
#' @param data The data as data frame.
#' @param template The object template.
#'
#' @returns A inhibition_data object.
#' @exportS3Method dplyr::dplyr_reconstruct
dplyr_reconstruct.induction_data <- function(data, template) {
  new_induction_data(data, attr(template, "precipitant"))
}


#' @importFrom vctrs vec_restore
#' @exportS3Method vctrs::vec_restore
vec_restore.induction_data <- function(x, to, ...) {
  new_induction_data(x, attr(to, "precipitant"))
}


#' Print method for induction_data objects
#'
#' @param x The inhibition
#' @param ... Further arguments.
#'
#' @returns Nothing.
#' @import dplyr
#' @import knitr
#' @exportS3Method base::print
print.induction_data <- function(x, ...) {
  # validate_induction_data(x)

  if (isTRUE(getOption("knitr.in.progress"))) {
    x |>
      mutate(source = case_when(
        !is.na(source) & !source == "" ~ source,
        .default = ""
      )) |>
      knitr::kable(caption = paste(
        "In vitro CYP induction data for precipitant",
        attr(x, "precipitant")
      ))
  } else {
    cat(paste0(hline(), " DDI induction data ", hline(), "\n"))

    cat(paste(
      "In vitro CYP induction data for precipitant", attr(x, "precipitant"), "\n\n"))
    out <- x |>
      mutate(source = case_when(
        !is.na(source) & !source == "" ~ source,
        .default = ""
      ))

    cat(df_to_string(out, colnames = TRUE))
  }
}
