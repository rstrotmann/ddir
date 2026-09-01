#' inner constructor for induction_experiment objects
#'
#' @param data A tibble.
#' @param precipitant The precipitant as character.
#'
#' @returns A induction_experiment object.
#' @noRd
new_induction_experiment <- function(data, precipitant = "") {
  structure(
    as_tibble(data),
    class = unique(c("induction_experiment", "tbl_df", "tbl", "data.frame")),
    precipitant = precipitant
  )
}


#' constructor for induction_experiment objects
#'
#' @param data A tibble.
#' @param precipitant The precipitant as character.
#'
#' @returns An induction_experiment object.
#' @export
induction_experiment <- function(data = NULL, precipitant = "") {
  if (is.null(data))
    data <- tibble(
      DONOR = character(),
      SAMPLE = character(),
      CONC = numeric(),
      OBJECT = character(),
      FOLD = numeric(),
      REL = numeric(),
      SOURCE = character()
    )
  out <- new_induction_experiment(data, precipitant)
  validate_induction_experiment(out)
}


#' Validate induction_experiment object
#'
#' @param obj A induction_experiment object.
#' @param allow_null Logical.
#' @param allowed_objects Character.
#'
#' @returns The induction_experiment object with unknown objects deleted.
#' @noRd
validate_induction_experiment <- function(
    obj, allow_null = FALSE, allowed_objects = NULL
  ) {
  # early fail on NULL, if applicable
  if (isTRUE(allow_null)) {
    if (is.null(obj))
      return(obj)
  }

  obj_name <- deparse(substitute(obj))

  if (!inherits(obj, "induction_experiment"))
    stop(paste(obj_name, "must be an induction_experiment object"))

  # check required fields
  expected_fields <- c("DONOR", "CONC", "OBJECT", "SOURCE")
  missing_fields <- setdiff(expected_fields, names(obj))
  if (length(missing_fields) > 0) {
    stop(paste0(
      "Missing columns in ", obj_name, ": ",
      nice_enumeration(missing_fields)
    ))
  }
  if (!any(c("FOLD", "REL") %in% names(obj))) {
    stop(paste("Either FOLD or REL must be in", obj_name))
  }

  # filter for allowed objects
  if (is.null(allowed_objects)) {
    allowed_objects <- unique(obj$OBJECT)
  } else {
    unexpected_objects <- setdiff(unique(obj$OBJECT), allowed_objects)
    if (length(unexpected_objects) > 0) {
      warning(paste0(
        "Unexpected objects in ", obj_name, " removed: ",
        nice_enumeration(unexpected_objects)
      ))
    }
  }

  filter(obj, .data$OBJECT %in% allowed_objects)
}


#' Print generic for induction_experiment
#'
#' @param x The induction_experiment object.
#'
#' @param ... Further parameters.
#'
#' @exportS3Method base::print
print.induction_experiment <- function(x, ...) {
  precipitant <- attr(x, "precipitant")

  if (isTRUE(getOption("knitr.in.progress"))) {
    x |>
      mutate(SOURCE = case_when(
        !is.na(SOURCE) & !SOURCE == "" ~ SOURCE,
        .default = ""
      )) |>
      knitr::kable(caption = paste(
        "CYP induction experiment for precipitant", precipitant
      ))
  } else {
    cat(paste(hline(), "CYP induction experiment", hline(), "\n"))
    cat(paste(
      "Experimental data for in vitro CYP induction by",
      precipitant, "\n\n"))

    out <- x |>
      mutate(SOURCE = case_when(
        !is.na(SOURCE) & SOURCE != "" ~ SOURCE,
        .default = ""
      ))

    cat(df_to_string(out, colnames = TRUE))
  }
}


#' Plot generic for induction_experiment
#'
#' @param x The induction_experiment object.
#'
#' @param type Columns to plot on y axis.
#' @param ... Further parameters.
#'
#' @exportS3Method graphics::plot
plot.induction_experiment <- function(x, type = "FOLD", ...) {
  induction_plot(x, type = type, ...) +
    labs(x = paste(attr(x, "precipitant"), "uM"),
         title = paste("CYP induction by", attr(x, "precipitant")))
}
