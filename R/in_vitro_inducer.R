#' In vitro induction data set
#'
#' @slot precipitant character.
#' @slot data data.frame.
#'
#' @returns An induction_experiment object
#' @export
setClass(Class = "induction_experiment",
  representation(
   precipitant = "character",
   data = "data.frame"
  ),
  prototype(
   precipitant = "",
   data = data.frame(
     DONOR = character(0),
     SAMPLE = character(0),
     CONC = character(0),
     OBJECT = character(0),
     FOLD = numeric(0),
     REL = numeric(0),
     SOURCE = character(0)
   )
  )
)


#' induction experiment data set constructor
#'
#' @param data A data frame.
#' @param precipitant The precipitant.
#'
#' @returns induction_experiment object.
#' @export
#'
#' @examples
#' induction_experiment(examplinib_in_vitro_ind, "examplinib")
induction_experiment <- function(data = NULL, precipitant = "") {
  if (is.null(data))
    data = data.frame(
      DONOR = character(0),
      SAMPLE = character(0),
      CONC = character(0),
      OBJECT = character(0),
      FOLD = numeric(0),
      REL = numeric(0),
      SOURCE = character(0)
    )
  new("induction_experiment", data = data, precipitant = precipitant)
}


setValidity("induction_experiment", function(object) {
  validate_argument(object@precipitant, "character", allow_empty = TRUE)

  # validate data
  expected_fields <- c("DONOR", "SAMPLE", "CONC", "OBJECT", "FOLD", "REL", "SOURCE")
  missing_fields <- setdiff(expected_fields, names(object@data))
  if (length(missing_fields) > 0)
    return(paste0("Missing fields: ", nice_enumeration(missing_fields)))

  known_objects <- c(
    "CYP1A1", "CYP1A2", "CYP2A6", "CYP2B6", "CYP2C8", "CYP2C9", "CYP2C18",
    "CYP2C19", "CYP2D6", "CYP2E1", "CYP2J2", "CYP3A4", "CYP3A5", "CYP3A7")

  unknown_objects <- setdiff(
    toupper(unique(object@data$OBJECT)), toupper(known_objects))

  if (length(unknown_objects) > 0)
    warning(paste0("Unknown objects: ", nice_enumeration(unknown_objects)))

  TRUE
})


#' Show method for induction_experiment objects
#'
#' @param object An `induction_experiment` object.
#'
#' @returns Nothing.
#' @export
setMethod(
  "show", "induction_experiment",
  function(object) {
    # line <- paste0(rep("=", 5), collapse="")
    line <- hline()
    cat(paste(line, "CYP induction experiment", line, "\n"))
    cat(paste(
      "Experiental data for in vitro CYP induction by",
      object@precipitant, "\n\n"))

    out <- object@data |>
      mutate(SOURCE = case_when(
        !is.na(SOURCE) & SOURCE != "" ~ SOURCE,
        .default = ""
      ))

    cat(df_to_string(out, colnames = TRUE))
  }
)


#' Pring method for induction_experiment objects
#'
#' @param x An `induction_experiment` object.
#'
#' @returns Nothing.
#' @export
setMethod(
  "print", "induction_experiment",
  function(x) {
    caption <- ifelse(
      x@precipitant != "",
      paste0("Experimental CYP induction data for ", x@precipitant),
      "")

    out <- x@data |>
      mutate(SOURCE = case_when(
        !is.na(SOURCE) & SOURCE != "" ~ SOURCE,
        .default = ""
      )) |>
      knitr::kable(caption = caption, col.names = make_labels(x@data))

    out
  }
)


#' Plot method for induction_experiment objects
#'
#' @param x An `induction_experiment` object.
#' @param type Metric to plot, `"FOLD"` or `"REL"`.
#' @param ... Further arguments.
#'
#' @returns A ggplot object.
#' @export
setMethod(
  "plot", "induction_experiment",
  function(x, type = "FOLD", ...) {
    caption <- ifelse(
      x@precipitant != "",
      paste0("Experimental CYP induction data for ", x@precipitant),
      "")

    induction_plot(x@data, type = type, ...) +
      labs(x = paste(x@precipitant, "uM"),
           title = paste("CYP induction by", x@precipitant))
  }
)
