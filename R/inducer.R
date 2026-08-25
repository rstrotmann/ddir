#' DDI CYP inducer data object
#'
#' @slot precipitant The name of the DDI perpetrator as character.
#' @slot data data.frame.
#'
#' @returns Inducer object
#' @export
setClass(Class = "inducer",
  representation(
    precipitant = "character",
    data = "data.frame"
  ),
  prototype(
    precipitant = "",
    data = data.frame(
      object = character(0),
      emax = numeric(0),
      ec50 = numeric(0),
      max_c = numeric(0),
      source = character(0)
    )
  )
)


#' DDI CYP inducer object constructor
#'
#' @param data Data frame with the columns:
#' * object: The CYP induction target. Must be one of "CYP1A1", "CYP1A2",
#'   "CYP2A6", "CYP2B6", "CYP2C8", "CYP2C9", "CYP2C18", "CYP2C19", "CYP2D6",
#'   "CYP2E1", "CYP2J2", "CYP3A4", "CYP3A5", or "CYP3A7".
#' * emax: The maximal induction effect as numeric.
#' * ec50: The EC50, i.e., the concentration causing the half-maximal induction
#'   effect, in uM.
#' * max_c: The maximal concentration tested in the assay in uM.
#' * source: Source information for the parameter.
#'
#' @param precipitant The name of the DDI perpetrator as character.
#'
#' @returns Inducer object
#' @export
#' @examples
#' inducer(data.frame(object = "CYP3A4", emax = 2, ec50 = .1, max_c = 100,
#' source = "test source"), precipitant = "test")
#'
#' inducer(data = NULL, precipitant = "test")
inducer <- function(data, precipitant = "") {
  if (is.null(data))
    data <- data.frame(
      object = character(0),
      emax = numeric(0),
      ec50 = numeric(0),
      max_c = numeric(0),
      source = character(0)
    )
  new("inducer", data = data, precipitant = precipitant)
}



setValidity("inducer", function(object) {
  validate_argument(object@precipitant, "character", allow_empty = TRUE)

  # validate data
  expected_fields <- c('object', "emax", "ec50", "max_c", 'source')
  missing_fields <- setdiff(expected_fields, names(object@data))
  if (length(missing_fields) > 0)
    return(paste0("Missing fields: ", nice_enumeration(missing_fields)))

  known_objects <- c(
    "CYP1A1", "CYP1A2", "CYP2A6", "CYP2B6", "CYP2C8", "CYP2C9", "CYP2C18",
    "CYP2C19", "CYP2D6", "CYP2E1", "CYP2J2", "CYP3A4", "CYP3A5", "CYP3A7")

  unknown_objects <- setdiff(
    toupper(unique(object@data$object)), toupper(known_objects))

  if (length(unknown_objects) > 0)
    warning(paste0("Unknown objects: ", nice_enumeration(unknown_objects)))

  TRUE
})


#' Show DDI inducer object
#'
#' @param object An `inducer` object.
#'
#' @returns Nothing.
#' @export
setMethod(
  "show", "inducer",
  function(object) {
    # line <- paste0(rep("-", 5), collapse="")
    # cat(paste0(line, " CYP induction data ", line, "\n"))
    cat(paste0(hline(), " CYP induction data ", hline(), "\n"))
    cat(paste(
      "In vitro induction data for precipitant",
      object@precipitant, "\n\n"))

    out <- object@data |>
      mutate(source = case_when(
        # !is.na(source) & source != "" ~ paste0("(", source, ")"),
        !is.na(source) & source != "" ~ source,
        .default = ""
      ))

    cat(df_to_string(out, colnames = TRUE))
  }
)


#' Print DDI inducer object
#'
#' @param x An `inducer` object.
#'
#' @returns Markdown-formatted text
#' @export
setMethod(
  "print", "inducer",
  function(x) {
    caption <- ifelse(
      x@precipitant != "",
      paste0("In vitro induction data for ", x@precipitant),
      "")

    out <- x@data |>
      mutate(source = case_when(
        # !is.na(source) & source != "" ~ paste0("(", source, ")"),
        !is.na(source) & source != "" ~ source,
        .default = ""
      )) |>
      knitr::kable(caption = caption, col.names = make_labels(x@data))

    out
  }
)
