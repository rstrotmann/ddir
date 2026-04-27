#' DDI CYP inducer data object
#'
#' @slot object The name of the DDI perpetrator as character.
#' @slot data data.frame.
#'
#' @returns Inducer object
#' @export
setClass(
  Class = "inducer",
  representation(
    object = "character",
    data = "data.frame"
  ),
  prototype(
    object = "",
    data = data.frame(
      target = character(0),
      emax = numeric(0),
      ec50 = numeric(0),
      max_c = numeric(0),
      source = character(0)
    )
  )
)


#' DDI CYP inducer object constructor
#'
#' This class constructor expects as the as the data argument a data frome with
#' the following columns:
#'
#' * target: The CYP induction target. Must be one of "CYP1A1", "CYP1A2",
#' "CYP2A6", "CYP2B6", "CYP2C8", "CYP2C9", "CYP2C18", "CYP2C19", "CYP2D6",
#' "CYP2E1", "CYP2J2", "CYP3A4", "CYP3A5", or "CYP3A7".
#' * emax: The maximal induction effect.
#' * ec50: The EC50, i.e., the concentration causing the half-maximal induction
#' effect.
#' * max_c: The maximal concentration tested in the assay in uM.
#' * source: Source information for the parameter.
#'
#' @param data Data frame with the columns 'target', 'emax', 'ec50', 'max_c' and
#'   'source'.
#' @param object The name of the DDI perpetrator as character.
#'
#' @returns Inducer object
#' @export
inducer <- function(data, object = "") {
  if (is.null(data))
    data <- data.frame(
      target = character(0),
      emax = numeric(0),
      ec50 = numeric(0),
      max_c = numeric(0),
      source = character(0)
    )
  new("inducer", data = data, object = object)
}



setValidity(
  "inducer",
  function(object) {
    validate_argument(object@object, "character", allow_empty = TRUE)

    # validate data
    expected_fields <- c('target', "emax", "ec50", "max_c", 'source')
    missing_fields <- setdiff(expected_fields, names(object@data))
    if (length(missing_fields) > 0)
      return(paste0("Missing fields: ", nice_enumeration(missing_fields)))

    known_targets <- c(
      "CYP1A1", "CYP1A2", "CYP2A6", "CYP2B6", "CYP2C8", "CYP2C9", "CYP2C18",
      "CYP2C19", "CYP2D6", "CYP2E1", "CYP2J2", "CYP3A4", "CYP3A5", "CYP3A7")

    unknown_targets <- setdiff(
      toupper(unique(object@data$target)), toupper(known_targets))
    if (length(unknown_targets) > 0)
      warning(paste0("Unknown targets: ", nice_enumeration(unknown_targets)))

    TRUE
})


#' Show DDI inducer object
#'
#' @param inducer DDI inducer object
#'
#' @returns Nothing.
#' @export
setMethod(
  "show", "inducer",
  function(object) {
    line <- paste0(rep("-", 5), collapse="")
    cat(paste0(line, " CYP induction data ", line, "\n"))
    cat(paste("In vitro induction data for object", object@object, "\n\n"))
    out <- object@data |>
      mutate(source = case_when(
        !is.na(source) & source != "" ~ paste0("(", source, ")"),
        .default = ""
      ))

    cat(df_to_string(out, colnames = TRUE))
  }
)


#' Print DDI inducer object
#'
#' @param inducer DDI inducer object
#'
#' @returns Markdown-formatted text
#' @export
setMethod(
  "print", "inducer",
  function(x) {
    caption <- ifelse(
      x@object != "",
      paste0("In vitro induction data for ", x@object),
      "")

    out <- x@data |>
      mutate(source = case_when(
        !is.na(source) & source != "" ~ paste0("(", source, ")"),
        .default = ""
      )) |>
      kable(caption = caption, col.names = make_labels(x@data))

    out
  }
)
