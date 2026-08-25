#' Inhibitor class definition
#'
#' @slot precipitant Name of DDI precipitant.
#' @slot data Data frame with the fields 'object', 'ki' and 'source'.
#'
#' @returns An inhibitor class object.
#' @export
setClass(
  Class = "inhibitor",
  representation(
    precipitant = "character",
    data = "data.frame"
  ),
  prototype(
    precipitant = "",
    data = data.frame(
      object = character(0),
      ki = numeric(0),
      source = character(0)
    )
  )
)


#' DDI inhibition data object constructor.
#'
#' @param data Data frame with the following fields:
#'
#' ### Direct CYP or UGT inhibition
#'
#' * object: The CYP or UGT enzyme target as character.
#' * ki: The ki in uM.
#' * source: Source information (e.g., study name) as character.
#'
#' ### Time-dependent CYP inhibition
#'
#' * object: The CYP or UGT enzyme target as character.
#' * ki: The ki in uM.
#' * kinact: The kinact in 1/h.
#' * source: Source information (e.g., study name) as character.
#'
#' ### Transporter inhibition
#'
#' * object: The transporter target as character.
#' * ic50: The IC50 in uM.
#' * source: Source information (e.g., study name) as character.
#' @param precipitant The name of the precipitant.
#'
#' @returns Inhibitor object
#' @export
#' @examples
#' inhibitor(
#'   data = tibble::tribble(
#'   ~object, ~ki, ~source,
#'   "CYP2C8", 11, "study 001",
#'   "CYP2C9", 0.6, "study 001",
#'   "CYP2C19", 0.25, "study 001"
#'   ),
#'   precipitant = "examplinib"
#' )
#'
inhibitor <- function(data, precipitant = "") {
  if (is.null(data))
    data <- data.frame(
      object = character(0),
      ki = numeric(0),
      kinact = numeric(0),
      source = character(0)
    )
  new("inhibitor", data = data, precipitant = precipitant)
}


#' Validity check for inhibitor class
#'
#' @param object Inhibitor object.
#'
#' @returns Logical.
#' @noRd
setValidity("inhibitor", function(object) {
  validate_argument(object@precipitant, "character", allow_empty = TRUE)

  # validate data
  expected_fields <- c('object', 'source')
  missing_fields <- setdiff(expected_fields, names(object@data))
  if (length(missing_fields) > 0)
    return(paste0("Missing fields: ", nice_enumeration(missing_fields)))

  expected_ki <- c("ic50", "ki") %in% names(object@data)
  if (length(which(expected_ki == TRUE)) != 1)
    return(paste0("Exactly one of ki or ic50 fields expected"))

  known_objects <- c(
    c("CYP1A1", "CYP1A2", "CYP2A6", "CYP2B6", "CYP2C8", "CYP2C9", "CYP2C18",
      "CYP2C19", "CYP2D6", "CYP2E1", "CYP2J2", "CYP3A4", "CYP3A5", "CYP3A7"),
    c("UGT1A1", "UGT1A3", "UGT1A4", "UGT1A6", "UGT1A9", "UGT2B7", "UGT2B15",
      "UGT2B17"),
    c("Pgp", "BCRP", "OATP1B1", "OATP1B3", "OAT1", "OAT3", "BSEP", "OCT1",
      "OCT2", "MATE1", "MATE2k")
  )
  unknown_objects <- setdiff(
    toupper(unique(object@data$object)), toupper(known_objects))
  if (length(unknown_objects) > 0)
    warning(paste0("Unknown objects: ", nice_enumeration(unknown_objects)))
  TRUE
})


#' Show DDI inhibitor object.
#'
#' @param inhibitor Inhibitor object.
#'
#' @returns Nothing.
#' @export
setMethod(
  "show", "inhibitor",
  function(object) {
    cat(paste0(hline(), " DDI inhibition data ", hline(), "\n"))

    cat(paste(
      "In vitro inhibition data for precipitant", object@precipitant, "\n\n"))
    out <- object@data |>
      mutate(source = case_when(
        # !is.na(source) & !source == "" ~ paste0("(", source, ")"),
        !is.na(source) & !source == "" ~ source,
        .default = ""
    ))

    cat(df_to_string(out, colnames = TRUE))
  }
)


#' Print DDI inhibitor object.
#'
#' @param inhibitor Inhibitor object.
#'
#' @returns Markdown-formatted text.
#' @export
setMethod(
  "print", "inhibitor",
  function(x) {
    caption <- ifelse(
      x@precipitant != "",
      paste0("In vitro inhibition data for ", x@precipitant),
      ""
    )

    col_names <- make_labels(x@data)
    out <- x@data |>
      mutate(source = case_when(
        # !is.na(source) & source != "" ~ paste0("(", source, ")"),
        !is.na(source) & source != "" ~ source,
        .default = ""
      )) |>
      knitr::kable(caption = caption, col.names = col_names)

    out
  }
)
