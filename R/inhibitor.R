#' Inhibitor class definition
#'
#' @slot object Name of perpetrator.
#' @slot data Data frame with the fields 'object', 'ki' and 'source'.
#'
#' @returns An inhibitor class object.
#' @export
setClass(
  Class = "inhibitor",
  representation(
    object = "character",
    data = "data.frame"
  ),
  prototype(
    object = "",
    data = data.frame(
      target = character(0),
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
#' * target: The CYP or UGT enzyme target as character.
#' * ki: The ki in uM.
#' * source: Source information (e.g., study name) as character.
#'
#' ### Time-dependent CYP inhibition
#'
#' * target: The CYP or UGT enzyme target as character.
#' * ki: The ki in uM.
#' * kinact: The kinact in 1/h.
#' * source: Source information (e.g., study name) as character.
#'
#' ### Transporter inhibition
#'
#' * target: The transporter target as character.
#' * ic50: The IC50 in uM.
#' * source: Source information (e.g., study name) as character.
#'
#' @param object Character
#'
#' @returns Inhibitor object
#' @export
inhibitor <- function(data, object = "") {
  if (is.null(data))
    data <- data.frame(
      target = character(0),
      ki = numeric(0),
      kinact = numeric(0),
      source = character(0)
    )
  new("inhibitor", data = data, object = object)
}



#' Validity check for inhibitor class
#'
#' @param object Inhibitor object.
#'
#' @returns Logical.
#' @noRd
setValidity("inhibitor", function(object) {
  validate_argument(object@object, "character", allow_empty = TRUE)

  # validate data
  expected_fields <- c('target', 'source')
  missing_fields <- setdiff(expected_fields, names(object@data))
  if (length(missing_fields) > 0)
    return(paste0("Missing fields: ", nice_enumeration(missing_fields)))
  expected_ki <- c("ic50", "ki") %in% names(object@data)
  # print(expected_ki)
  if (length(which(expected_ki == TRUE)) != 1)
    return(paste0("Exactly one of ki or ic50 fields expected"))

  known_targets <- c(
    c("CYP1A1", "CYP1A2", "CYP2A6", "CYP2B6", "CYP2C8", "CYP2C9", "CYP2C18",
      "CYP2C19", "CYP2D6", "CYP2E1", "CYP2J2", "CYP3A4", "CYP3A5", "CYP3A7"),
    c("UGT1A1", "UGT1A3", "UGT1A4", "UGT1A6", "UGT1A9", "UGT2B7", "UGT2B15",
      "UGT2B17"),
    c("Pgp", "BCRP", "OATP1B1", "OATP1B3", "OAT1", "OAT3", "BSEP", "OCT1",
      "OCT2", "MATE1", "MATE2k")
  )
  unknown_targets <- setdiff(
    toupper(unique(object@data$target)), toupper(known_targets))
  if (length(unknown_targets) > 0)
    warning(paste0("Unknown targets: ", nice_enumeration(unknown_targets)))
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
    line <- paste0(rep("-", 5), collapse="")
    cat(paste0(line, " DDI inhibition set ", line, "\n"))
    cat(paste("In vitro inhibition data for object", object@object, "\n\n"))
    out <- object@data |>
      mutate(source = case_when(
        !is.na(source) ~ paste0("(", source, ")"),
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
      x@object != "",
      paste0("In vitro inhibition data for ", x@object),
      "")

    col_names <- make_labels(x@data)

    out <- x@data |>
      mutate(source = case_when(
        !is.na(source) ~ paste0("(", source, ")"),
        .default = ""
      )) |>
      # kable(caption = caption, col.names = c("Target", "$K_{i}$", "Source"))
      kable(caption = caption, col.names = col_names)

    out
  }
)
