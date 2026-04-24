# cyp_inhibitor_data <- function(
#     target, ki, source = character(0)
# ){
#   # input validation
#   validate_argument(target, "character")
#   validate_named_vector(ki, allowed_names = c(
#     "CYP1A2", "CYP2B6", "CYP2C8", "CYP2C9", "CYP2C19", "CYP2D6", "CYP3A4"))
#   validate_argument(source, "character", allow_multiple = TRUE)
#
#   # source must either be a named vector or a singleton
#   if (!has_names(source)) {
#     if (length(source) > 1)
#       stop("Source must be a named vector or an unnamed scalar")
#     if (length(source) == 0)
#       source <- ""
#     source_table <- data.frame(
#       cyp = names(ki),
#       source = source
#     )
#   } else {
#     unknown <- setdiff(names(source), names(ki))
#     if (length(unknown) > 0) {
#       stop(paste0(
#         "source contains unknown parameter name(s): ",
#         nice_enumeration(unknown)
#       ))
#     }
#     source_table <- enframe(source, name = "cyp", value = "source") |>
#       mutate(cyp = as.character(cyp))
#   }
#
#   out <- data.frame(
#     name = target,
#     cyp = names(ki),
#     ki = ki
#   )
#
#   out <- out |>
#     left_join(source_table, by = "cyp") |>
#     mutate(source = case_when(is.na(.data$source) ~ "", .default = .data$source))
#
#   out
# }


# direct inhibition

# di <- tibble::tribble(
#        ~object,   ~target,  ~ki,     ~source,
#   "examplinib",  "CYP3A4",    1,  "source 1",
#   "examplinib",  "CYP1A2",    2,  "source 2",
#   "examplinib",  "CYP2D6",    3,          "",
#   "examplinib",  "UGT1A1",   15, "study 009",
#   "examplinib",  "UGT1A3",   15, "study 009",
#   "examplinib",  "UGT1A4",   15, "study 009",
#   "examplinib",  "UGT1A6",   15, "study 009",
#   "examplinib",  "UGT1A9",  3.8, "study 009",
#   "examplinib",  "UGT2B7",   15, "study 009",
#   "examplinib", "UGT2B15",   15, "study 009",
#   "examplinib", "UGT2B17",  6.1, "study 009",
#   "examplinib",     "Pgp", 0.41, "study 005",
#   "examplinib",    "BCRP",  1.9, "study 005",
#   "examplinib",    "OCT1",  2.3, "study 006",
#   "examplinib", "OATP1B1",  177, "study 006",
#   "examplinib", "OATP1B3",   35, "study 006",
#   "examplinib",    "OAT1",  271,          "",
#   "examplinib",    "OAT3",  300,          "",
#   "examplinib",    "BSEP", 12.8,          "",
#   "examplinib",    "OCT2",   67, "study 006",
#   "examplinib",   "MATE1",  3.6, "study 006",
#   "examplinib",  "MATE2k",  1.1, "study 006"
# )
#
#
# di <- tibble::tribble(
#     ~target,  ~ki,     ~source,
#    "CYP3A4",    1,  "source 1",
#    "CYP1A2",    2,  "source 2",
#    "CYP2D6",    3,          "",
#    "UGT1A1",   15, "study 009",
#    "UGT1A3",   15, "study 009",
#    "UGT1A4",   15, "study 009",
#    "UGT1A6",   15, "study 009",
#    "UGT1A9",  3.8, "study 009",
#    "UGT2B7",   15, "study 009",
#   "UGT2B15",   15, "study 009",
#   "UGT2B17",  6.1, "study 009",
#       "Pgp", 0.41, "study 005",
#      "bcrp",  1.9, "study 005",
#      "OCT1",  2.3, "study 006",
#   "OATP1B1",  177, "study 006",
#   "OATP1B3",   35, "study 006",
#      "OAT1",  271,          "",
#      "OAT3",  300,          "",
#      "BSEP", 12.8,          "",
#      "OCT2",   67, "study 006",
#     "MATE1",  3.6, "study 006",
#    "MATE2k",  1.1, "study 006"
# )


#' Direct inhibitor class definition
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


#' Title
#'
#' @param data
#'
#' @returns
#' @export
#'
#' @examples
inhibitor <- function(data, object = "") {
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


#' Title
#'
#' @param inhibitor
#'
#' @returns
#' @export
#'
#' @examples
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


#' Title
#'
#' @param inhibitor
#'
#' @returns
#' @export
#'
#' @examples
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
