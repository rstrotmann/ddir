#' validate function argument
#'
#' @param param The argument.
#' @param type The expected parameter type (one of 'character', 'logical',
#' 'numeric' or 'date').
#' @param allow_null Allow NULL values.
#' @param allow_empty Allow empty values.
#' @param allow_multiple Allow multiple values.
#' @param allow_na Allow NA values.
#' @param values Allowed values.
#' @param expect_positive Numerical values must be positive.
#'
#' @returns Nothing or stop.
#' @noRd
validate_argument <- function(
    param,
    type = c("character", "logical", "numeric", "date"),
    allow_null = FALSE,
    allow_empty = FALSE,
    allow_multiple = FALSE,
    allow_na = FALSE,
    expect_positive = FALSE,
    values = NULL
) {
  # Validate type parameter
  type <- match.arg(type)

  param_name <- deparse(substitute(param))

  # Check for NULL first
  if (is.null(param)) {
    if (allow_null) {
      return(invisible(NULL))
    } else {
      stop(paste0(param_name, " must not be NULL"))
    }
  }

  # Check for NA values
  if (!allow_na && any(is.na(param))) {
    stop(paste0(param_name, " must not contain NA"))
  }

  # Type checking
  if ((type == "character" && !is.character(param)) ||
      (type == "logical" && !is.logical(param)) ||
      (type == "numeric" && !is.numeric(param)) ||
      (type == "date" && !is.Date(param))) {
    stop(paste0(param_name, " must be a ", type, " value"))
  }

  # positive
  if (type == "numeric" && expect_positive) {
    if (any(param < 0))
      stop(paste0(param_name, " must be positive"))
  }

  # Length checking
  if (length(param) != 1 && !allow_multiple) {
    stop(paste0(param_name, " must be a single value"))
  }

  # Empty string check (only for character types)
  if (
    type == "character" &&
    !allow_empty &&
    length(param) > 0 &&
    any(nchar(param) == 0)
  ) {
    stop(paste0(param_name, " must be a non-empty string"))
  }

  if (!is.null(values)) {
    if (!all(param %in% values))
      stop(paste0(
        param_name, " must be ",
        nice_enumeration(values, conjunction = "or"), "!"
      ))
  }

  invisible(NULL)
}


#' Validate fraction argument
#'
#' @param param Numeric.
#' @param allow_multiple Logical.
#'
#' @returns Nothing
#' @noRd
validate_fraction <- function(param, allow_multiple = FALSE) {
  # Validate type parameter
  param_name <- deparse(substitute(param))
  if (!is.numeric(param))
    stop(paste0(param_name, " must be numeric!"))

  if (any(param < 0))
    stop(paste0(param_name, " must be positive"))

  if (any(param < 0 | param > 1))
    stop(paste0(param_name, " must be between 0 and 1!"))

  # Length checking
  if (length(param) != 1 && !allow_multiple) {
    stop(paste0(param_name, " must be a single value"))
  }
}
