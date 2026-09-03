
#' Validate precipitant object
#'
#' @param obj A precipitant object
#'
#' @returns Nothing.
#' @noRd
validate_precipitant <- function(object) {
  obj_name <- deparse(substitute(object))

  if (!inherits(object, "precipitant"))
    stop(paste(obj_name, "must be a precipitant object"))

  validate_argument(object$name, "character", allow_empty = TRUE)
  validate_argument(object$oral, "logical")
  validate_argument(object$mw, "numeric", expect_positive = TRUE)
  validate_argument(object$dose, "numeric", expect_positive = TRUE, allow_na = TRUE)
  validate_argument(object$imaxss, "numeric", expect_positive = TRUE)
  validate_argument(object$fu, "fraction")
  validate_argument(object$fumic, "fraction")
  validate_argument(object$rb, "numeric", expect_positive = TRUE)
  validate_argument(object$fa, "fraction", allow_na = TRUE)
  validate_argument(object$fg, "fraction", allow_na = TRUE)
  validate_argument(object$ka, "numeric", expect_positive = TRUE, allow_na = TRUE)
  validate_argument(object$solubility, "numeric", expect_positive = TRUE)
  validate_argument(object$source, "character", allow_multiple = TRUE, allow_empty = TRUE)

  validate_named_vector(object$source, allowed_names = names(object))

  invisible(NULL)
}


#' Validate inhibition_data object
#'
#' @param allow_null Logical.
#' @param obj Inhibitor object
#' @param allowed_objects Objects that will be retained in the output, defaults
#' to all, if NULL.
#'
#' @returns The inhibition_data object, filtered to the allowed objects.
#' @noRd
validate_inhibition_data <- function(
    obj,
    allowed_objects = NULL,
    allow_null = FALSE
  ){
  # early fail on NULL, if applicable
  if (isTRUE(allow_null)) {
    if (is.null(obj))
      return(obj)
  }

  obj_name <- deparse(substitute(obj))

  if (!inherits(obj, "inhibition_data"))
    stop(paste(obj_name, "must be an inhibition_data object"))

  # check required fields
  expected_fields <- c("object", "source")
  missing_fields <- setdiff(expected_fields, names(obj))
  if (length(missing_fields) > 0) {
    stop(paste0(
      "Missing columns in ", obj_name, ": ",
      nice_enumeration(missing_fields)
    ))
  }
  if (!any(c("ki", "ic50") %in% names(obj))) {
    stop(paste("Either ki or ic50 must be in", obj_name))
  }

  # filter for allowed objects
  if (!is.null(allowed_objects)) {
    unexpected_objects <- setdiff(unique(obj$object), allowed_objects)
    if (length(unexpected_objects) > 0) {
      warning(paste0(
        "Unexpected objects in ", obj_name, " removed: ",
        nice_enumeration(unexpected_objects)
      ))
    }

    out <- filter(obj, .data$object %in% allowed_objects)
  } else {
    out <- obj
  }

  out
}


#' Validate induction_data object
#'
#' @param obj Induction_data object.
#' @param allowed_objects Objects that will be retained in the output, defaults
#' to all, if NULL.
#' @param allow_null Logical.
#'
#' @returns Nothing.
#' @noRd
validate_induction_data <- function(
    obj,
    allowed_objects = NULL,
    allow_null = FALSE
) {
  if (isTRUE(allow_null)) {
    if (is.null(obj))
      # return(invisible(NULL))
      return(obj)
  }

  obj_name <- deparse(substitute(obj))

  if (!inherits(obj, "induction_data"))
    stop(paste(obj_name, "must be an induction_data object"))

  expected_fields <- c("object", "emax", "ec50", "max_c", "source")
  missing_fields <- setdiff(expected_fields, names(obj))
  if (length(missing_fields) > 0) {
    stop(paste0(
      "Missing columns in ", obj_name, ": ",
      nice_enumeration(missing_fields)
    ))
  }

  if (!is.null(allowed_objects)) {
    unexpected_objects <- setdiff(unique(obj$object), allowed_objects)
    if (length(unexpected_objects) > 0) {
      warning(paste0(
        "Unexpected objects in ", obj_name, " removed: ",
        nice_enumeration(unexpected_objects)
      ))
    }
    out <- filter(obj, .data$object %in% allowed_objects)
  } else {
    out <- obj
  }

  out
}


#' validate function argument
#'
#' @param param The argument.
#' @param type The expected parameter type (one of 'character', 'logical',
#' 'numeric', 'fraction' or 'date').
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
    type = c("character", "logical", "numeric", "fraction", "date"),
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
      (type == "fraction" && !is.numeric(param)) ||
      (type == "date" && !inherits(param, "Date"))) {
    stop(paste0(param_name, " must be a ", type, " value"))
  }

  # fraction
  if (type == "fraction") {
    if (any(param < 0 | param > 1, na.rm = TRUE))
      stop(paste0(param_name, " must be between 0 and 1!"))
  }

  # positive
  if (type == "numeric" && expect_positive) {
    if (any(param < 0, na.rm = TRUE))
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
    any(nchar(param) == 0, na.rm = TRUE)
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


#' Validate function argument as data frame
#'
#' @param param The argument.
#' @param expected_fields Expected columns in the data frame.
#' @param allow_null Allow NULL.
#'
#' @returns Nothing or stop.
#' @noRd
validate_df_argument <- function(
    param,
    expected_fields = NULL,
    allow_null = FALSE
) {
  param_name <- deparse(substitute(param))

  # Check for NULL first
  if (is.null(param)) {
    if (allow_null) {
      return(invisible(NULL))
    } else {
      stop(paste0(param_name, " must not be NULL"))
    }
  }

  if (!inherits(param, "data.frame"))
    stop(paste0(param_name, " must be a data.frame"))

  missing_fields <- setdiff(expected_fields, names(param))
  if (length(missing_fields) > 0) {
    stop(paste0(
      "Missing columns in ", param_name, ": ",
      nice_enumeration(missing_fields)
    ))
  }
  invisible(NULL)
}


#' Test whether object is named
#'
#' @param obj Scalar or vector.
#'
#' @returns Logical.
#' @noRd
has_names <- function(obj) {
  obj_names <- names(obj)
  !is.null(obj_names)
}


#' Validate named vector
#'
#' @param param Named vector
#' @param allowed_names Allowed names as character.
#'
#' @returns Nothing.
#' @noRd
validate_named_vector <- function(
    param, allowed_names = NULL
) {
  # input validation
  validate_argument(allowed_names, "character", allow_null = TRUE,
                    allow_multiple = TRUE)

  param_name <- deparse(substitute(param))

  if (length(param) > 0) {
    if (is.null(names(param)) || any(names(param) == "")) {
      stop(paste0(param_name, " must be a named vector"))
    }

    if (!is.null(allowed_names)) {
      unknown <- setdiff(names(param), allowed_names)
      if (length(unknown) > 0) {
        stop(paste0(
          "allowed_names contains unknown parameter name(s): ",
          nice_enumeration(unknown)
        ))
      }
    }
  }
}

