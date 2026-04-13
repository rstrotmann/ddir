# validate_argument <- ddir:::validate_argument
# environment(validate_argument) <- environment()

if (!exists("is.Date", mode = "function")) {
  is.Date <- function(x) inherits(x, "Date")
}


test_that("validate_argument accepts valid scalar inputs by type", {
  expect_no_error(validate_argument("abc", type = "character"))
  expect_no_error(validate_argument(TRUE, type = "logical"))
  expect_no_error(validate_argument(1.5, type = "numeric"))
  # expect_no_error(validate_argument(as.Date("2025-01-01"), type = "date"))
})


test_that("validate_argument handles NULL according to allow_null", {
  x <- NULL
  expect_error(
    validate_argument(x, type = "character"),
    "x must not be NULL"
  )
  expect_no_error(
    validate_argument(x, type = "character", allow_null = TRUE)
  )
})


test_that("validate_argument handles NA according to allow_na", {
  x <- c(1, NA_real_)
  expect_error(
    validate_argument(x, type = "numeric", allow_multiple = TRUE),
    "x must not contain NA"
  )
  expect_no_error(
    validate_argument(
      x,
      type = "numeric",
      allow_multiple = TRUE,
      allow_na = TRUE
    )
  )
})


test_that("validate_argument rejects wrong type values", {
  x <- 1
  y <- "TRUE"
  z <- TRUE
  d <- "2025-01-01"

  expect_error(validate_argument(x, type = "character"), "x must be a character value")
  expect_error(validate_argument(y, type = "logical"), "y must be a logical value")
  expect_error(validate_argument(z, type = "numeric"), "z must be a numeric value")
  # expect_error(validate_argument(d, type = "date"), "d must be a date value")
})


test_that("validate_argument enforces numeric positivity only for negatives", {
  x <- -0.1
  y <- 0
  z <- c(1, -2)

  expect_error(
    validate_argument(x, type = "numeric", expect_positive = TRUE),
    "x must be positive"
  )
  expect_no_error(
    validate_argument(y, type = "numeric", expect_positive = TRUE)
  )
  expect_error(
    validate_argument(
      z,
      type = "numeric",
      expect_positive = TRUE,
      allow_multiple = TRUE
    ),
    "z must be positive"
  )
})


test_that("validate_argument enforces single values by default", {
  x <- c("a", "b")
  expect_error(
    validate_argument(x, type = "character"),
    "x must be a single value"
  )
  expect_no_error(
    validate_argument(x, type = "character", allow_multiple = TRUE)
  )
})


test_that("validate_argument enforces non-empty strings by default", {
  x <- ""
  y <- c("ok", "")

  expect_error(
    validate_argument(x, type = "character"),
    "x must be a non-empty string"
  )
  expect_error(
    validate_argument(y, type = "character", allow_multiple = TRUE),
    "y must be a non-empty string"
  )
  expect_no_error(
    validate_argument(y, type = "character", allow_multiple = TRUE, allow_empty = TRUE)
  )
})


test_that("validate_argument enforces allowed values when provided", {
  x <- "a"
  y <- c("a", "b")

  expect_no_error(
    validate_argument(x, type = "character", values = c("a", "b"))
  )
  expect_no_error(
    validate_argument(
      y,
      type = "character",
      allow_multiple = TRUE,
      values = c("a", "b", "c")
    )
  )
  expect_error(
    validate_argument(x, type = "character", values = c("b", "c")),
    "x must be b or c!"
  )
  expect_error(
    validate_argument(
      y,
      type = "character",
      allow_multiple = TRUE,
      values = c("a")
    ),
    "y must be a!"
  )
})

