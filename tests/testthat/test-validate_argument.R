test_that("validate_argument accepts valid scalar inputs by type", {
  expect_no_error(validate_argument("abc", type = "character"))
  expect_no_error(validate_argument(TRUE, type = "logical"))
  expect_no_error(validate_argument(FALSE, type = "logical"))
  expect_no_error(validate_argument(1.5, type = "numeric"))
  expect_no_error(validate_argument(1L, type = "numeric"))
  expect_no_error(validate_argument(0, type = "fraction"))
  expect_no_error(validate_argument(1, type = "fraction"))
  expect_no_error(validate_argument(0.5, type = "fraction"))
  expect_no_error(validate_argument(as.Date("2025-01-01"), type = "date"))
})


test_that("validate_argument defaults type to character via match.arg", {
  expect_no_error(validate_argument("abc"))
  x <- 1
  expect_error(validate_argument(x), "x must be a character value")
})


test_that("validate_argument rejects unknown types", {
  expect_error(validate_argument("abc", type = "foo"), "arg")
})


test_that("validate_argument returns invisible NULL", {
  expect_null(validate_argument("abc", type = "character"))
  expect_invisible(validate_argument("abc", type = "character"))
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
  expect_no_error(
    validate_argument(x, type = "numeric", allow_null = TRUE)
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

  y <- NA_character_
  expect_error(
    validate_argument(y, type = "character"),
    "y must not contain NA"
  )
  expect_no_error(
    validate_argument(y, type = "character", allow_na = TRUE)
  )

  z <- NA
  expect_error(
    validate_argument(z, type = "logical"),
    "z must not contain NA"
  )
  expect_no_error(
    validate_argument(z, type = "logical", allow_na = TRUE)
  )

  d <- as.Date(NA)
  expect_error(
    validate_argument(d, type = "date"),
    "d must not contain NA"
  )
  expect_no_error(
    validate_argument(d, type = "date", allow_na = TRUE)
  )
})


test_that("validate_argument treats NaN as NA", {
  x <- NaN
  expect_error(
    validate_argument(x, type = "numeric"),
    "x must not contain NA"
  )
  expect_no_error(
    validate_argument(x, type = "numeric", allow_na = TRUE)
  )
})


test_that("validate_argument reports NA before a type mismatch", {
  x <- NA_character_
  expect_error(
    validate_argument(x, type = "numeric"),
    "x must not contain NA"
  )
})


test_that("validate_argument rejects wrong type values", {
  x <- 1
  y <- "TRUE"
  z <- TRUE
  d <- "2025-01-01"
  f <- factor("a")
  p <- as.POSIXct("2025-01-01")
  lst <- list("a")

  expect_error(validate_argument(x, type = "character"), "x must be a character value")
  expect_error(validate_argument(y, type = "logical"), "y must be a logical value")
  expect_error(validate_argument(z, type = "numeric"), "z must be a numeric value")
  expect_error(validate_argument(d, type = "date"), "d must be a date value")
  expect_error(validate_argument(f, type = "character"), "f must be a character value")
  expect_error(validate_argument(p, type = "date"), "p must be a date value")
  expect_error(validate_argument(lst, type = "character"), "lst must be a character value")
  expect_error(validate_argument(y, type = "numeric"), "y must be a numeric value")
  expect_error(validate_argument(z, type = "fraction"), "z must be a fraction value")
})


test_that("validate_argument enforces fraction range", {
  x <- -0.01
  y <- 1.01
  z <- c(0.5, 1.2)
  inf <- Inf

  expect_error(
    validate_argument(x, type = "fraction"),
    "x must be between 0 and 1!"
  )
  expect_error(
    validate_argument(y, type = "fraction"),
    "y must be between 0 and 1!"
  )
  expect_error(
    validate_argument(z, type = "fraction", allow_multiple = TRUE),
    "z must be between 0 and 1!"
  )
  expect_error(
    validate_argument(inf, type = "fraction"),
    "inf must be between 0 and 1!"
  )
  expect_no_error(
    validate_argument(c(0, 0.5, 1), type = "fraction", allow_multiple = TRUE)
  )
})


test_that("validate_argument allows NA in fractions when allow_na is TRUE", {
  x <- NA_real_
  y <- c(0.2, NA_real_)

  expect_no_error(
    validate_argument(x, type = "fraction", allow_na = TRUE)
  )
  expect_no_error(
    validate_argument(y, type = "fraction", allow_multiple = TRUE, allow_na = TRUE)
  )
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


test_that("validate_argument treats Inf as positive and -Inf as negative", {
  expect_no_error(
    validate_argument(Inf, type = "numeric", expect_positive = TRUE)
  )

  x <- -Inf
  expect_error(
    validate_argument(x, type = "numeric", expect_positive = TRUE),
    "x must be positive"
  )
})


test_that("validate_argument allows NA with expect_positive when allow_na is TRUE", {
  x <- NA_real_
  y <- c(1, NA_real_)

  expect_no_error(
    validate_argument(x, type = "numeric", expect_positive = TRUE, allow_na = TRUE)
  )
  expect_no_error(
    validate_argument(
      y,
      type = "numeric",
      expect_positive = TRUE,
      allow_multiple = TRUE,
      allow_na = TRUE
    )
  )
})


test_that("validate_argument ignores expect_positive unless type is numeric", {
  x <- -0.1
  expect_error(
    validate_argument(x, type = "fraction", expect_positive = TRUE),
    "x must be between 0 and 1!"
  )
  expect_no_error(
    validate_argument(-1, type = "numeric", expect_positive = FALSE)
  )
})


test_that("validate_argument enforces single values by default", {
  x <- c("a", "b")
  y <- c(TRUE, FALSE)
  z <- c(1, 2)
  d <- as.Date(c("2025-01-01", "2025-01-02"))

  expect_error(
    validate_argument(x, type = "character"),
    "x must be a single value"
  )
  expect_error(
    validate_argument(y, type = "logical"),
    "y must be a single value"
  )
  expect_error(
    validate_argument(z, type = "numeric"),
    "z must be a single value"
  )
  expect_error(
    validate_argument(d, type = "date"),
    "d must be a single value"
  )
  expect_no_error(
    validate_argument(x, type = "character", allow_multiple = TRUE)
  )
  expect_no_error(
    validate_argument(y, type = "logical", allow_multiple = TRUE)
  )
  expect_no_error(
    validate_argument(z, type = "numeric", allow_multiple = TRUE)
  )
  expect_no_error(
    validate_argument(d, type = "date", allow_multiple = TRUE)
  )
})


test_that("validate_argument rejects empty vectors unless allow_multiple is TRUE", {
  x <- character(0)
  y <- numeric(0)

  expect_error(
    validate_argument(x, type = "character"),
    "x must be a single value"
  )
  expect_no_error(
    validate_argument(x, type = "character", allow_multiple = TRUE)
  )
  expect_no_error(
    validate_argument(y, type = "numeric", allow_multiple = TRUE)
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
    validate_argument(x, type = "character", allow_empty = TRUE)
  )
  expect_no_error(
    validate_argument(y, type = "character", allow_multiple = TRUE, allow_empty = TRUE)
  )
})


test_that("validate_argument treats whitespace as a non-empty string", {
  expect_no_error(validate_argument(" ", type = "character"))
})


test_that("validate_argument allows NA in character vectors when allow_na is TRUE", {
  x <- c("ok", NA_character_)

  expect_no_error(
    validate_argument(x, type = "character", allow_multiple = TRUE, allow_na = TRUE)
  )
  z <- c("ok", NA_character_, "")
  expect_error(
    validate_argument(z, type = "character", allow_multiple = TRUE, allow_na = TRUE),
    "z must be a non-empty string"
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
  expect_no_error(
    validate_argument(x, type = "character", values = c("a", "b", "c"))
  )
  z <- "d"
  expect_error(
    validate_argument(z, type = "character", values = c("a", "b", "c")),
    "z must be a, b or c!"
  )
})


test_that("validate_argument enforces allowed numeric values", {
  x <- 1
  y <- 3

  expect_no_error(
    validate_argument(x, type = "numeric", values = c(1, 2))
  )
  expect_error(
    validate_argument(y, type = "numeric", values = c(1, 2)),
    "y must be 1 or 2!"
  )
})


test_that("validate_argument checks values after allowing NA", {
  x <- NA_character_
  expect_error(
    validate_argument(
      x,
      type = "character",
      allow_na = TRUE,
      values = c("a", "b")
    ),
    "x must be a or b!"
  )
  expect_no_error(
    validate_argument(
      x,
      type = "character",
      allow_na = TRUE,
      values = c("a", NA_character_)
    )
  )
})


test_that("validate_argument uses the argument expression in error messages", {
  expect_error(
    validate_argument(c("a", "b"), type = "character"),
    'c\\("a", "b"\\) must be a single value'
  )
  expect_error(
    validate_argument(123, type = "character"),
    "123 must be a character value"
  )
})
