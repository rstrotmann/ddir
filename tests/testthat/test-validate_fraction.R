# validate_fraction <- ddir:::validate_fraction


test_that("validate_fraction accepts valid scalar boundary values", {
  expect_no_error(validate_fraction(0))
  expect_no_error(validate_fraction(1))
  expect_no_error(validate_fraction(0.5))
})


test_that("validate_fraction rejects non-numeric inputs", {
  x <- "0.5"
  y <- TRUE

  expect_error(validate_fraction(x), "x must be numeric!")
  expect_error(validate_fraction(y), "y must be numeric!")
})


test_that("validate_fraction rejects negative values", {
  x <- -0.01
  y <- c(0.2, -0.1)

  expect_error(validate_fraction(x), "x must be positive")
  expect_error(
    validate_fraction(y, allow_multiple = TRUE),
    "y must be positive"
  )
})


test_that("validate_fraction rejects values above one", {
  x <- 1.01
  y <- c(0.5, 1.2)

  expect_error(
    validate_fraction(x),
    "x must be between 0 and 1!"
  )
  expect_error(
    validate_fraction(y, allow_multiple = TRUE),
    "y must be between 0 and 1!"
  )
})


test_that("validate_fraction enforces single values by default", {
  x <- c(0.2, 0.8)

  expect_error(
    validate_fraction(x),
    "x must be a single value"
  )
  expect_no_error(
    validate_fraction(x, allow_multiple = TRUE)
  )
})


test_that("validate_fraction accepts valid vectors when allow_multiple is TRUE", {
  x <- c(0, 0.1, 0.5, 1)

  expect_no_error(
    validate_fraction(x, allow_multiple = TRUE)
  )
})


test_that("validate_fraction propagates NA behavior for numeric comparisons", {
  x <- NA_real_

  expect_error(
    validate_fraction(x),
    "missing value where TRUE/FALSE needed"
  )
})

