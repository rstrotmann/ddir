with_knitr <- function(expr) {
  old <- options(knitr.in.progress = TRUE)
  on.exit(options(old), add = TRUE)
  force(expr)
}


make_cyp_inh <- function(data = NULL, precipitant = "testdrug") {
  if (is.null(data)) {
    data <- tibble::tribble(
      ~object ,  ~ki, ~source,
      "CYP3A4",   10, "study"
    )
  }
  inhibition_data(data, precipitant = precipitant)
}


test_that("inhibition_data() with NULL is an empty valid object", {
  x <- inhibition_data(NULL)

  expect_s3_class(x, "inhibition_data")
  expect_s3_class(x, "data.frame")
  expect_true(is.data.frame(x))
  expect_equal(nrow(x), 0)
  expect_named(x, c("object", "ki", "source"))
  expect_identical(attr(x, "precipitant"), "")
})


test_that("inhibition_data constructor stores CYP ki data and precipitant", {
  data <- tibble::tribble(
    ~object ,  ~ki,    ~source,
    "CYP2C9",  0.6, "study 001",
    "CYP3A4", 12.5, "study 001"
  )
  x <- inhibition_data(data, precipitant = "examplinib")

  expect_s3_class(x, "inhibition_data")
  expect_equal(x$object, c("CYP2C9", "CYP3A4"))
  expect_equal(x$ki, c(0.6, 12.5))
  expect_equal(x$source, c("study 001", "study 001"))
  expect_identical(attr(x, "precipitant"), "examplinib")
})


test_that("inhibition_data constructor accepts transporter ic50 data without deriving ki", {
  data <- tibble::tribble(
    ~object  , ~ic50, ~source,
    "OATP1B1",   0.5, "study",
    "Pgp"    ,    10, "study"
  )

  expect_no_warning(x <- inhibition_data(data, precipitant = "testdrug"))
  expect_s3_class(x, "inhibition_data")
  expect_equal(x$object, c("OATP1B1", "Pgp"))
  expect_equal(x$ic50, c(0.5, 10))
  expect_false("ki" %in% names(x))
})


test_that("inhibition_data constructor keeps extra TDI columns", {
  data <- tibble::tribble(
    ~object ,  ~ki, ~kinact, ~source,
    "CYP3A4", 0.17,    0.04, "study"
  )
  x <- inhibition_data(data, precipitant = "testdrug")

  expect_equal(x$kinact, 0.04)
  expect_named(x, c("object", "ki", "kinact", "source"))
})


test_that("inhibition_data constructor accepts both ki and ic50 columns", {
  data <- tibble::tribble(
    ~object , ~ki, ~ic50, ~source,
    "CYP3A4",  10,    20, "study"
  )

  expect_no_error(x <- inhibition_data(data, precipitant = "testdrug"))
  expect_equal(x$ki, 10)
  expect_equal(x$ic50, 20)
})


test_that("inhibition_data constructor allows NA ki and empty source", {
  data <- tibble::tribble(
    ~object , ~ki, ~source,
    "CYP1A2",  NA,      ""
  )
  x <- inhibition_data(data, precipitant = "testdrug")

  expect_true(is.na(x$ki))
  expect_identical(x$source, "")
})


test_that("inhibition_data constructor rejects missing required columns", {
  no_object <- tibble::tribble(
    ~ki, ~source,
     10, "study"
  )
  no_source <- tibble::tribble(
    ~object , ~ki,
    "CYP3A4",  10
  )
  no_potency <- tibble::tribble(
    ~object , ~source,
    "CYP3A4", "study"
  )

  expect_error(inhibition_data(no_object), "Missing columns")
  expect_error(inhibition_data(no_source), "Missing columns")
  expect_error(inhibition_data(no_potency), "Either ki or ic50")
})


test_that("inhibition_data constructor rejects input without inhibition columns", {
  expect_error(inhibition_data(1), "Missing columns")
  expect_error(inhibition_data(data.frame(a = 1)), "Missing columns")
})


test_that("ensure_ki derives ki from ic50/2 and warns", {
  data <- tibble::tribble(
    ~object , ~ic50, ~source,
    "CYP3A4",    20, "study"
  )

  expect_message(
    out <- ddir:::ensure_ki(data),
    "ki derived from ic50/2 assuming substrate concentration is close to KM"
  )
  expect_equal(out$ki, 10)
  expect_equal(out$ic50, 20)
})


test_that("ensure_ki leaves an existing ki column unchanged", {
  data <- tibble::tribble(
    ~object , ~ki, ~ic50, ~source,
    "CYP3A4",  10,    20, "study"
  )

  expect_no_warning(out <- ddir:::ensure_ki(data))
  expect_equal(out$ki, 10)
  expect_equal(out$ic50, 20)
})


test_that("ensure_ki errors when neither ki nor ic50 is present", {
  data <- tibble::tribble(
    ~object , ~source,
    "CYP3A4", "study"
  )

  expect_error(
    ddir:::ensure_ki(data),
    "either a ki or an ic50"
  )
})


test_that("inhibition_data constructor keeps unknown objects without warning", {
  data <- tibble::tribble(
    ~object , ~ki, ~source,
    "CYP9Z9",  10, "study",
    "CYP3A4",   5, "study"
  )

  expect_no_warning(x <- inhibition_data(data, precipitant = "testdrug"))
  expect_equal(x$object, c("CYP9Z9", "CYP3A4"))
})


test_that("inhibition_data constructor treats known objects case-insensitively", {
  data <- tibble::tribble(
    ~object , ~ki, ~source,
    "cyp3a4",  10, "study",
    "ugt1a1",   2, "study"
  )

  expect_no_warning(inhibition_data(data, precipitant = "testdrug"))
})


test_that("dplyr_reconstruct.inhibition_data restores class and precipitant", {
  x <- make_cyp_inh()
  out <- dplyr:::dplyr_reconstruct(as.data.frame(x)[, "object", drop = FALSE], x)

  expect_s3_class(out, "inhibition_data")
  expect_named(out, "object")
  expect_identical(attr(out, "precipitant"), "testdrug")
})


test_that("dplyr row verbs reconstruct inhibition_data and keep precipitant", {
  x <- make_cyp_inh()
  filtered <- dplyr::filter(x, object == "CYP3A4")

  expect_s3_class(filtered, "inhibition_data")
  expect_identical(attr(filtered, "precipitant"), "testdrug")
  expect_equal(nrow(filtered), 1)
})


test_that("dplyr column verbs keep the inhibition_data class without re-validating", {
  x <- make_cyp_inh()

  expect_no_error(mutated <- dplyr::mutate(x, kiu = ki * 2))
  expect_no_error(slim <- dplyr::distinct(x, object))

  expect_s3_class(mutated, "inhibition_data")
  expect_s3_class(slim, "inhibition_data")
  expect_equal(mutated$kiu, 20)
  expect_named(slim, "object")
})


test_that("print.inhibition_data renders a console table", {
  data <- tibble::tribble(
    ~object , ~ki, ~source,
    "CYP3A4",  10, "study",
    "CYP2C9",  NA,      ""
  )
  x <- inhibition_data(data, precipitant = "testdrug")
  text <- paste(capture.output(print.inhibition_data(x)), collapse = "\n")

  expect_match(text, "DDI inhibition data")
  expect_match(text, "In vitro inhibition data for precipitant testdrug")
  expect_match(text, "CYP3A4")
  expect_match(text, "CYP2C9")
  expect_match(text, "study")
})


test_that("print.inhibition_data in knitr returns a kable with caption", {
  x <- make_cyp_inh()
  out <- with_knitr(print.inhibition_data(x))
  text <- paste(capture.output(print(out)), collapse = "\n")

  expect_s3_class(out, "knitr_kable")
  expect_match(text, "In vitro inhibition data for precipitant testdrug")
  expect_match(text, "CYP3A4")
})


test_that("package examplinib_cyp_inhibition is a valid inhibition_data object", {
  expect_s3_class(examplinib_cyp_inhibition, "inhibition_data")
  expect_identical(attr(examplinib_cyp_inhibition, "precipitant"), "examplinib")
  expect_true("CYP3A4" %in% examplinib_cyp_inhibition$object)
  expect_true("ki" %in% names(examplinib_cyp_inhibition))
  expect_no_error(invisible(capture.output(print.inhibition_data(examplinib_cyp_inhibition))))
})
