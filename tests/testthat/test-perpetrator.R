make_test_perpetrator <- function(...) {
  defaults <- list(
    name = "examplinib",
    oral = TRUE,
    mw = 492.6,
    dose = 450,
    imaxss = 3530,
    fu = 0.023,
    fumic = 1,
    rb = 1,
    fa = 0.81,
    fg = 1,
    ka = 0.00267,
    solubility = Inf,
    source = c(dose = "clinical dose", imaxss = "study 001", fu = "study 002")
  )
  do.call(perpetrator, utils::modifyList(defaults, list(...)))
}


test_that("perpetrator constructor populates all slots", {
  x <- make_test_perpetrator()

  expect_s4_class(x, "perpetrator")
  expect_identical(x@name, "examplinib")
  expect_identical(x@oral, TRUE)
  expect_equal(x@mw, 492.6)
  expect_equal(x@dose, 450)
  expect_equal(x@imaxss, 3530)
  expect_equal(x@fu, 0.023)
  expect_equal(x@fumic, 1)
  expect_equal(x@rb, 1)
  expect_equal(x@fa, 0.81)
  expect_equal(x@fg, 1)
  expect_equal(x@ka, 0.00267)
  expect_identical(x@solubility, Inf)
  expect_identical(
    unname(x@source[c("dose", "imaxss", "fu")]),
    c("clinical dose", "study 001", "study 002")
  )
})


test_that("perpetrator constructor applies defaults", {
  x <- perpetrator(
    name = "default-test",
    oral = TRUE,
    mw = 500,
    dose = 100,
    imaxss = 1000
  )

  expect_equal(x@fu, 1)
  expect_equal(x@fumic, 1)
  expect_equal(x@rb, 1)
  expect_equal(x@fa, 1)
  expect_equal(x@fg, 1)
  expect_equal(x@ka, 0.1)
  expect_identical(x@solubility, Inf)
  expect_length(x@source, 0)
})


test_that("perpetrator constructor validates argument types and ranges", {
  expect_error(make_test_perpetrator(name = 1), "invalid class|must be a character value")
  expect_error(make_test_perpetrator(oral = 1), "invalid class|must be a logical value")
  expect_error(make_test_perpetrator(mw = -1), "must be positive")
  expect_error(make_test_perpetrator(dose = -1), "must be positive")
  expect_error(make_test_perpetrator(imaxss = -1), "must be positive")
  expect_error(make_test_perpetrator(fu = 1.1), "must be between 0 and 1")
  expect_error(make_test_perpetrator(fumic = 1.1), "must be between 0 and 1")
  expect_error(make_test_perpetrator(rb = 1.1), "must be between 0 and 1")
  expect_error(make_test_perpetrator(fa = 1.1), "must be between 0 and 1")
  expect_error(make_test_perpetrator(fg = 1.1), "must be between 0 and 1")
  expect_error(make_test_perpetrator(ka = -0.1), "must be positive")
  expect_error(make_test_perpetrator(solubility = -1), "must be positive")
  expect_error(make_test_perpetrator(source = 1), "invalid class|must be a character value")
})


test_that("perpetrator constructor validates source names", {
  expect_error(
    make_test_perpetrator(source = c(bad_name = "some source")),
    "unxpected source"
  )
})


test_that("perpetrator concentration functions reject non-perpetrator objects", {
  expect_error(igut("not-a-perpetrator"), "object must be a perpetrotor object")
  expect_error(imaxssu("not-a-perpetrator"), "object must be a perpetrotor object")
  expect_error(imaxinletu("not-a-perpetrator"), "object must be a perpetrotor object")
  expect_error(imaxintest("not-a-perpetrator"), "object must be a perpetrotor object")
  expect_error(key_conc("not-a-perpetrator"), "object must be a perpetrotor object")
})


test_that("igut returns expected concentrations and handles non-oral branch", {
  oral_x <- make_test_perpetrator()
  iv_x <- make_test_perpetrator(oral = FALSE)
  expected_mass <- 450 / 250 * 1e6
  expected_molar <- expected_mass / 492.6

  expect_equal(igut(oral_x, molar = FALSE), expected_mass)
  expect_equal(igut(oral_x, molar = TRUE), expected_molar)
  expect_equal(igut(iv_x, molar = FALSE), 0)
  expect_equal(igut(iv_x, molar = TRUE), 0)
})


test_that("igut is limited by solubility and warns", {
  x <- make_test_perpetrator(solubility = 100)

  expect_warning(
    igut_value <- igut(x, molar = FALSE),
    "Igut is limited by solubility"
  )
  expect_equal(igut_value, 100 * 1000)
})


test_that("imaxssu returns expected mass and molar concentrations", {
  x <- make_test_perpetrator()
  expected_mass <- 3530 * 0.023
  expected_molar <- expected_mass / 492.6

  expect_equal(imaxssu(x, molar = FALSE), expected_mass)
  expect_equal(imaxssu(x, molar = TRUE), expected_molar)
})


test_that("portal_term includes oral input and qh scaling", {
  x <- make_test_perpetrator()
  qh <- 1.616
  expected <- 450 * 0.81 * 1 * 0.00267 / qh / 1 * 1000

  expect_equal(portal_term(x, qh = qh), expected)
  expect_equal(portal_term(make_test_perpetrator(oral = FALSE), qh = qh), 0)
})


test_that("imaxinletu includes portal term for oral and omits for non-oral", {
  oral_x <- make_test_perpetrator()
  iv_x <- make_test_perpetrator(oral = FALSE)
  qh <- 1.616
  portal <- 450 * 0.81 * 1 * 0.00267 / qh / 1 * 1000

  oral_mass <- (3530 + portal) * 0.023
  iv_mass <- 3530 * 0.023

  expect_equal(imaxinletu(oral_x, qh = qh, molar = FALSE), oral_mass)
  expect_equal(imaxinletu(oral_x, qh = qh, molar = TRUE), oral_mass / 492.6)
  expect_equal(imaxinletu(iv_x, qh = qh, molar = FALSE), iv_mass)
  expect_equal(imaxinletu(iv_x, qh = qh, molar = TRUE), iv_mass / 492.6)
})


test_that("imaxintest returns expected values for oral and non-oral compounds", {
  oral_x <- make_test_perpetrator()
  iv_x <- make_test_perpetrator(oral = FALSE)
  qent <- 18 / 60

  oral_mass <- 450 * 0.81 * 0.00267 / qent * 1000
  iv_nonmolar <- imaxssu(iv_x, molar = TRUE)
  iv_molar <- iv_nonmolar / iv_x@mw

  expect_equal(imaxintest(oral_x, qent = qent, molar = FALSE), oral_mass)
  expect_equal(imaxintest(oral_x, qent = qent, molar = TRUE), oral_mass / 492.6)
  expect_equal(imaxintest(iv_x, qent = qent, molar = FALSE), iv_nonmolar)
  expect_equal(imaxintest(iv_x, qent = qent, molar = TRUE), iv_molar)
})


test_that("key_conc returns a formatted table with key concentration labels", {
  x <- make_test_perpetrator()
  output <- capture.output(key_conc(x, round = 2))

  expect_true(any(grepl("Key perpetrator concentrations for examplinib", output)))
  expect_true(any(grepl("\\$I_\\{gut\\}\\$", output)))
  expect_true(any(grepl("\\$I_\\{max,ss,u\\}\\$", output)))
  expect_true(any(grepl("\\$I_\\{max,inlet,u\\}\\$", output)))
  expect_true(any(grepl("\\$I_\\{max,intestinal\\}\\$", output)))
})


test_that("show and print methods render expected sections", {
  x <- make_test_perpetrator()
  show_text <- capture.output(show(x))
  print_text <- capture.output(print(x))

  expect_true(any(grepl("DDI perpetrator object", show_text)))
  expect_true(any(grepl("examplinib", show_text)))
  expect_true(any(grepl("Perpetrator compound parameters for examplinib", print_text)))
  expect_true(any(grepl("\\$C_\\{max,ss\\}\\$", print_text)))
  expect_true(any(grepl("\\$f_\\{u,mic\\}\\$", print_text)))
})

