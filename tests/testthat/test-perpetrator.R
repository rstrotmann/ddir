make_test_perpetrator <- function(
    oral = TRUE,
    solubility = Inf,
    source = c(dose = "clinical dose", imaxss = "study 001", fu = "study 002")
) {
  perpetrator(
    name = "examplinib",
    oral = oral,
    mw = 492.6,
    dose = 450,
    imaxss = 3530,
    fu = 0.023,
    fumic = 1,
    rb = 1,
    fa = 0.81,
    fg = 1,
    ka = 0.00267,
    solubility = solubility,
    source = source
  )
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
  expect_identical(unname(x@source[c("dose", "imaxss", "fu")]),
                   c("clinical dose", "study 001", "study 002"))
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


test_that("perpetrator constructor validates source names", {
  expect_error(
    make_test_perpetrator(source = c(bad_name = "some source")),
    "unxpected source"
  )
})


test_that("igut returns zero for non-oral perpetrator", {
  x <- make_test_perpetrator(oral = FALSE)

  expect_equal(igut(x, molar = FALSE), 0)
  expect_equal(igut(x, molar = TRUE), 0)
})


test_that("igut is limited by solubility and warns", {
  x <- make_test_perpetrator(solubility = 100)

  expect_warning(
    igut_value <- igut(x, molar = FALSE),
    "Igut is limited by solubility"
  )
  expect_equal(igut_value, 100 * 1000)
})


test_that("igut returns expected mass and molar concentrations", {
  x <- make_test_perpetrator()
  expected_mass <- 450 / 250 * 1e6
  expected_molar <- expected_mass / 492.6

  expect_equal(igut(x, molar = FALSE), expected_mass)
  expect_equal(igut(x, molar = TRUE), expected_molar)
})


test_that("imaxssu returns expected concentrations", {
  x <- make_test_perpetrator()
  expected_mass <- 3530 * 0.023
  expected_molar <- expected_mass / 492.6

  expect_equal(imaxssu(x, molar = FALSE), expected_mass)
  expect_equal(imaxssu(x, molar = TRUE), expected_molar)
})


test_that("imaxinletu includes portal term for oral perpetrator", {
  x <- make_test_perpetrator()
  qh <- 1.616
  portal <- 450 * 0.81 * 1 * 0.00267 / qh / 1 * 1000
  expected_mass <- (3530 + portal) * 0.023
  expected_molar <- expected_mass / 492.6

  expect_equal(imaxinletu(x, qh = qh, molar = FALSE), expected_mass)
  expect_equal(imaxinletu(x, qh = qh, molar = TRUE), expected_molar)
})


test_that("imaxinletu excludes portal term for non-oral perpetrator", {
  x <- make_test_perpetrator(oral = FALSE)
  expected_mass <- 3530 * 0.023
  expected_molar <- expected_mass / 492.6

  expect_equal(imaxinletu(x, molar = FALSE), expected_mass)
  expect_equal(imaxinletu(x, molar = TRUE), expected_molar)
})


test_that("imaxintest returns expected concentrations", {
  x <- make_test_perpetrator()
  qent <- 18 / 60
  expected_mass <- 450 * 0.81 * 0.00267 / qent * 1000
  expected_molar <- expected_mass / 492.6

  expect_equal(imaxintest(x, qent = qent, molar = FALSE), expected_mass)
  expect_equal(imaxintest(x, qent = qent, molar = TRUE), expected_molar)
})


test_that("imaxintest uses imaxssu for non-oral perpetrator", {
  x <- make_test_perpetrator(oral = FALSE)

  # current implementation computes imaxssu() with default molar = TRUE first
  expected_nonmolar <- imaxssu(x, molar = TRUE)
  expected_molar <- expected_nonmolar / x@mw

  expect_equal(imaxintest(x, molar = FALSE), expected_nonmolar)
  expect_equal(imaxintest(x, molar = TRUE), expected_molar)
})


test_that("kc returns all key concentrations in expected order", {
  x <- make_test_perpetrator()
  out <- kc(x, molar = TRUE)

  expect_named(out, c("igut", "imaxssu", "imaxinletu", "imaxintest"))
  expect_equal(out[["igut"]], igut(x, molar = TRUE))
  expect_equal(out[["imaxssu"]], imaxssu(x, molar = TRUE))
  expect_equal(out[["imaxinletu"]], imaxinletu(x, molar = TRUE))
  expect_equal(out[["imaxintest"]], imaxintest(x, molar = TRUE))
})


test_that("prop returns a kable and omits oral-only fields for non-oral", {
  oral_x <- make_test_perpetrator(oral = TRUE)
  nonoral_x <- make_test_perpetrator(oral = FALSE)

  oral_tbl <- prop(oral_x)
  nonoral_tbl <- prop(nonoral_x)

  expect_s3_class(oral_tbl, "knitr_kable")
  expect_s3_class(nonoral_tbl, "knitr_kable")

  nonoral_text <- as.character(nonoral_tbl)
  expect_false(any(grepl("\\$F_a\\$", nonoral_text)))
  expect_false(any(grepl("\\$F_g\\$", nonoral_text)))
  expect_false(any(grepl("\\$k_a\\$", nonoral_text)))
})


test_that("conctbl returns rounded concentration table as kable", {
  x <- make_test_perpetrator()
  tbl <- conctbl(x, round = 3)

  expect_s3_class(tbl, "knitr_kable")
  expect_true(any(grepl("Key perpetrator concentrations for examplinib",
                        as.character(tbl))))
})


test_that("show and print methods run without error and emit content", {
  x <- make_test_perpetrator()

  show_text <- capture.output(show(x))
  print_text <- capture.output(print(x))

  expect_true(any(grepl("DDI perpetrator object", show_text)))
  expect_true(any(grepl("Perpetrator compound parameters for examplinib",
                        print_text)))
})
