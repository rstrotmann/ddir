portal_term <- ddir:::portal_term


make_test_precipitant <- function(...) {
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
  do.call(precipitant, utils::modifyList(defaults, list(...)))
}


with_knitr <- function(expr) {
  old <- options(knitr.in.progress = TRUE)
  on.exit(options(old), add = TRUE)
  force(expr)
}


expect_table_contains_signif <- function(text, x, digits) {
  v <- signif(x, digits)
  sci <- formatC(v, format = "e", digits = max(digits - 1L, 0L))
  expect_true(
    grepl(as.character(v), text, fixed = TRUE) ||
      grepl(sci, text, fixed = TRUE)
  )
}


test_that("precipitant constructor populates all fields", {
  x <- make_test_precipitant()

  expect_s3_class(x, "precipitant")
  expect_identical(x$name, "examplinib")
  expect_identical(x$oral, TRUE)
  expect_equal(x$mw, 492.6)
  expect_equal(x$dose, 450)
  expect_equal(x$imaxss, 3530)
  expect_equal(x$fu, 0.023)
  expect_equal(x$fumic, 1)
  expect_equal(x$rb, 1)
  expect_equal(x$fa, 0.81)
  expect_equal(x$fg, 1)
  expect_equal(x$ka, 0.00267)
  expect_identical(x$solubility, Inf)
  expect_identical(
    unname(x$source[c("dose", "imaxss", "fu")]),
    c("clinical dose", "study 001", "study 002")
  )
})


test_that("precipitant constructor applies defaults", {
  x <- precipitant(
    name = "default-test",
    oral = TRUE,
    mw = 500,
    dose = 100,
    imaxss = 1000
  )

  expect_s3_class(x, "precipitant")
  expect_equal(x$fu, 1)
  expect_equal(x$fumic, 1)
  expect_equal(x$rb, 1)
  expect_equal(x$fa, 1)
  expect_equal(x$fg, 1)
  expect_equal(x$ka, 0.1)
  expect_identical(x$solubility, Inf)
  expect_length(x$source, 0)
})


test_that("empty precipitant() constructor is valid", {
  expect_no_error(x <- precipitant())
  expect_s3_class(x, "precipitant")
  expect_identical(x$name, "")
})


test_that("precipitant constructor validates argument types and ranges", {
  expect_error(make_test_precipitant(name = 1), "must be a character value")
  expect_error(
    make_test_precipitant(oral = c(TRUE, FALSE)),
    "must be a single value"
  )
  expect_error(make_test_precipitant(mw = -1), "must be positive")
  expect_error(make_test_precipitant(dose = -1), "must be positive")
  expect_error(make_test_precipitant(imaxss = -1), "must be positive")
  expect_error(make_test_precipitant(fu = 1.1), "must be between 0 and 1")
  expect_error(make_test_precipitant(fumic = 1.1), "must be between 0 and 1")
  expect_error(make_test_precipitant(fa = 1.1), "must be between 0 and 1")
  expect_error(make_test_precipitant(fg = 1.1), "must be between 0 and 1")
  expect_error(make_test_precipitant(ka = -0.1), "must be positive")
  expect_error(make_test_precipitant(solubility = -1), "must be positive")
  expect_error(make_test_precipitant(source = 1), "must be a character value")
})


test_that("precipitant constructor validates source names", {
  expect_error(
    make_test_precipitant(source = c(bad_name = "some source")),
    "unknown parameter name"
  )
  expect_error(
    make_test_precipitant(source = c(foo = "a", bar = "b")),
    "unknown parameter name"
  )
  expect_error(
    make_test_precipitant(source = "study note"),
    "must be a named vector"
  )
})


test_that("precipitant concentration functions reject non-precipitant objects", {
  expect_error(igut("not-a-precipitant"), "must be a precipitant object")
  expect_error(imaxssu("not-a-precipitant"), "must be a precipitant object")
  expect_error(imaxinletu("not-a-precipitant"), "must be a precipitant object")
  expect_error(imaxintest("not-a-precipitant"), "must be a precipitant object")
  expect_error(key_conc_table("not-a-precipitant"), "must be a precipitant object")
  expect_error(portal_term("not-a-precipitant"), "must be a precipitant object")
})


test_that("igut returns expected concentrations and handles non-oral branch", {
  oral_x <- make_test_precipitant()
  iv_x <- make_test_precipitant(oral = FALSE)
  expected_mass <- 450 / 250 * 1e6
  expected_molar <- expected_mass / 492.6

  expect_equal(igut(oral_x, molar = FALSE), expected_mass)
  expect_equal(igut(oral_x, molar = TRUE), expected_molar)
  expect_equal(igut(iv_x, molar = FALSE), 0)
  expect_equal(igut(iv_x, molar = TRUE), 0)
})


test_that("igut is limited by solubility and warns", {
  x <- make_test_precipitant(solubility = 100)

  expect_warning(
    igut_value <- igut(x, molar = FALSE),
    "Igut is limited by solubility"
  )
  expect_equal(igut_value, 100 * 1000)
})


test_that("imaxssu returns expected mass and molar concentrations", {
  x <- make_test_precipitant()
  expected_mass <- 3530 * 0.023
  expected_molar <- expected_mass / 492.6

  expect_equal(imaxssu(x, molar = FALSE), expected_mass)
  expect_equal(imaxssu(x, molar = TRUE), expected_molar)
})


test_that("portal_term includes oral input and qh scaling", {
  x <- make_test_precipitant()
  qh <- 1.616
  expected <- 450 * 0.81 * 1 * 0.00267 / qh / 1 * 1000

  expect_equal(portal_term(x, qh = qh), expected)
  expect_equal(portal_term(make_test_precipitant(oral = FALSE), qh = qh), 0)
})


test_that("imaxinletu includes portal term for oral and omits for non-oral", {
  oral_x <- make_test_precipitant()
  iv_x <- make_test_precipitant(oral = FALSE)
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
  oral_x <- make_test_precipitant()
  iv_x <- make_test_precipitant(oral = FALSE)
  qent <- 18 / 60

  oral_mass <- 450 * 0.81 * 0.00267 / qent * 1000
  iv_mass <- 3530 * 0.023

  expect_equal(imaxintest(oral_x, qent = qent, molar = FALSE), oral_mass)
  expect_equal(imaxintest(oral_x, qent = qent, molar = TRUE), oral_mass / 492.6)
  expect_equal(imaxintest(iv_x, qent = qent, molar = FALSE), iv_mass)
  expect_equal(imaxintest(iv_x, qent = qent, molar = TRUE), iv_mass / 492.6)
  expect_equal(imaxintest(iv_x, qent = qent, molar = FALSE), imaxssu(iv_x, molar = FALSE))
  expect_equal(imaxintest(iv_x, qent = qent, molar = TRUE), imaxssu(iv_x, molar = TRUE))
})


test_that("key_conc_table returns a formatted table with key concentration labels", {
  x <- make_test_precipitant()
  output <- capture.output(key_conc_table(x, signif = 3))

  expect_true(any(grepl("Key precipitant concentrations for examplinib", output)))
  expect_true(any(grepl("\\$I_\\{gut\\}\\$", output)))
  expect_true(any(grepl("\\$I_\\{max,ss,u\\}\\$", output)))
  expect_true(any(grepl("\\$I_\\{max,inlet,u\\}\\$", output)))
  expect_true(any(grepl("\\$I_\\{max,intestinal\\}\\$", output)))
})


test_that("print renders a console table with source notes", {
  x <- make_test_precipitant()
  text <- paste(capture.output(print.precipitant(x)), collapse = "\n")

  expect_match(text, "DDI precipitant")
  expect_match(text, "examplinib")
  expect_match(text, "492.6")
  expect_match(text, "0.023")
  expect_match(text, "clinical dose")
  expect_match(text, "study 001")
  expect_match(text, "study 002")
})


test_that("print without source notes does not add empty source labels", {
  x <- make_test_precipitant(source = character(0))
  text <- paste(capture.output(print.precipitant(x)), collapse = "\n")

  expect_false(grepl("clinical dose", text))
  expect_false(grepl("\\(\\)", text))
})


test_that("print in knitr returns a kable of oral parameters", {
  x <- make_test_precipitant()
  text <- with_knitr(paste(capture.output(print.precipitant(x)), collapse = "\n"))

  expect_match(text, "Precipitant compound parameters for examplinib")
  expect_match(text, "\\$F_a\\$")
  expect_match(text, "\\$F_g\\$")
  expect_match(text, "\\$k_a\\$")
  expect_match(text, "\\$R_B\\$")
  expect_match(text, "\\$C_\\{max,ss\\}\\$")
  expect_match(text, "\\$f_\\{u,mic\\}\\$")
  expect_match(text, "clinical dose")
  expect_match(text, "study 001")
})


test_that("print in knitr omits fa, fg, and ka for non-oral compounds", {
  x <- make_test_precipitant(oral = FALSE)
  text <- with_knitr(paste(capture.output(print.precipitant(x)), collapse = "\n"))

  expect_false(grepl("\\$F_a\\$", text))
  expect_false(grepl("\\$F_g\\$", text))
  expect_false(grepl("\\$k_a\\$", text))
  expect_match(text, "\\$MW\\$")
  expect_match(text, "\\$C_\\{max,ss\\}\\$")
})


test_that("precipitant allows empty name and NA dose, fa, fg, and ka", {
  x <- make_test_precipitant(
    name = "",
    dose = NA_real_,
    fa = NA_real_,
    fg = NA_real_,
    ka = NA_real_
  )

  expect_identical(x$name, "")
  expect_true(is.na(x$dose))
  expect_true(is.na(x$fa))
  expect_true(is.na(x$fg))
  expect_true(is.na(x$ka))
})


test_that("precipitant allows fraction boundaries, rb above 1, and zero numeric fields", {
  x <- make_test_precipitant(
    mw = 0,
    imaxss = 0,
    fu = 0,
    fumic = 1,
    rb = 1.4,
    fa = 0,
    fg = 1,
    ka = 0,
    solubility = 0
  )

  expect_equal(x$fu, 0)
  expect_equal(x$fumic, 1)
  expect_equal(x$rb, 1.4)
  expect_equal(x$fa, 0)
  expect_equal(x$fg, 1)
  expect_equal(x$mw, 0)
  expect_equal(x$imaxss, 0)
  expect_equal(x$ka, 0)
  expect_equal(x$solubility, 0)
})


test_that("precipitant rejects NA in fields that do not allow it", {
  expect_error(make_test_precipitant(name = NA_character_), "must not contain NA")
  expect_error(make_test_precipitant(oral = NA), "must not contain NA")
  expect_error(make_test_precipitant(mw = NA_real_), "must not contain NA")
  expect_error(make_test_precipitant(imaxss = NA_real_), "must not contain NA")
  expect_error(make_test_precipitant(fu = NA_real_), "must not contain NA")
  expect_error(make_test_precipitant(fumic = NA_real_), "must not contain NA")
  expect_error(make_test_precipitant(rb = NA_real_), "must not contain NA")
  expect_error(make_test_precipitant(solubility = NA_real_), "must not contain NA")
})


test_that("precipitant rejects a length-2 name", {
  expect_error(
    make_test_precipitant(name = c("a", "b")),
    "must be a single value"
  )
})


test_that("precipitant accepts source names that match field names", {
  source <- c(
    name = "label",
    oral = "route",
    mw = "certificate",
    dose = "clinical dose",
    imaxss = "study 001",
    fu = "study 002",
    fumic = "default",
    rb = "study 003",
    fa = "study 003",
    fg = "default",
    ka = "unknown",
    solubility = "default",
    source = "meta"
  )
  x <- make_test_precipitant(source = source)

  expect_identical(x$source, source)
})


test_that("igut defaults to mass concentration", {
  x <- make_test_precipitant()

  expect_equal(igut(x), igut(x, molar = FALSE))
  expect_false(isTRUE(all.equal(igut(x), igut(x, molar = TRUE))))
})


test_that("igut does not warn when solubility is Inf or equal to Igut", {
  unlimited <- make_test_precipitant(solubility = Inf)
  at_limit <- make_test_precipitant(solubility = (450 / 250 * 1e6) / 1000)

  expect_no_warning(igut(unlimited, molar = FALSE))
  expect_equal(igut(unlimited, molar = FALSE), 450 / 250 * 1e6)
  expect_no_warning(igut(at_limit, molar = FALSE))
  expect_equal(igut(at_limit, molar = FALSE), 450 / 250 * 1e6)
})


test_that("igut molar output is also limited by solubility", {
  x <- make_test_precipitant(solubility = 100)

  expect_warning(
    molar_igut <- igut(x, molar = TRUE),
    "Igut is limited by solubility"
  )
  expect_equal(molar_igut, 100 * 1000 / 492.6)
})


test_that("igut for IV is zero regardless of dose and solubility", {
  x <- make_test_precipitant(oral = FALSE, dose = 1000, solubility = 1)

  expect_no_warning(expect_equal(igut(x, molar = FALSE), 0))
  expect_equal(igut(x, molar = TRUE), 0)
})


test_that("imaxssu defaults to molar concentration", {
  x <- make_test_precipitant()

  expect_equal(imaxssu(x), imaxssu(x, molar = TRUE))
  expect_equal(imaxssu(x, molar = FALSE), 3530 * 0.023)
  expect_equal(imaxssu(x, molar = TRUE) * x$mw, imaxssu(x, molar = FALSE))
})


test_that("imaxssu is zero when fu or imaxss is zero", {
  expect_equal(imaxssu(make_test_precipitant(fu = 0), molar = FALSE), 0)
  expect_equal(imaxssu(make_test_precipitant(imaxss = 0), molar = FALSE), 0)
})


test_that("portal_term uses default qh and scales with rb", {
  x <- make_test_precipitant()
  expected_default <- 450 * 0.81 * 1 * 0.00267 / 1.616 / 1 * 1000
  expected_rb <- 450 * 0.81 * 1 * 0.00267 / 1.616 / 1.4 * 1000

  expect_equal(portal_term(x), expected_default)
  expect_equal(portal_term(x, qh = 1.616), expected_default)
  expect_equal(
    portal_term(make_test_precipitant(rb = 1.4), qh = 1.616),
    expected_rb
  )
})


test_that("portal_term scales inversely with qh and includes fa and fg", {
  x <- make_test_precipitant()
  half_qh <- portal_term(x, qh = 0.808)
  default_qh <- portal_term(x, qh = 1.616)
  no_fa <- portal_term(make_test_precipitant(fa = 0), qh = 1.616)
  half_fg <- portal_term(make_test_precipitant(fg = 0.5), qh = 1.616)

  expect_equal(half_qh, 2 * default_qh)
  expect_equal(no_fa, 0)
  expect_equal(half_fg, 0.5 * default_qh)
})


test_that("imaxinletu defaults to molar and matches imaxssu for IV", {
  oral_x <- make_test_precipitant()
  iv_x <- make_test_precipitant(oral = FALSE)

  expect_equal(imaxinletu(oral_x), imaxinletu(oral_x, molar = TRUE))
  expect_equal(imaxinletu(iv_x, molar = FALSE), imaxssu(iv_x, molar = FALSE))
  expect_equal(imaxinletu(iv_x, molar = TRUE), imaxssu(iv_x, molar = TRUE))
  expect_gt(imaxinletu(oral_x, molar = FALSE), imaxssu(oral_x, molar = FALSE))
})


test_that("imaxinletu uses custom qh", {
  x <- make_test_precipitant()
  portal <- 450 * 0.81 * 1 * 0.00267 / 0.808 / 1 * 1000
  expected_mass <- (3530 + portal) * 0.023

  expect_equal(imaxinletu(x, qh = 0.808, molar = FALSE), expected_mass)
  expect_equal(imaxinletu(x, qh = 0.808, molar = TRUE), expected_mass / 492.6)
})


test_that("imaxintest defaults to molar and is independent of fu, rb, and fg", {
  x <- make_test_precipitant()
  qent <- 18 / 60
  expected_mass <- 450 * 0.81 * 0.00267 / qent * 1000

  expect_equal(imaxintest(x), imaxintest(x, molar = TRUE))
  expect_equal(imaxintest(x, qent = qent, molar = FALSE), expected_mass)
  expect_equal(
    imaxintest(make_test_precipitant(fu = 0.5), qent = qent, molar = FALSE),
    expected_mass
  )
  expect_equal(
    imaxintest(make_test_precipitant(rb = 1.4), qent = qent, molar = FALSE),
    expected_mass
  )
  expect_equal(
    imaxintest(make_test_precipitant(fg = 0.5), qent = qent, molar = FALSE),
    expected_mass
  )
})


test_that("imaxintest scales inversely with qent and is zero when fa is zero", {
  x <- make_test_precipitant()
  default <- imaxintest(x, qent = 18 / 60, molar = FALSE)
  half <- imaxintest(x, qent = 9 / 60, molar = FALSE)

  expect_equal(half, 2 * default)
  expect_equal(
    imaxintest(make_test_precipitant(fa = 0), molar = FALSE),
    0
  )
})


test_that("key_conc_table returns knitr_kable with matching helper values", {
  x <- make_test_precipitant()
  digits <- 3
  out <- key_conc_table(x, signif = digits)
  text <- paste(capture.output(out), collapse = "\n")

  expect_s3_class(out, "knitr_kable")
  expect_match(text, "value \\(\\$ng/ml\\$\\)")
  expect_match(text, "value \\(\\$\\\\mu M\\$\\)")
  expect_table_contains_signif(text, igut(x, molar = FALSE), digits)
  expect_table_contains_signif(text, igut(x, molar = TRUE), digits)
  expect_table_contains_signif(text, imaxssu(x, molar = FALSE), digits)
  expect_table_contains_signif(text, imaxssu(x, molar = TRUE), digits)
  expect_table_contains_signif(text, imaxinletu(x, molar = FALSE), digits)
  expect_table_contains_signif(text, imaxinletu(x, molar = TRUE), digits)
  expect_table_contains_signif(text, imaxintest(x, molar = FALSE), digits)
  expect_table_contains_signif(text, imaxintest(x, molar = TRUE), digits)
})


test_that("key_conc_table honors signif, qh, and qent", {
  x <- make_test_precipitant()
  coarse <- paste(capture.output(key_conc_table(x, signif = 1)), collapse = "\n")
  custom <- paste(
    capture.output(key_conc_table(x, signif = 3, qh = 0.808, qent = 9 / 60)),
    collapse = "\n"
  )

  expect_table_contains_signif(coarse, imaxssu(x, molar = FALSE), 1)
  expect_table_contains_signif(
    custom,
    imaxinletu(x, qh = 0.808, molar = FALSE),
    3
  )
  expect_table_contains_signif(
    custom,
    imaxintest(x, qent = 9 / 60, molar = FALSE),
    3
  )
})


test_that("key_conc_table uses IV branches for gut and intestinal terms", {
  x <- make_test_precipitant(oral = FALSE)
  digits <- 3
  text <- paste(capture.output(key_conc_table(x, signif = digits)), collapse = "\n")

  expect_table_contains_signif(text, igut(x, molar = FALSE), digits)
  expect_table_contains_signif(text, imaxintest(x, molar = FALSE), digits)
  expect_table_contains_signif(text, imaxssu(x, molar = FALSE), digits)
  expect_equal(igut(x, molar = FALSE), 0)
  expect_equal(imaxintest(x, molar = FALSE), imaxssu(x, molar = FALSE))
})


test_that("package examplinib precipitant object is valid", {
  expect_s3_class(examplinib, "precipitant")
  expect_true(examplinib$oral)
  expect_identical(examplinib$name, "examplinib")
  expect_no_error(igut(examplinib))
  expect_no_error(imaxssu(examplinib))
  expect_no_error(imaxinletu(examplinib))
  expect_no_error(imaxintest(examplinib))
  expect_no_error(key_conc_table(examplinib))
})
