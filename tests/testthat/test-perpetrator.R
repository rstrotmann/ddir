portal_term <- ddir:::portal_term


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
  expect_error(key_conc_table("not-a-perpetrator"), "object must be a perpetrotor object")
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
  iv_mass <- 3530 * 0.023

  expect_equal(imaxintest(oral_x, qent = qent, molar = FALSE), oral_mass)
  expect_equal(imaxintest(oral_x, qent = qent, molar = TRUE), oral_mass / 492.6)
  expect_equal(imaxintest(iv_x, qent = qent, molar = FALSE), iv_mass)
  expect_equal(imaxintest(iv_x, qent = qent, molar = TRUE), iv_mass / 492.6)
  expect_equal(imaxintest(iv_x, qent = qent, molar = FALSE), imaxssu(iv_x, molar = FALSE))
  expect_equal(imaxintest(iv_x, qent = qent, molar = TRUE), imaxssu(iv_x, molar = TRUE))
})


test_that("key_conc_table returns a formatted table with key concentration labels", {
  x <- make_test_perpetrator()
  output <- capture.output(key_conc_table(x, round = 2))

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

  expect_true(any(grepl("DDI perpetrator", show_text)))
  expect_true(any(grepl("examplinib", show_text)))
  expect_true(any(grepl("Perpetrator compound parameters for examplinib", print_text)))
  expect_true(any(grepl("\\$C_\\{max,ss\\}\\$", print_text)))
  expect_true(any(grepl("\\$f_\\{u,mic\\}\\$", print_text)))
})


test_that("validObject accepts a complete perpetrator", {
  expect_true(validObject(make_test_perpetrator()))
  expect_true(validObject(new("perpetrator")))
})


test_that("perpetrator allows empty name and NA dose, fa, fg, and ka", {
  x <- make_test_perpetrator(
    name = "",
    dose = NA_real_,
    fa = NA_real_,
    fg = NA_real_,
    ka = NA_real_
  )

  expect_identical(x@name, "")
  expect_true(is.na(x@dose))
  expect_true(is.na(x@fa))
  expect_true(is.na(x@fg))
  expect_true(is.na(x@ka))
  expect_true(validObject(x))
})


test_that("perpetrator allows fraction boundaries, rb above 1, and zero numeric slots", {
  x <- make_test_perpetrator(
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

  expect_equal(x@fu, 0)
  expect_equal(x@fumic, 1)
  expect_equal(x@rb, 1.4)
  expect_equal(x@fa, 0)
  expect_equal(x@fg, 1)
  expect_equal(x@mw, 0)
  expect_equal(x@imaxss, 0)
  expect_equal(x@ka, 0)
  expect_equal(x@solubility, 0)
})


test_that("perpetrator rejects NA in slots that do not allow it", {
  expect_error(make_test_perpetrator(name = NA_character_), "must not contain NA")
  expect_error(make_test_perpetrator(oral = NA), "must not contain NA")
  expect_error(make_test_perpetrator(mw = NA_real_), "must not contain NA")
  expect_error(make_test_perpetrator(imaxss = NA_real_), "must not contain NA")
  expect_error(make_test_perpetrator(fu = NA_real_), "must not contain NA")
  expect_error(make_test_perpetrator(fumic = NA_real_), "must not contain NA")
  expect_error(make_test_perpetrator(rb = NA_real_), "must not contain NA")
  expect_error(make_test_perpetrator(solubility = NA_real_), "must not contain NA")
})


test_that("perpetrator rejects length-2 name and unnamed extra source is allowed", {
  expect_error(
    make_test_perpetrator(name = c("a", "b")),
    "must be a single value"
  )

  unnamed <- make_test_perpetrator(source = "study note")
  expect_identical(unnamed@source, "study note")
})


test_that("perpetrator accepts source names that match slot names", {
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
  x <- make_test_perpetrator(source = source)

  expect_identical(x@source, source)
})


test_that("perpetrator rejects multiple unexpected source names", {
  expect_error(
    make_test_perpetrator(source = c(foo = "a", bar = "b")),
    "unxpected source"
  )
})


test_that("show includes source annotations, units, and all slots", {
  x <- make_test_perpetrator()
  show_text <- paste(capture.output(show(x)), collapse = "\n")

  expect_match(show_text, "g/mol")
  expect_match(show_text, "ng/ml")
  expect_match(show_text, "mg/ml")
  expect_match(show_text, "/min")
  expect_match(show_text, "\\(clinical dose\\)")
  expect_match(show_text, "\\(study 001\\)")
  expect_match(show_text, "\\(study 002\\)")
  expect_match(show_text, "492.6")
  expect_match(show_text, "0.023")
})


test_that("show without source notes does not add parentheses", {
  x <- make_test_perpetrator(source = character(0))
  show_text <- paste(capture.output(show(x)), collapse = "\n")

  expect_false(grepl("\\(clinical dose\\)", show_text))
  expect_false(grepl("\\(\\)", show_text))
})


test_that("print returns a knitr_kable table of oral parameters", {
  x <- make_test_perpetrator()
  out <- print(x)
  text <- paste(capture.output(print(x)), collapse = "\n")

  expect_s3_class(out, "knitr_kable")
  expect_match(text, "\\$F_a\\$")
  expect_match(text, "\\$F_g\\$")
  expect_match(text, "\\$k_a\\$")
  expect_match(text, "\\$R_B\\$")
  expect_match(text, "mg/l")
  expect_match(text, "clinical dose")
  expect_match(text, "study 001")
})


test_that("print omits fa, fg, and ka for non-oral compounds", {
  x <- make_test_perpetrator(oral = FALSE)
  text <- paste(capture.output(print(x)), collapse = "\n")

  expect_false(grepl("\\$F_a\\$", text))
  expect_false(grepl("\\$F_g\\$", text))
  expect_false(grepl("\\$k_a\\$", text))
  expect_match(text, "\\$MW\\$")
  expect_match(text, "\\$C_\\{max,ss\\}\\$")
})


test_that("igut defaults to mass concentration", {
  x <- make_test_perpetrator()

  expect_equal(igut(x), igut(x, molar = FALSE))
  expect_false(isTRUE(all.equal(igut(x), igut(x, molar = TRUE))))
})


test_that("igut does not warn when solubility is Inf or equal to Igut", {
  unlimited <- make_test_perpetrator(solubility = Inf)
  at_limit <- make_test_perpetrator(solubility = (450 / 250 * 1e6) / 1000)

  expect_no_warning(igut(unlimited, molar = FALSE))
  expect_equal(igut(unlimited, molar = FALSE), 450 / 250 * 1e6)
  expect_no_warning(igut(at_limit, molar = FALSE))
  expect_equal(igut(at_limit, molar = FALSE), 450 / 250 * 1e6)
})


test_that("igut molar output is also limited by solubility", {
  x <- make_test_perpetrator(solubility = 100)

  expect_warning(
    molar_igut <- igut(x, molar = TRUE),
    "Igut is limited by solubility"
  )
  expect_equal(molar_igut, 100 * 1000 / 492.6)
})


test_that("igut for IV is zero regardless of dose and solubility", {
  x <- make_test_perpetrator(oral = FALSE, dose = 1000, solubility = 1)

  expect_no_warning(expect_equal(igut(x, molar = FALSE), 0))
  expect_equal(igut(x, molar = TRUE), 0)
})


test_that("imaxssu defaults to molar concentration", {
  x <- make_test_perpetrator()

  expect_equal(imaxssu(x), imaxssu(x, molar = TRUE))
  expect_equal(imaxssu(x, molar = FALSE), 3530 * 0.023)
  expect_equal(imaxssu(x, molar = TRUE) * x@mw, imaxssu(x, molar = FALSE))
})


test_that("imaxssu is zero when fu or imaxss is zero", {
  expect_equal(imaxssu(make_test_perpetrator(fu = 0), molar = FALSE), 0)
  expect_equal(imaxssu(make_test_perpetrator(imaxss = 0), molar = FALSE), 0)
})


test_that("portal_term uses default qh and scales with rb", {
  x <- make_test_perpetrator()
  expected_default <- 450 * 0.81 * 1 * 0.00267 / 1.616 / 1 * 1000
  expected_rb <- 450 * 0.81 * 1 * 0.00267 / 1.616 / 1.4 * 1000

  expect_equal(portal_term(x), expected_default)
  expect_equal(portal_term(x, qh = 1.616), expected_default)
  expect_equal(
    portal_term(make_test_perpetrator(rb = 1.4), qh = 1.616),
    expected_rb
  )
})


test_that("portal_term scales inversely with qh and includes fa and fg", {
  x <- make_test_perpetrator()
  half_qh <- portal_term(x, qh = 0.808)
  default_qh <- portal_term(x, qh = 1.616)
  no_fa <- portal_term(make_test_perpetrator(fa = 0), qh = 1.616)
  half_fg <- portal_term(make_test_perpetrator(fg = 0.5), qh = 1.616)

  expect_equal(half_qh, 2 * default_qh)
  expect_equal(no_fa, 0)
  expect_equal(half_fg, 0.5 * default_qh)
})


test_that("portal_term rejects non-perpetrator objects", {
  expect_error(portal_term("not-a-perpetrator"), "object must be a perpetrotor object")
})


test_that("imaxinletu defaults to molar and matches imaxssu for IV", {
  oral_x <- make_test_perpetrator()
  iv_x <- make_test_perpetrator(oral = FALSE)

  expect_equal(imaxinletu(oral_x), imaxinletu(oral_x, molar = TRUE))
  expect_equal(imaxinletu(iv_x, molar = FALSE), imaxssu(iv_x, molar = FALSE))
  expect_equal(imaxinletu(iv_x, molar = TRUE), imaxssu(iv_x, molar = TRUE))
  expect_gt(imaxinletu(oral_x, molar = FALSE), imaxssu(oral_x, molar = FALSE))
})


test_that("imaxinletu uses custom qh", {
  x <- make_test_perpetrator()
  portal <- 450 * 0.81 * 1 * 0.00267 / 0.808 / 1 * 1000
  expected_mass <- (3530 + portal) * 0.023

  expect_equal(imaxinletu(x, qh = 0.808, molar = FALSE), expected_mass)
  expect_equal(imaxinletu(x, qh = 0.808, molar = TRUE), expected_mass / 492.6)
})


test_that("imaxintest defaults to molar and is independent of fu, rb, and fg", {
  x <- make_test_perpetrator()
  qent <- 18 / 60
  expected_mass <- 450 * 0.81 * 0.00267 / qent * 1000

  expect_equal(imaxintest(x), imaxintest(x, molar = TRUE))
  expect_equal(imaxintest(x, qent = qent, molar = FALSE), expected_mass)
  expect_equal(
    imaxintest(make_test_perpetrator(fu = 0.5), qent = qent, molar = FALSE),
    expected_mass
  )
  expect_equal(
    imaxintest(make_test_perpetrator(rb = 1.4), qent = qent, molar = FALSE),
    expected_mass
  )
  expect_equal(
    imaxintest(make_test_perpetrator(fg = 0.5), qent = qent, molar = FALSE),
    expected_mass
  )
})


test_that("imaxintest scales inversely with qent and is zero when fa is zero", {
  x <- make_test_perpetrator()
  default <- imaxintest(x, qent = 18 / 60, molar = FALSE)
  half <- imaxintest(x, qent = 9 / 60, molar = FALSE)

  expect_equal(half, 2 * default)
  expect_equal(
    imaxintest(make_test_perpetrator(fa = 0), molar = FALSE),
    0
  )
})


test_that("key_conc_table returns knitr_kable with matching helper values", {
  x <- make_test_perpetrator()
  out <- key_conc_table(x, round = 2)
  text <- paste(capture.output(out), collapse = "\n")

  expect_s3_class(out, "knitr_kable")
  expect_match(text, "value \\(\\$ng/ml\\$\\)")
  expect_match(text, "value \\(\\$\\\\mu M\\$\\)")
  expect_match(text, as.character(round(igut(x, molar = FALSE), 2)))
  expect_match(text, as.character(round(igut(x, molar = TRUE), 2)))
  expect_match(text, as.character(round(imaxssu(x, molar = FALSE), 2)))
  expect_match(text, as.character(round(imaxssu(x, molar = TRUE), 2)))
  expect_match(text, as.character(round(imaxinletu(x, molar = FALSE), 2)))
  expect_match(text, as.character(round(imaxinletu(x, molar = TRUE), 2)))
  expect_match(text, as.character(round(imaxintest(x, molar = FALSE), 2)))
  expect_match(text, as.character(round(imaxintest(x, molar = TRUE), 2)))
})


test_that("key_conc_table honors round, qh, and qent", {
  x <- make_test_perpetrator()
  rounded <- paste(capture.output(key_conc_table(x, round = 0)), collapse = "\n")
  custom <- paste(
    capture.output(key_conc_table(x, round = 2, qh = 0.808, qent = 9 / 60)),
    collapse = "\n"
  )

  expect_match(rounded, as.character(round(imaxssu(x, molar = FALSE), 0)))
  expect_match(
    custom,
    as.character(round(imaxinletu(x, qh = 0.808, molar = FALSE), 2))
  )
  expect_match(
    custom,
    as.character(round(imaxintest(x, qent = 9 / 60, molar = FALSE), 2))
  )
})


test_that("key_conc_table uses IV branches for gut and intestinal terms", {
  x <- make_test_perpetrator(oral = FALSE)
  text <- paste(capture.output(key_conc_table(x, round = 2)), collapse = "\n")

  expect_match(text, as.character(round(igut(x, molar = FALSE), 2)))
  expect_match(text, as.character(round(imaxintest(x, molar = FALSE), 2)))
  expect_match(text, as.character(round(imaxssu(x, molar = FALSE), 2)))
  expect_equal(igut(x, molar = FALSE), 0)
  expect_equal(imaxintest(x, molar = FALSE), imaxssu(x, molar = FALSE))
})


test_that("package examplinib perpetrator objects are valid", {
  expect_s4_class(examplinib_parent, "perpetrator")
  expect_s4_class(examplinib_metabolite, "perpetrator")
  expect_true(validObject(examplinib_parent))
  expect_true(validObject(examplinib_metabolite))
  expect_true(examplinib_parent@oral)
  expect_false(examplinib_metabolite@oral)
  expect_true(is.na(examplinib_metabolite@dose))
  expect_no_error(igut(examplinib_parent))
  expect_no_error(imaxssu(examplinib_parent))
  expect_no_error(imaxinletu(examplinib_parent))
  expect_no_error(imaxintest(examplinib_parent))
  expect_no_error(key_conc_table(examplinib_parent))
  expect_equal(igut(examplinib_metabolite, molar = FALSE), 0)
  expect_equal(
    imaxintest(examplinib_metabolite, molar = FALSE),
    imaxssu(examplinib_metabolite, molar = FALSE)
  )
})

