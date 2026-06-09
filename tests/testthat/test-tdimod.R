emax_nls_start <- ddir:::emax_nls_start


make_tdi_activity_data <- function(
    kinact,
    kI,
    conc_levels = c(10, 50, 100),
    times = c(0, 10, 20),
    noise_sd = 0.5
) {
  rows <- lapply(conc_levels, function(conc) {
    kobs <- kinact * conc / (kI + conc)
    act <- 100 * exp(-kobs * times)
    if (noise_sd > 0)
      act <- pmax(act + stats::rnorm(length(times), sd = noise_sd), 1)
    data.frame(
      TIME = times,
      CONC = conc,
      ACT  = act
    )
  })
  do.call(rbind, rows)
}


test_that("tdimod returns the expected list structure on valid input", {
  set.seed(1)
  x <- make_tdi_activity_data(kinact = 0.04, kI = 50)

  res <- suppressWarnings(tdimod(x))

  expect_type(res, "list")
  expect_named(
    res,
    c("data", "kobs", "kobs_plot", "tdi_param", "tdi_plot")
  )
  expect_s3_class(res$data, "data.frame")
  expect_s3_class(res$kobs, "data.frame")
  expect_s3_class(res$kobs_plot, "ggplot")
  expect_s3_class(res$tdi_plot, "ggplot")
  expect_s3_class(res$tdi_param, "data.frame")
})


test_that("tdimod rejects non-data.frame input", {
  expect_error(
    tdimod(list(TIME = 0, CONC = 1, ACT = 100)),
    "x must be a data frame"
  )
})


test_that("tdimod rejects missing required columns", {
  x <- tibble::tribble(
    ~TIME, ~CONC,
        0,     1,
       10,     1
  )

  expect_error(
    tdimod(x),
    "Missing column\\(s\\) in input: ACT"
  )
})


test_that("tdimod accepts extra columns in input", {
  set.seed(1)
  x <- make_tdi_activity_data(kinact = 0.04, kI = 50)
  x$batch <- "A"

  expect_no_error(suppressWarnings(tdimod(x)))
})


test_that("tdimod validates target_name", {
  set.seed(1)
  x <- make_tdi_activity_data(kinact = 0.04, kI = 50)
  expect_no_error(suppressWarnings(tdimod(x, target_name = NULL)))
  expect_no_error(suppressWarnings(tdimod(x, target_name = "CYP3A4")))
  expect_error(
    tdimod(x, target_name = 1),
    "target_name must be a character value"
  )
})


test_that("tdimod rejects non-positive ACT values", {
  x <- tibble::tribble(
    ~TIME, ~CONC, ~ACT,
        0,    10,   100,
       10,    10,     0,
       20,    10,    25,
        0,    50,   100,
       10,    50,    50,
       20,    50,    25,
        0,   100,   100,
       10,   100,    50,
       20,   100,    25
  )

  expect_error(
    tdimod(x),
    "NA/NaN/Inf in 'y'"
  )
})


test_that("tdimod derives kobs as the negative slope of log activity over time", {
  x <- tibble::tribble(
    ~TIME, ~CONC, ~ACT,
        0,    10,   100,
       10,    10,    50,
       20,    10,    25,
        0,    50,   100,
       10,    50,    50,
       20,    50,    25,
        0,   100,   100,
       10,   100,    50,
       20,   100,    25
  )

  res <- suppressWarnings(tdimod(x))

  expect_named(res$kobs, c("CONC", "kobs", "std.error", "p.value"))
  expect_length(unique(res$kobs$CONC), 3)
  expect_equal(
    res$kobs$kobs[res$kobs$CONC == 10],
    log(2) / 10,
    tolerance = 1e-4
  )
})


test_that("tdimod produces identical kobs for sorted and shuffled input", {
  set.seed(1)
  x <- make_tdi_activity_data(kinact = 0.04, kI = 50)

  res_sorted <- suppressWarnings(tdimod(x[order(x$CONC, x$TIME), ]))
  res_shuffled <- suppressWarnings(tdimod(x[sample(nrow(x)), ]))

  expect_equal(res_sorted$kobs, res_shuffled$kobs)
})


test_that("tdimod fits synthetic Emax data to recover Kinact and KI", {
  true_kinact <- 0.04
  true_kI <- 50
  set.seed(2)
  x <- make_tdi_activity_data(kinact = true_kinact, kI = true_kI)

  res <- suppressWarnings(tdimod(x))

  kinact_est <- res$tdi_param$estimate[res$tdi_param$term == "kinact"]
  kI_est <- res$tdi_param$estimate[res$tdi_param$term == "kI"]

  expect_equal(kinact_est, true_kinact, tolerance = 0.01)
  expect_equal(kI_est, true_kI, tolerance = 8)
})


test_that("tdimod fits bundled examplinib in vitro TDI data", {
  res <- tdimod(examplinib_in_vitro_tdi)

  expect_setequal(res$tdi_param$term, c("kinact", "kI"))
  expect_equal(
    res$tdi_param$estimate[res$tdi_param$term == "kinact"],
    0.0412,
    tolerance = 0.001
  )
  expect_equal(
    res$tdi_param$estimate[res$tdi_param$term == "kI"],
    43.6,
    tolerance = 1
  )
  expect_length(unique(res$kobs$CONC), 7)
})


test_that("tdimod skips Emax fit with fewer than 3 concentration levels", {
  set.seed(3)
  x <- make_tdi_activity_data(
    kinact      = 0.04,
    kI          = 50,
    conc_levels = c(10, 50)
  )

  expect_warning(
    res <- tdimod(x),
    "At least 3 concentration levels are needed to fit Kinact and KI"
  )
  expect_null(res$tdi_param)
  expect_s3_class(res$tdi_plot, "ggplot")
  expect_length(ggplot2::ggplot_build(res$tdi_plot)$data, 2)
})


test_that("tdimod warns and skips Emax fit when kobs cannot be estimated", {
  x <- tibble::tribble(
    ~TIME, ~CONC, ~ACT,
        0,    10,   100,
        0,    50,   100,
        0,   100,   100
  )

  expect_warning(
    res <- tdimod(x),
    "Emax fit of kobs over CONC failed"
  )
  expect_null(res$tdi_param)
  expect_true(all(is.na(res$kobs$kobs)))
})


test_that("tdimod adds Emax fit line to tdi_plot when fit succeeds", {
  set.seed(4)
  x <- make_tdi_activity_data(kinact = 0.04, kI = 50)

  res <- suppressWarnings(tdimod(x))

  expect_length(ggplot2::ggplot_build(res$tdi_plot)$data, 3)
})


test_that("emax_nls_start derives starting values from linearization", {
  kobs_df <- tibble::tribble(
    ~CONC, ~kobs,
       20,  0.017,
       50,  0.020,
      100,  0.029
  )

  start <- emax_nls_start(kobs_df)

  expect_type(start, "list")
  expect_named(start, c("kinact", "kI"))
  expect_gt(start$kinact, 0)
  expect_gt(start$kI, 0)
  expect_false(identical(start, list(kinact = 0.03, kI = 10)))
})


test_that("emax_nls_start falls back when too few positive kobs values exist", {
  kobs_df <- tibble::tribble(
    ~CONC, ~kobs,
        1, -0.01,
       10, -0.02
  )

  start <- emax_nls_start(kobs_df)

  expect_equal(start, list(kinact = 0.03, kI = 10))
})


test_that("emax_nls_start falls back when linearization yields non-positive slope", {
  kobs_df <- tibble::tribble(
    ~CONC, ~kobs,
       10,  0.10,
       50,  0.05,
      100,  0.02
  )

  start <- emax_nls_start(kobs_df)

  expect_equal(start, list(kinact = 0.03, kI = 10))
})
