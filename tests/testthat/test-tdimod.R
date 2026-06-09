emax_nls_start <- ddir:::emax_nls_start


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
