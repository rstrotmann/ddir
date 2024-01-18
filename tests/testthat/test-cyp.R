test_that("cyp inhibitor data loads without error", {
  expect_no_error(
    cyp_inh <- read_inhibitor_data(test_path("fixtures",
                                           "examplinib_cyp_inhibition.csv")))
})

test_that("cyp inhibitor data loads without error", {
  expect_no_error(
    cyp_ind <- read_inducer_data(test_path("fixtures",
                                         "examplinib_cyp_induction.csv")))
})
