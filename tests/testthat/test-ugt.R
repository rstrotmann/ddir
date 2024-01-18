test_that("ugt inhibitor data loads without error", {
  expect_no_error(
    read_inhibitor_data(test_path("fixtures",
                                           "examplinib_ugt_inhibition.csv")))
})
