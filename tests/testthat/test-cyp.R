test_that("cyp inhibitor data loads", {
  cyp_inh <- load_dmpk_data(test_path("fixtures", "cyp-inhibition.csv"))
})

test_that("cyp inhibitor data loads", {
  cyp_ind <- load_dmpk_data(test_path("fixtures", "cyp-induction.csv"))
})
