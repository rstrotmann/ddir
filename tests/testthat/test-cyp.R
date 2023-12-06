test_that("cyp inhibitor data loads", {
  cyp_inh <- read_dmpk_data(test_path("fixtures", "cyp-inhibition.csv"))
})

test_that("cyp inhibitor data loads", {
  cyp_ind <- read_cyp_inducer_data(test_path("fixtures", "cyp-induction.csv"))
})
