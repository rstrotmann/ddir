test_that("compounds load without error", {
  expect_no_error(
    read_perpetrators(test_path("fixtures", "examplinib_compounds.csv")))
})


test_that("new_perpetrator", {
  temp <- tribble(
    ~param,    ~value,     ~source,
    "name",     "test",      "",
    "oral",     "TRUE",      "",
    "mw",       "492.6",     "",
    "dose",     "450",       "clinical dose",
    "imaxss",   "3530",      "study 001",
    "fu",       "0.023",     "study 002",
    "fumic",    "1",         "default",
    "rb",       "1",         "study 003",
    "fa",       "0.81",      "study 003",
    "fg",       "1",         "default",
    "ka",       "0.00267",   "unknown")
  new_perpetrator(temp)
})


test_that("solubility limitation works", {
  test_compound <- examplinib_parent %>% as.data.frame()
  # test_compound <- examplinib_parent
  expect_false(is_igut_solubility_limited(test_compound))

  test_compound[which(test_compound$param=="solubility"), "value"] <- 100
  test_compound <- new_perpetrator(test_compound)

  expect_true(is_igut_solubility_limited(test_compound))
})


test_that("key concentrations", {
  expect_no_error(key_concentrations(examplinib_parent))
})
