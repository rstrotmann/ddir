test_that("compound loading", {
  cpds <- read_perpetrators(test_path("fixtures", "compounds.csv"))
})



test_that("solubility limitation", {
  test_compound <- examplinib_compounds[[1]] %>%
    as.data.frame() %>%
    add_row(name="examplinib", param="solubility", value="1000",
            source="test")
  row.names(test_compound) <- test_compound$param
  test_compound <- new_perpetrator(test_compound)

  expect_true(is_igut_solubility_limited(test_compound))
})
