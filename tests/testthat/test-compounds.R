test_that("compound loads", {
  cpds <- load_perpetrators(test_path("fixtures", "compounds.csv"))
  print(cpds)
})
