test_that("multiplication works", {
  expect_equal(2 * 2, 4)
})


test_that("compound imports from text", {
  compounds <- read_perpetrators(textConnection(examplinib_compounds_string))
  expect_equal(length(compounds), 2)
})
