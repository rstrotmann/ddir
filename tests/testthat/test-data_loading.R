
test_that("compound imports from textwithout error", {
  compounds <- read_perpetrators(textConnection(examplinib_compounds_string))
  expect_equal(length(compounds), 2)
})


test_that("new_perpetrator - compatibility with old 'type' field", {
  test_p <- tribble(
    ~param,     ~value,      ~source,
    "name",     "test",      "",
    "type",     "parent",    "",
    "mw",       "492.6",     "",
    "dose",     "450",       "clinical dose",
    "imaxss",   "3530",      "study 001",
    "fu",       "0.023",     "study 002",
    "fumic",    "1",         "default",
    "rb",       "1",         "study 003",
    "fa",       "0.81",      "study 003",
    "fg",       "1",         "default",
    "ka",       "0.00267",   "unknown")

  test_m <- tribble(
    ~param,     ~value,      ~source,
    "name",     "test",      "",
    "type",     "metabolite","",
    "mw",       "492.6",     "",
    "dose",     "450",       "clinical dose",
    "imaxss",   "3530",      "study 001",
    "fu",       "0.023",     "study 002",
    "fumic",    "1",         "default",
    "rb",       "1",         "study 003",
    "fa",       "0.81",      "study 003",
    "fg",       "1",         "default",
    "ka",       "0.00267",   "unknown")

  new_perpetrator(test_p)
  new_perpetrator(test_m)
})


test_that("read CYP inhibitor data without error", {
  test <- "
    examplinib, CYP1A2,  NA, \n
    examplinib, CYP2B6,  NA,\n
    examplinib, CYP2C8,  11,   study 001\n
    examplinib, CYP2C9,  13.5, study 001\n
    examplinib, CYP2C19, 15,   study 001\n
    examplinib, CYP2D6,  NA,\n
    examplinib, CYP3A4,  12.5, study 001\n
    \n
    # METABOLITE\n
    M1,         CYP2C9,  4.4,  study 002\n"

  read_inhibitor_data(textConnection(test))
})
