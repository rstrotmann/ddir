# Transporter inhibition cases from Gomez-Mantilla et al., Clin Pharmacokinet 2023,
# Table 17 (baricitinib / OAT3 victims) and Posada et al., Clin Transl Sci 2017.
#
# Workflow: take the published MSM AUCR (ft = 0.44), back-calculate Iu/IC50,
# and check that transporter_inh_risk() recovers that ratio.


infer_i_over_ic50 <- function(aucr, ft = 0.44) {
  d <- 1 / aucr
  ft / (d - (1 - ft)) - 1
}


make_perp_from_Iu <- function(name, mw, Iu, fu = 1) {
  perpetrator(
    name       = name,
    oral       = FALSE,
    mw         = mw,
    dose       = 1,
    imaxss     = Iu * mw / fu,
    fu         = fu,
    fumic      = 1,
    rb         = 1,
    fa         = 1,
    fg         = 1,
    ka         = 0.1,
    solubility = Inf
  )
}


test_that("probenecid OAT3 I/IC50 recovered from Gomez-Mantilla 2023 AUCR 1.67 is a risk", {
  ic50 <- 4.4
  paper_aucr <- 1.67
  i_over_ic50 <- infer_i_over_ic50(paper_aucr, ft = 0.44)
  Iu <- i_over_ic50 * ic50
  perp <- make_perp_from_Iu("probenecid", mw = 285.36, Iu = Iu)
  inh <- inhibition_data(
    tibble::tribble(
      ~object, ~ic50,                      ~source,
      "OAT3" ,   4.4, "Gomez-Mantilla 2023 Table 17"
    ),
    precipitant = "probenecid"
  )
  tbl <- transporter_inh_risk(perp, inh)
  oat3 <- tbl[tbl$object == "OAT3", ]

  expect_equal(imaxssu(perp, molar = TRUE), Iu)
  expect_equal(oat3$r, i_over_ic50)
  expect_equal(oat3$threshold, 0.1)
  expect_gt(oat3$r, 0.1)
  expect_true(oat3$risk)
})


test_that("ibuprofen OAT3 I/IC50 recovered from Gomez-Mantilla 2023 AUCR 1.11 is a risk", {
  ic50 <- 4.4
  paper_aucr <- 1.11
  i_over_ic50 <- infer_i_over_ic50(paper_aucr, ft = 0.44)
  Iu <- i_over_ic50 * ic50
  perp <- make_perp_from_Iu("ibuprofen", mw = 206.28, Iu = Iu)
  inh <- inhibition_data(
    tibble::tribble(
      ~object, ~ic50,                      ~source,
      "OAT3" ,   4.4, "Gomez-Mantilla 2023 Table 17"
    ),
    precipitant = "ibuprofen"
  )
  tbl <- transporter_inh_risk(perp, inh)
  oat3 <- tbl[tbl$object == "OAT3", ]

  expect_equal(oat3$r, i_over_ic50)
  expect_gt(oat3$r, 0.1)
  expect_true(oat3$risk)
})


test_that("diclofenac OAT3 I/IC50 recovered from Gomez-Mantilla 2023 AUCR 1.00 is not a risk", {
  ic50 <- 3.8
  paper_aucr <- 1.00
  i_over_ic50 <- infer_i_over_ic50(paper_aucr, ft = 0.44)
  Iu <- i_over_ic50 * ic50
  perp <- make_perp_from_Iu("diclofenac", mw = 296.15, Iu = Iu)
  inh <- inhibition_data(
    tibble::tribble(
      ~object, ~ic50,                      ~source,
      "OAT3" ,   3.8, "Gomez-Mantilla 2023 Table 17"
    ),
    precipitant = "diclofenac"
  )
  tbl <- transporter_inh_risk(perp, inh)
  oat3 <- tbl[tbl$object == "OAT3", ]

  expect_equal(oat3$r, i_over_ic50)
  expect_equal(oat3$r, 0, tolerance = 1e-10)
  expect_false(oat3$risk)
})
