# Literature-based CYP perpetrator examples.
# Sources:
# Gomez-Mantilla et al., Clin Pharmacokinet 2023 (doi:10.1007/s40262-022-01204-4)
# Bergagnini-Kolev et al., AAPS J 2023 (doi:10.1208/s12248-023-00828-z)
# Almond et al., Drug Metab Dispos 2016, Table 2 (doi:10.1124/dmd.115.066845)
# Desch et al., Clin Transl Sci 2023 (doi:10.1111/cts.13470)


liver_only_midazolam <- tibble::tribble(
  ~object,  ~substrate, ~fgut,  ~fm, ~fmcyp,
  "CYP3A4", "midazolam",     1, 0.96,      1
)


test_that("ivosidenib CYP3A4 MSM without gut matches Gomez-Mantilla 2023 AUCR 0.46", {
  # Cavg,ss 9.54 uM, fu 0.057, Emax 21.25, EC50 8.6 uM.
  # Paper MSM midazolam AUCR is 0.46 when gut contribution is omitted.
  perp <- perpetrator(
    name       = "ivosidenib",
    oral       = FALSE,
    mw         = 583,
    dose       = 500,
    imaxss     = 9.54 * 583,
    fu         = 0.057,
    fumic      = 1,
    rb         = 0.56,
    fa         = 0.69,
    fg         = 1,
    ka         = 1.38 / 60,
    solubility = Inf
  )
  inh <- inhibitor(
    tibble::tribble(
      ~object, ~ki, ~source,
      "CYP3A4",  NA, "none"
    ),
    precipitant = "ivosidenib"
  )
  ind <- inducer(
    tibble::tribble(
      ~object, ~emax, ~ec50, ~max_c, ~source,
      "CYP3A4", 21.25,   8.6,  24.59, "Gomez-Mantilla 2023"
    ),
    precipitant = "ivosidenib"
  )

  tbl <- mech_stat_cyp_risk(
    perp,
    inh,
    ind,
    cyp_tdi = NULL,
    include_induction = TRUE,
    substr = liver_only_midazolam
  )@table
  Iu <- imaxssu(perp, molar = TRUE)
  Ch <- 1 + 21.25 * Iu / (Iu + 8.6)
  expected_aucr <- 1 / (Ch * 0.96 + (1 - 0.96))

  expect_equal(Iu, 9.54 * 0.057)
  expect_equal(tbl$Ch, Ch)
  expect_equal(tbl$aucr, expected_aucr)
  expect_equal(tbl$aucr, 0.46, tolerance = 0.02)
  expect_lt(tbl$aucr, 0.8)
  expect_true(tbl$risk)
})


test_that("ivosidenib is a kinetic and static CYP3A4 induction risk", {
  perp <- perpetrator(
    name       = "ivosidenib",
    oral       = FALSE,
    mw         = 583,
    dose       = 500,
    imaxss     = 9.54 * 583,
    fu         = 0.057,
    fumic      = 1,
    rb         = 0.56,
    fa         = 0.69,
    fg         = 1,
    ka         = 1.38 / 60,
    solubility = Inf
  )
  ind <- inducer(
    tibble::tribble(
      ~object, ~emax, ~ec50, ~max_c, ~source,
      "CYP3A4", 21.25,   8.6,  24.59, "Gomez-Mantilla 2023"
    ),
    precipitant = "ivosidenib"
  )

  static <- static_cyp_induction_risk(perp, ind)@table
  kinetic <- kinetic_cyp_induction_risk(perp, ind)@table
  Iu <- imaxssu(perp, molar = TRUE)
  expected_r3 <- round(
    1 / (1 + 21.25 * 10 * Iu / (8.6 + 10 * Iu)),
    3
  )

  expect_true(static$risk)
  expect_equal(kinetic$r, expected_r3)
  expect_lt(kinetic$r, 0.8)
  expect_true(kinetic$risk)
})


test_that("voxelotor CYP3A4 MSM without gut is close to Gomez-Mantilla scenario 2 AUCR 2.27", {
  perp <- perpetrator(
    name       = "voxelotor",
    oral       = FALSE,
    mw         = 337,
    dose       = 1500,
    imaxss     = 44.25 * 337,
    fu         = 0.002,
    fumic      = 1,
    rb         = 1,
    fa         = 1,
    fg         = 1,
    ka         = 0.1,
    solubility = Inf
  )
  inh <- inhibitor(
    tibble::tribble(
      ~object,  ~ki, ~source,
      "CYP3A4", 0.06, "scenario 2"
    ),
    precipitant = "voxelotor"
  )
  ind <- inducer(
    tibble::tribble(
      ~object, ~emax, ~ec50, ~max_c, ~source,
      "CYP3A4",    NA,    NA,     NA, "none"
    ),
    precipitant = "voxelotor"
  )

  tbl <- mech_stat_cyp_risk(
    perp,
    inh,
    ind,
    cyp_tdi = NULL,
    include_induction = FALSE,
    substr = liver_only_midazolam
  )@table
  Iu <- imaxssu(perp, molar = TRUE)
  Ah <- 1 / (1 + Iu / 0.06)
  expected_aucr <- 1 / (Ah * 0.96 + (1 - 0.96))

  expect_equal(Iu, 44.25 * 0.002)
  expect_equal(tbl$Ah, Ah)
  expect_equal(tbl$aucr, expected_aucr)
  expect_equal(tbl$aucr, 2.27, tolerance = 0.1)
  expect_gt(tbl$aucr, 1.25)
  expect_true(tbl$risk)
})


test_that("PUR1900 itraconazole basic R1 with OH-itraconazole matches Bergagnini-Kolev 2023", {
  perp <- perpetrator(
    name       = "PUR1900_itraconazole",
    oral       = TRUE,
    mw         = 705.6,
    dose       = 35,
    imaxss     = 0.0215 * 705.6,
    fu         = 0.016,
    fumic      = 1,
    rb         = 0.58,
    fa         = 0.159,
    fg         = 0.58,
    ka         = 0.1,
    solubility = Inf
  )
  inh <- inhibitor(
    tibble::tribble(
      ~object,    ~ki, ~source,
      "CYP3A4", 0.0013, "Bergagnini-Kolev 2023"
    ),
    precipitant = "PUR1900_itraconazole"
  )
  tbl <- basic_cyp_inhibition_risk(perp, inh)@table
  parent_r1 <- 1 + tbl$r
  oh_term <- (0.0120 * 0.016) / 0.0023
  combined_r1 <- parent_r1 + oh_term

  expect_equal(imaxssu(perp, molar = TRUE), 0.0215 * 0.016)
  expect_equal(tbl$r, round((0.0215 * 0.016) / 0.0013, 4))
  expect_true(tbl$risk_hep)
  expect_equal(combined_r1, 1.35, tolerance = 0.01)
})


test_that("rifampicin mRNA parameters from Almond 2016 are a strong kinetic CYP3A4 inducer", {
  # Indmax is fold-change; Fahmi / ddir kinetic Emax is fold-increase.
  # Unbound Cmax ~ 1.7 uM is a typical 600 mg QD value (Cmax 7 ug/ml, fu 0.2).
  perp <- perpetrator(
    name       = "rifampicin",
    oral       = TRUE,
    mw         = 822.94,
    dose       = 600,
    imaxss     = 7000,
    fu         = 0.2,
    fumic      = 1,
    rb         = 0.9,
    fa         = 1,
    fg         = 1,
    ka         = 0.01,
    solubility = Inf
  )
  ind <- inducer(
    tibble::tribble(
      ~object, ~emax, ~ec50, ~max_c, ~source,
      "CYP3A4",  28.9,  0.71,     10, "Almond 2016 Table 2"
    ),
    precipitant = "rifampicin"
  )
  static <- static_cyp_induction_risk(perp, ind)@table
  kinetic <- kinetic_cyp_induction_risk(perp, ind)@table
  Iu <- imaxssu(perp, molar = TRUE)
  expected_r3 <- round(
    1 / (1 + 28.9 * 10 * Iu / (0.71 + 10 * Iu)),
    3
  )

  expect_true(static$risk)
  expect_equal(kinetic$r, expected_r3)
  expect_lt(kinetic$r, 0.1)
  expect_true(kinetic$risk)
})


test_that("BI 425809 CYP3A4 mRNA Emax from Desch 2023 is a static induction risk", {
  ind <- inducer(
    tibble::tribble(
      ~object, ~emax, ~ec50, ~max_c, ~source,
      "CYP3A4",  18.7,  1.96,     50, "Desch 2023"
    ),
    precipitant = "BI 425809"
  )
  perp <- perpetrator(
    name       = "BI 425809",
    oral       = TRUE,
    mw         = 500,
    dose       = 25,
    imaxss     = 1000,
    fu         = 0.1,
    fumic      = 1,
    rb         = 1,
    fa         = 1,
    fg         = 1,
    ka         = 0.01,
    solubility = Inf
  )
  static <- static_cyp_induction_risk(perp, ind)@table

  expect_true(static$risk)
  expect_equal(static$emax, 18.7)
})
