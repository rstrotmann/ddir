make_risk_perp <- function(...) {
  defaults <- list(
    name       = "testdrug",
    oral       = TRUE,
    mw         = 100,
    dose       = 250,
    imaxss     = 1000,
    fu         = 0.1,
    fumic      = 1,
    rb         = 1,
    fa         = 1,
    fg         = 1,
    ka         = 0.1,
    solubility = Inf,
    source     = character(0)
  )
  do.call(perpetrator, utils::modifyList(defaults, list(...)))
}


test_that("basic_cyp_inhibition_risk returns a ddi_risk object with expected structure", {
  perp <- make_risk_perp()
  inh <- inhibitor(
    tibble::tribble(
      ~object, ~ki, ~source,
      "CYP2C9",  40, "study",
      "CYP3A4", 500, "study"
    ),
    precipitant = "testdrug"
  )

  res <- basic_cyp_inhibition_risk(perp, inh)

  expect_s4_class(res, "ddi_risk")
  expect_identical(res@precipitant, perp)
  expect_identical(res@title, "Direct CYP inhibition risk for testdrug")
  expect_named(
    res@table,
    c("object", "ki", "kiu", "source", "r", "risk_hep", "r_gut", "risk_intest")
  )
  expect_equal(res@table$object, c("CYP2C9", "CYP3A4"))
})


test_that("basic_cyp_inhibition_risk computes hepatic and intestinal R from first principles", {
  perp <- make_risk_perp()
  inh <- inhibitor(
    tibble::tribble(
      ~object ,  ~ki, ~source,
      "CYP2C9",    40, "study",
      "CYP2D6",   100, "study",
      "CYP3A4",  2000, "study"
    ),
    precipitant = "testdrug"
  )
  Iu <- imaxssu(perp, molar = TRUE)
  Igut <- igut(perp, molar = TRUE)

  res <- basic_cyp_inhibition_risk(perp, inh)
  tbl <- res@table

  expect_equal(Iu, 1)
  expect_equal(Igut, 10000)
  expect_equal(tbl$kiu, tbl$ki * perp@fumic)
  expect_equal(tbl$r, round(Iu / tbl$kiu, 4))
  expect_equal(tbl$risk_hep, (Iu / tbl$kiu) > 0.02)
  expect_equal(
    tbl$r_gut,
    c(NA_real_, NA_real_, round(Igut / 2000, 4))
  )
  expect_equal(tbl$risk_intest, c(NA, NA, (Igut / 2000) > 10))
  expect_equal(tbl$risk_hep, c(TRUE, FALSE, FALSE))
  expect_false(tbl$risk_intest[3])
})


test_that("basic_cyp_inhibition_risk flags intestinal CYP3A4 risk above the threshold of 10", {
  perp <- make_risk_perp()
  inh <- inhibitor(
    tibble::tribble(
      ~object, ~ki, ~source,
      "CYP3A4", 2000, "study"
    ),
    precipitant = "testdrug"
  )
  no_risk <- basic_cyp_inhibition_risk(perp, inh)@table
  yes_risk <- basic_cyp_inhibition_risk(
    perp,
    inhibitor(
      tibble::tribble(
        ~object, ~ki, ~source,
        "CYP3A4", 500, "study"
      ),
      precipitant = "testdrug"
    )
  )@table

  expect_false(no_risk$risk_intest)
  expect_equal(no_risk$r_gut, 5)
  expect_true(yes_risk$risk_intest)
  expect_equal(yes_risk$r_gut, 20)
})


test_that("basic_cyp_inhibition_risk scales Ki,u with fumic", {
  perp <- make_risk_perp(fumic = 0.5)
  inh <- inhibitor(
    tibble::tribble(
      ~object, ~ki, ~source,
      "CYP2C9",  40, "study"
    ),
    precipitant = "testdrug"
  )
  tbl <- basic_cyp_inhibition_risk(perp, inh)@table

  expect_equal(tbl$kiu, 20)
  expect_equal(tbl$r, round(1 / 20, 4))
})


test_that("basic_cyp_inhibition_risk sets intestinal R to zero for IV perpetrators", {
  perp <- make_risk_perp(oral = FALSE)
  inh <- inhibitor(
    tibble::tribble(
      ~object, ~ki, ~source,
      "CYP3A4",  10, "study"
    ),
    precipitant = "testdrug"
  )
  tbl <- basic_cyp_inhibition_risk(perp, inh)@table

  expect_equal(igut(perp, molar = TRUE), 0)
  expect_equal(tbl$r_gut, 0)
  expect_false(tbl$risk_intest)
})


test_that("basic_cyp_inhibition_risk keeps NA Ki rows and does not assign intestinal R to other CYPs", {
  perp <- make_risk_perp()
  inh <- inhibitor(
    tibble::tribble(
      ~object ,   ~ki, ~source,
      "CYP1A2",    NA,      "",
      "CYP3A4", 12.50, "study"
    ),
    precipitant = "testdrug"
  )
  tbl <- basic_cyp_inhibition_risk(perp, inh)@table

  expect_true(is.na(tbl$r[tbl$object == "CYP1A2"]))
  expect_true(is.na(tbl$risk_hep[tbl$object == "CYP1A2"]))
  expect_true(is.na(tbl$r_gut[tbl$object == "CYP1A2"]))
  expect_false(is.na(tbl$r_gut[tbl$object == "CYP3A4"]))
})


test_that("basic_cyp_inhibition_risk rejects wrong input classes", {
  perp <- make_risk_perp()
  inh <- inhibitor(
    tibble::tribble(
      ~object, ~ki, ~source,
      "CYP3A4",  10, "study"
    ),
    precipitant = "testdrug"
  )

  expect_error(
    basic_cyp_inhibition_risk("not-a-perp", inh),
    "perp must be a perpetrator object"
  )
  expect_error(
    basic_cyp_inhibition_risk(perp, "not-an-inhibitor"),
    "cyp_inh must be an inhibitor object"
  )
})


test_that("basic_cyp_inhibition_risk stops when no allowed CYP enzymes are present", {
  perp <- make_risk_perp()
  inh <- inhibitor(
    tibble::tribble(
      ~object, ~ki, ~source,
      "CYP3A5",  10, "study"
    ),
    precipitant = "testdrug"
  )

  expect_error(
    basic_cyp_inhibition_risk(perp, inh),
    "No inhibition data for known CYP enzymes found"
  )
})


test_that("basic_cyp_inhibition_risk warns and drops non-CYP objects", {
  perp <- make_risk_perp()
  inh <- inhibitor(
    tibble::tribble(
      ~object , ~ki, ~source,
      "CYP3A4",  10, "study",
      "UGT1A1",  20, "study"
    ),
    precipitant = "testdrug"
  )

  expect_warning(
    res <- basic_cyp_inhibition_risk(perp, inh),
    "Non-CYP data were excluded \\(UGT1A1\\)"
  )
  expect_equal(res@table$object, "CYP3A4")
})


test_that("basic_ugt_inhibition_risk treats IC50 as 2 * Ki and uses Cmax,ss,u", {
  perp <- make_risk_perp()
  inh <- inhibitor(
    as.data.frame(tibble::tribble(
      ~object, ~ic50, ~source,
      "UGT1A1",  100, "study",
      "UGT1A9",   40, "study"
    )),
    precipitant = "testdrug"
  )
  tbl <- basic_ugt_inhibition_risk(perp, inh)@table
  Iu <- imaxssu(perp, molar = TRUE)

  expect_s4_class(basic_ugt_inhibition_risk(perp, inh), "ddi_risk")
  expect_named(tbl, c("object", "ki", "kiu", "source", "r", "risk"))
  expect_equal(tbl$ki, c(50, 20))
  expect_equal(tbl$kiu, tbl$ki * perp@fumic)
  expect_equal(tbl$r, round(Iu / tbl$kiu, 4))
  expect_equal(tbl$risk, c(FALSE, TRUE))
  expect_identical(
    basic_ugt_inhibition_risk(perp, inh)@title,
    "UGT inhibition risk for testdrug"
  )
})


test_that("basic_ugt_inhibition_risk rejects wrong input classes and empty UGT data", {
  perp <- make_risk_perp()
  ugt <- inhibitor(
    as.data.frame(tibble::tribble(
      ~object, ~ic50, ~source,
      "UGT1A1",   15, "study"
    )),
    precipitant = "testdrug"
  )
  cyp_only <- inhibitor(
    as.data.frame(tibble::tribble(
      ~object, ~ic50, ~source,
      "CYP3A4",   10, "study"
    )),
    precipitant = "testdrug"
  )

  expect_error(
    basic_ugt_inhibition_risk("not-a-perp", ugt),
    "perp must be a perpetrator object"
  )
  expect_error(
    basic_ugt_inhibition_risk(perp, "not-an-inhibitor"),
    "ugt_inh must be an inhibitor object"
  )
  expect_error(
    basic_ugt_inhibition_risk(perp, cyp_only),
    "No inhibition data for known UGT enzymes found"
  )
})


test_that("basic_ugt_inhibition_risk drops CYP rows without warning when the data use object not target", {
  perp <- make_risk_perp()
  inh <- inhibitor(
    as.data.frame(tibble::tribble(
      ~object, ~ic50, ~source,
      "UGT1A1",    40, "study",
      "CYP3A4",    10, "study"
    )),
    precipitant = "testdrug"
  )

  expect_no_warning(res <- basic_ugt_inhibition_risk(perp, inh))
  expect_equal(res@table$object, "UGT1A1")
})


test_that("transporter_inh_risk expands Pgp and BCRP into intestinal and systemic rows", {
  perp <- make_risk_perp()
  inh <- inhibitor(
    tibble::tribble(
      ~object, ~ic50, ~source,
      "Pgp"  ,   100, "study",
      "BCRP" ,  2000, "study",
      "OAT1" ,    20, "study"
    ),
    precipitant = "testdrug"
  )
  tbl <- transporter_inh_risk(perp, inh)@table

  expect_s4_class(transporter_inh_risk(perp, inh), "ddi_risk")
  expect_false(any(tbl$object %in% c("Pgp", "BCRP")))
  expect_setequal(
    tbl$object,
    c("Pgp_int", "Pgp_sys", "BCRP_int", "BCRP_sys", "OAT1")
  )
  expect_named(
    tbl,
    c("object", "ic50", "source", "i", "r", "threshold", "risk")
  )
  expect_identical(
    transporter_inh_risk(perp, inh)@title,
    "Transporter inhibition risk for testdrug"
  )
})


test_that("transporter_inh_risk uses the ICH concentration metric and threshold for each row", {
  perp <- make_risk_perp()
  inh <- inhibitor(
    tibble::tribble(
      ~object , ~ic50, ~source,
      "Pgp"   ,   100, "study",
      "OATP1B1",   200, "study",
      "MATE1" ,    25, "study"
    ),
    precipitant = "testdrug"
  )
  tbl <- transporter_inh_risk(perp, inh)@table
  Igut <- igut(perp, molar = TRUE)
  Iu <- imaxssu(perp, molar = TRUE)
  Ih <- imaxinletu(perp, molar = TRUE)

  pgp_int <- tbl[tbl$object == "Pgp_int", ]
  pgp_sys <- tbl[tbl$object == "Pgp_sys", ]
  oatp <- tbl[tbl$object == "OATP1B1", ]
  mate <- tbl[tbl$object == "MATE1", ]

  expect_equal(pgp_int$i, "igut")
  expect_equal(pgp_int$threshold, 10)
  expect_equal(pgp_int$r, Igut / 100)
  expect_true(pgp_int$risk)

  expect_equal(pgp_sys$i, "imaxssu")
  expect_equal(pgp_sys$threshold, 0.02)
  expect_equal(pgp_sys$r, Iu / 100)
  expect_false(pgp_sys$risk)

  expect_equal(oatp$i, "imaxinletu")
  expect_equal(oatp$threshold, 0.1)
  expect_equal(oatp$r, Ih / 200)

  expect_equal(mate$i, "imaxssu")
  expect_equal(mate$threshold, 0.02)
  expect_equal(mate$r, Iu / 25)
  expect_true(mate$risk)
})


test_that("transporter_inh_risk leaves r as NA when IC50 is missing", {
  perp <- make_risk_perp()
  inh <- inhibitor(
    tibble::tribble(
      ~object, ~ic50, ~source,
      "OAT1" ,    NA, "study"
    ),
    precipitant = "testdrug"
  )
  tbl <- transporter_inh_risk(perp, inh)@table

  expect_true(is.na(tbl$r))
  expect_true(is.na(tbl$risk))
})


test_that("transporter_inh_risk passes qh through to hepatic inlet concentration", {
  perp <- make_risk_perp()
  inh <- inhibitor(
    tibble::tribble(
      ~object , ~ic50, ~source,
      "OATP1B1",   10, "study"
    ),
    precipitant = "testdrug"
  )
  default <- transporter_inh_risk(perp, inh, qh = 1.616)@table
  custom <- transporter_inh_risk(perp, inh, qh = 0.808)@table

  expect_equal(default$r, imaxinletu(perp, qh = 1.616, molar = TRUE) / 10)
  expect_equal(custom$r, imaxinletu(perp, qh = 0.808, molar = TRUE) / 10)
  expect_gt(custom$r, default$r)
})


test_that("transporter_inh_risk uses a custom reference table when supplied", {
  perp <- make_risk_perp()
  inh <- inhibitor(
    tibble::tribble(
      ~object, ~ic50, ~source,
      "OCT2" ,    10, "study"
    ),
    precipitant = "testdrug"
  )
  ref <- tibble::tribble(
    ~object, ~threshold,        ~i,
    "OCT2" ,          1, "imaxssu"
  )
  tbl <- transporter_inh_risk(perp, inh, transporter_ref = ref)@table

  expect_equal(tbl$threshold, 1)
  expect_equal(tbl$r, 1 / 10)
  expect_false(tbl$risk)
})


test_that("transporter_inh_risk rejects wrong classes and empty transporter data", {
  perp <- make_risk_perp()
  inh <- inhibitor(
    tibble::tribble(
      ~object, ~ic50, ~source,
      "OAT1" ,    10, "study"
    ),
    precipitant = "testdrug"
  )
  cyp_only <- inhibitor(
    tibble::tribble(
      ~object, ~ic50, ~source,
      "CYP3A4",   10, "study"
    ),
    precipitant = "testdrug"
  )

  expect_error(
    transporter_inh_risk("not-a-perp", inh),
    "perp must be a perpetrator object"
  )
  expect_error(
    transporter_inh_risk(perp, "not-an-inhibitor"),
    "transporter_inh must be an inhibitor object"
  )
  expect_error(
    transporter_inh_risk(perp, cyp_only),
    "No inhibition data for known transporters found"
  )
})


test_that("transporter_inh_risk warns and drops non-transporter objects", {
  perp <- make_risk_perp()
  inh <- inhibitor(
    tibble::tribble(
      ~object , ~ic50, ~source,
      "OAT1"  ,    10, "study",
      "CYP3A4",    20, "study"
    ),
    precipitant = "testdrug"
  )

  expect_warning(
    res <- transporter_inh_risk(perp, inh),
    "Non-transporter data were excluded \\(CYP3A4\\)"
  )
  expect_equal(res@table$object, "OAT1")
})


test_that("basic_cyp_tdi_risk computes R from kobs and hepatic kdeg", {
  perp <- make_risk_perp()
  tdi <- inhibitor(
    tibble::tribble(
      ~object,  ~ki, ~kinact, ~source,
      "CYP3A4",   1,     0.1, "study"
    ),
    precipitant = "testdrug"
  )
  Iu <- imaxssu(perp, molar = TRUE)
  kobs <- 0.1 * 5 * Iu / (1 * perp@fumic + 5 * Iu)
  kdeg <- cyp_turnover$kdeg_hepatic[cyp_turnover$object == "CYP3A4"]
  expected_r <- (kobs + kdeg) / kdeg
  tbl <- basic_cyp_tdi_risk(perp, tdi)@table

  expect_s4_class(basic_cyp_tdi_risk(perp, tdi), "ddi_risk")
  expect_named(
    tbl,
    c("object", "ki", "fu", "kinact", "kdeg", "source", "r", "risk")
  )
  expect_equal(tbl$fu, perp@fu)
  expect_equal(tbl$kdeg, kdeg)
  expect_equal(tbl$r, expected_r)
  expect_true(tbl$risk)
  expect_identical(
    basic_cyp_tdi_risk(perp, tdi)@title,
    "Time-dependent CYP inhibition risk for testdrug"
  )
})


test_that("basic_cyp_tdi_risk is negative when kobs is small relative to kdeg", {
  perp <- make_risk_perp()
  tdi <- inhibitor(
    tibble::tribble(
      ~object,  ~ki, ~kinact, ~source,
      "CYP3A4",   1,   0.001, "study"
    ),
    precipitant = "testdrug"
  )
  tbl <- basic_cyp_tdi_risk(perp, tdi)@table

  expect_lt(tbl$r, 1.25)
  expect_false(tbl$risk)
})


test_that("basic_cyp_tdi_risk uses fumic in the kobs denominator", {
  perp_full <- make_risk_perp(fumic = 1)
  perp_half <- make_risk_perp(fumic = 0.5)
  tdi <- inhibitor(
    tibble::tribble(
      ~object,  ~ki, ~kinact, ~source,
      "CYP3A4",   1,     0.1, "study"
    ),
    precipitant = "testdrug"
  )

  r_full <- basic_cyp_tdi_risk(perp_full, tdi)@table$r
  r_half <- basic_cyp_tdi_risk(perp_half, tdi)@table$r

  expect_gt(r_half, r_full)
})


test_that("basic_cyp_tdi_risk can use a custom turnover table", {
  perp <- make_risk_perp()
  tdi <- inhibitor(
    tibble::tribble(
      ~object,  ~ki, ~kinact, ~source,
      "CYP3A4",   1,     0.1, "study"
    ),
    precipitant = "testdrug"
  )
  kdeg <- tibble::tribble(
    ~object, ~kdeg_hepatic, ~kdeg_intestinal,
    "CYP3A4",          0.5,             0.03
  )
  Iu <- imaxssu(perp, molar = TRUE)
  kobs <- 0.1 * 5 * Iu / (1 * 1 + 5 * Iu)
  tbl <- basic_cyp_tdi_risk(perp, tdi, cyp_kdeg = kdeg)@table

  expect_equal(tbl$kdeg, 0.5)
  expect_equal(tbl$r, (kobs + 0.5) / 0.5)
})


test_that("basic_cyp_tdi_risk rejects wrong classes, missing columns, and empty CYP TDI", {
  perp <- make_risk_perp()
  tdi <- inhibitor(
    tibble::tribble(
      ~object,  ~ki, ~kinact, ~source,
      "CYP3A4",   1,     0.1, "study"
    ),
    precipitant = "testdrug"
  )
  no_kinact <- inhibitor(
    tibble::tribble(
      ~object, ~ki, ~source,
      "CYP3A4",  1, "study"
    ),
    precipitant = "testdrug"
  )
  only_3a5 <- inhibitor(
    tibble::tribble(
      ~object,  ~ki, ~kinact, ~source,
      "CYP3A5",   1,     0.1, "study"
    ),
    precipitant = "testdrug"
  )

  expect_error(
    basic_cyp_tdi_risk("not-a-perp", tdi),
    "perp must be a perpetrotor object"
  )
  expect_error(
    basic_cyp_tdi_risk(perp, "not-an-inhibitor"),
    "cyp_tdi must be an inhibitor object"
  )
  expect_error(
    basic_cyp_tdi_risk(perp, no_kinact),
    "Missing columns in cyp_tdi: kinact"
  )
  expect_error(
    basic_cyp_tdi_risk(perp, only_3a5),
    "No TDI data for known CYP enzymes found"
  )
})


test_that("basic_cyp_tdi_risk warns about CYP3A5 but still includes it in the fitted table", {
  perp <- make_risk_perp()
  tdi <- inhibitor(
    tibble::tribble(
      ~object,  ~ki, ~kinact, ~source,
      "CYP3A4",   1,     0.1, "study",
      "CYP3A5",   1,     0.1, "study"
    ),
    precipitant = "testdrug"
  )

  expect_warning(
    res <- basic_cyp_tdi_risk(perp, tdi),
    "Non-CYP data were excluded \\(CYP3A5\\)"
  )
  expect_equal(res@table$object, c("CYP3A4", "CYP3A5"))
})


test_that("static_cyp_induction_risk flags Emax of 2 or more and notes insufficient max_c", {
  perp <- make_risk_perp()
  ind <- inducer(
    tibble::tribble(
      ~object , ~emax, ~ec50, ~max_c, ~source,
      "CYP1A2",   1.5,   1.0,    100, "study",
      "CYP3A4",   2.0,   1.0,     10, "study"
    ),
    precipitant = "testdrug"
  )
  tbl <- static_cyp_induction_risk(perp, ind)@table

  expect_s4_class(static_cyp_induction_risk(perp, ind), "ddi_risk")
  expect_false("ec50" %in% names(tbl))
  expect_false("maxc_imaxssu" %in% names(tbl))
  expect_equal(tbl$risk, c(FALSE, TRUE))
  expect_equal(
    tbl$note,
    c("", "Not tested up to 50-fold Cmax,u")
  )
  expect_identical(
    static_cyp_induction_risk(perp, ind)@title,
    "Static CYP induction risk for testdrug"
  )
})


test_that("static_cyp_induction_risk leaves the note empty when max_c covers 50-fold Cmax,u", {
  perp <- make_risk_perp()
  ind <- inducer(
    tibble::tribble(
      ~object, ~emax, ~ec50, ~max_c, ~source,
      "CYP3A4",   3.0,   1.0,     50, "study"
    ),
    precipitant = "testdrug"
  )
  tbl <- static_cyp_induction_risk(perp, ind)@table

  expect_equal(tbl$note, "")
  expect_true(tbl$risk)
})


test_that("static_cyp_induction_risk rejects wrong classes and empty allowed CYP data", {
  perp <- make_risk_perp()
  ind <- inducer(
    tibble::tribble(
      ~object, ~emax, ~ec50, ~max_c, ~source,
      "CYP3A4",    2,     1,    100, "study"
    ),
    precipitant = "testdrug"
  )
  only_3a5 <- inducer(
    tibble::tribble(
      ~object, ~emax, ~ec50, ~max_c, ~source,
      "CYP3A5",    3,     1,    100, "study"
    ),
    precipitant = "testdrug"
  )

  expect_error(
    static_cyp_induction_risk("not-a-perp", ind),
    "perp must be a perpetrotor object"
  )
  expect_error(
    static_cyp_induction_risk(perp, "not-an-inducer"),
    "cyp_ind must be an inducer object"
  )

  expect_error(
    static_cyp_induction_risk(perp, only_3a5),
    "No CYP induction data for known CYP enzymes found"
  )
})


test_that("static_cyp_induction_risk warns about CYP3A5 but still includes it in the table", {
  perp <- make_risk_perp()
  ind <- inducer(
    tibble::tribble(
      ~object, ~emax, ~ec50, ~max_c, ~source,
      "CYP3A4",    2,     1,    100, "study",
      "CYP3A5",    3,     1,    100, "study"
    ),
    precipitant = "testdrug"
  )

  expect_warning(
    res <- static_cyp_induction_risk(perp, ind),
    "Non-CYP data were excluded \\(CYP3A5\\)"
  )
  expect_equal(res@table$object, c("CYP3A4", "CYP3A5"))
})


test_that("kinetic_cyp_induction_risk uses the R3-style formula and 0.8 cutoff", {
  perp <- make_risk_perp()
  ind <- inducer(
    tibble::tribble(
      ~object, ~emax, ~ec50, ~max_c, ~source,
      "CYP3A4",     5,     1,    100, "study"
    ),
    precipitant = "testdrug"
  )
  Iu <- imaxssu(perp, molar = TRUE)
  expected_r <- round(1 / (1 + 1 * 5 * 10 * Iu / (1 + 10 * Iu)), 3)
  tbl <- kinetic_cyp_induction_risk(perp, ind)@table

  expect_s4_class(kinetic_cyp_induction_risk(perp, ind), "ddi_risk")
  expect_equal(tbl$r, expected_r)
  expect_true(tbl$risk)
  expect_identical(
    kinetic_cyp_induction_risk(perp, ind)@title,
    "Kinetic CYP induction risk for testdrug"
  )
})


test_that("kinetic_cyp_induction_risk is negative for weak induction and when d is 0", {
  perp <- make_risk_perp()
  weak <- inducer(
    tibble::tribble(
      ~object, ~emax, ~ec50, ~max_c, ~source,
      "CYP3A4",   0.1,   100,    100, "study"
    ),
    precipitant = "testdrug"
  )
  strong <- inducer(
    tibble::tribble(
      ~object, ~emax, ~ec50, ~max_c, ~source,
      "CYP3A4",     5,     1,    100, "study"
    ),
    precipitant = "testdrug"
  )

  expect_false(kinetic_cyp_induction_risk(perp, weak)@table$risk)
  expect_gt(kinetic_cyp_induction_risk(perp, weak)@table$r, 0.8)
  expect_equal(kinetic_cyp_induction_risk(perp, strong, d = 0)@table$r, 1)
  expect_false(kinetic_cyp_induction_risk(perp, strong, d = 0)@table$risk)
})


test_that("kinetic_cyp_induction_risk rejects wrong classes and empty allowed CYP data", {
  perp <- make_risk_perp()
  ind <- inducer(
    tibble::tribble(
      ~object, ~emax, ~ec50, ~max_c, ~source,
      "CYP3A4",    2,     1,    100, "study"
    ),
    precipitant = "testdrug"
  )
  only_3a5 <- inducer(
    tibble::tribble(
      ~object, ~emax, ~ec50, ~max_c, ~source,
      "CYP3A5",    3,     1,    100, "study"
    ),
    precipitant = "testdrug"
  )

  expect_error(
    kinetic_cyp_induction_risk("not-a-perp", ind),
    "perp must be a perpetrator object"
  )
  expect_error(
    kinetic_cyp_induction_risk(perp, "not-an-inducer"),
    "cyp_ind must be an inducer object"
  )
  expect_error(
    kinetic_cyp_induction_risk(perp, only_3a5),
    "No CYP induction data for known CYP enzymes found"
  )
})


test_that("kinetic_cyp_induction_risk warns about CYP3A5 but still includes it in the table", {
  perp <- make_risk_perp()
  ind <- inducer(
    tibble::tribble(
      ~object, ~emax, ~ec50, ~max_c, ~source,
      "CYP3A4",    2,     1,    100, "study",
      "CYP3A5",    3,     1,    100, "study"
    ),
    precipitant = "testdrug"
  )

  expect_warning(
    res <- kinetic_cyp_induction_risk(perp, ind),
    "Non-CYP data were excluded \\(CYP3A5\\)"
  )
  expect_equal(res@table$object, c("CYP3A4", "CYP3A5"))
})


test_that("mech_stat_cyp_risk returns AUCR from Ag, Bg, Cg, Ah, Bh and Ch", {
  perp <- make_risk_perp()
  inh <- inhibitor(
    tibble::tribble(
      ~object, ~ki, ~source,
      "CYP3A4",  10, "study"
    ),
    precipitant = "testdrug"
  )
  ind <- inducer(
    tibble::tribble(
      ~object, ~emax, ~ec50, ~max_c, ~source,
      "CYP3A4",     5,     1,    100, "study"
    ),
    precipitant = "testdrug"
  )
  tdi <- inhibitor(
    tibble::tribble(
      ~object,  ~ki, ~kinact, ~source,
      "CYP3A4",   2,    0.05, "study"
    ),
    precipitant = "testdrug"
  )
  Ig <- imaxintest(perp)
  Ih <- imaxinletu(perp)
  kiu <- 10 * perp@fumic
  kdeg_h <- cyp_turnover$kdeg_hepatic[cyp_turnover$object == "CYP3A4"]
  kdeg_g <- cyp_turnover$kdeg_intestinal[cyp_turnover$object == "CYP3A4"]
  Ag <- 1 / (1 + Ig / kiu)
  Ah <- 1 / (1 + Ih / kiu)
  Bg <- kdeg_g / (kdeg_g + (Ig * 0.05 / (Ig + 2)))
  Bh <- kdeg_h / (kdeg_h + (Ih * 0.05 / (Ih + 2)))
  Cg <- 1 + (1 * 5 * Ig / (Ig + 1))
  Ch <- 1 + (1 * 5 * Ih / (Ih + 1))
  fgut <- 0.57
  fm <- 0.96
  fmcyp <- 1
  expected_aucr <- 1 / (Ag * Bg * Cg * (1 - fgut) + fgut) *
    1 / (Ah * Bh * Ch * fm * fmcyp + (1 - fm * fmcyp))

  tbl <- mech_stat_cyp_risk(perp, inh, ind, tdi)@table

  expect_s4_class(mech_stat_cyp_risk(perp, inh, ind, tdi), "ddi_risk")
  expect_equal(tbl$object, "CYP3A4")
  expect_equal(tbl$substrate, "midazolam")
  expect_equal(tbl$kiu, kiu)
  expect_equal(tbl$Ag, Ag)
  expect_equal(tbl$Ah, Ah)
  expect_equal(tbl$Bg, Bg)
  expect_equal(tbl$Bh, Bh)
  expect_equal(tbl$Cg, Cg)
  expect_equal(tbl$Ch, Ch)
  expect_equal(tbl$aucr, expected_aucr)
  expect_equal(tbl$risk, expected_aucr > 1.25 | expected_aucr < 0.8)
  expect_identical(
    mech_stat_cyp_risk(perp, inh, ind, tdi)@title,
    "Mechanistic-static CYP inhibition risk for testdrug"
  )
})


test_that("mech_stat_cyp_risk sets inhibition terms to 1 when Ki is NA", {
  perp <- make_risk_perp()
  inh <- inhibitor(
    tibble::tribble(
      ~object, ~ki, ~source,
      "CYP3A4",  NA, "none"
    ),
    precipitant = "testdrug"
  )
  ind <- inducer(
    tibble::tribble(
      ~object, ~emax, ~ec50, ~max_c, ~source,
      "CYP3A4",    NA,    NA,     NA, "none"
    ),
    precipitant = "testdrug"
  )
  tbl <- mech_stat_cyp_risk(perp, inh, ind, cyp_tdi = NULL)@table

  expect_equal(tbl$Ag, 1)
  expect_equal(tbl$Ah, 1)
  expect_equal(tbl$Bg, 1)
  expect_equal(tbl$Bh, 1)
  expect_equal(tbl$Cg, 1)
  expect_equal(tbl$Ch, 1)
})


test_that("mech_stat_cyp_risk can omit TDI and turn induction off", {
  perp <- make_risk_perp()
  inh <- inhibitor(
    tibble::tribble(
      ~object, ~ki, ~source,
      "CYP3A4",  10, "study"
    ),
    precipitant = "testdrug"
  )
  ind <- inducer(
    tibble::tribble(
      ~object, ~emax, ~ec50, ~max_c, ~source,
      "CYP3A4",     5,     1,    100, "study"
    ),
    precipitant = "testdrug"
  )
  no_tdi <- mech_stat_cyp_risk(perp, inh, ind, cyp_tdi = NULL)@table
  no_ind <- mech_stat_cyp_risk(
    perp, inh, ind, cyp_tdi = NULL, include_induction = FALSE
  )@table

  expect_equal(no_tdi$Bg, 1)
  expect_equal(no_tdi$Bh, 1)
  expect_gt(no_tdi$Cg, 1)
  expect_gt(no_tdi$Ch, 1)
  expect_equal(no_ind$Cg, 1)
  expect_equal(no_ind$Ch, 1)
})


test_that("mech_stat_cyp_risk scales induction with d and concentrations with qh and qent", {
  perp <- make_risk_perp()
  inh <- inhibitor(
    tibble::tribble(
      ~object, ~ki, ~source,
      "CYP3A4",  10, "study"
    ),
    precipitant = "testdrug"
  )
  ind <- inducer(
    tibble::tribble(
      ~object, ~emax, ~ec50, ~max_c, ~source,
      "CYP3A4",     5,     1,    100, "study"
    ),
    precipitant = "testdrug"
  )
  d0 <- mech_stat_cyp_risk(perp, inh, ind, d = 0, cyp_tdi = NULL)@table
  d1 <- mech_stat_cyp_risk(perp, inh, ind, d = 1, cyp_tdi = NULL)@table
  qent <- mech_stat_cyp_risk(
    perp, inh, ind, cyp_tdi = NULL, include_induction = FALSE, qent = 9 / 60
  )@table
  default_qent <- mech_stat_cyp_risk(
    perp, inh, ind, cyp_tdi = NULL, include_induction = FALSE
  )@table

  expect_equal(d0$Cg, 1)
  expect_equal(d0$Ch, 1)
  expect_gt(d1$Cg, d0$Cg)
  expect_lt(qent$Ag, default_qent$Ag)
})


test_that("mech_stat_cyp_risk uses a custom substrate table", {
  perp <- make_risk_perp()
  inh <- inhibitor(
    tibble::tribble(
      ~object, ~ki, ~source,
      "CYP3A4",  10, "study"
    ),
    precipitant = "testdrug"
  )
  ind <- inducer(
    tibble::tribble(
      ~object, ~emax, ~ec50, ~max_c, ~source,
      "CYP3A4",    NA,    NA,     NA, "none"
    ),
    precipitant = "testdrug"
  )
  substr <- tibble::tribble(
    ~object,  ~substrate, ~fgut, ~fm, ~fmcyp,
    "CYP3A4", "testsub" ,     1,   1,      1
  )
  tbl <- mech_stat_cyp_risk(
    perp, inh, ind, cyp_tdi = NULL, include_induction = FALSE, substr = substr
  )@table

  expect_equal(tbl$substrate, "testsub")
  expect_equal(tbl$fgut, 1)
  expect_equal(tbl$fm, 1)
  expect_equal(tbl$fmcyp, 1)
})


test_that("mech_stat_cyp_risk flags AUCR above 1.25 or below 0.8", {
  perp <- make_risk_perp()
  strong_inh <- inhibitor(
    tibble::tribble(
      ~object, ~ki, ~source,
      "CYP3A4", 0.01, "study"
    ),
    precipitant = "testdrug"
  )
  no_ind <- inducer(
    tibble::tribble(
      ~object, ~emax, ~ec50, ~max_c, ~source,
      "CYP3A4",    NA,    NA,     NA, "none"
    ),
    precipitant = "testdrug"
  )
  strong_ind <- inducer(
    tibble::tribble(
      ~object, ~emax, ~ec50, ~max_c, ~source,
      "CYP3A4",    50,  0.01,    100, "study"
    ),
    precipitant = "testdrug"
  )
  weak_inh <- inhibitor(
    tibble::tribble(
      ~object, ~ki, ~source,
      "CYP3A4",  NA, "none"
    ),
    precipitant = "testdrug"
  )

  inh_tbl <- mech_stat_cyp_risk(
    perp, strong_inh, no_ind, cyp_tdi = NULL, include_induction = FALSE
  )@table
  ind_tbl <- mech_stat_cyp_risk(
    perp, weak_inh, strong_ind, cyp_tdi = NULL, include_induction = TRUE
  )@table

  expect_gt(inh_tbl$aucr, 1.25)
  expect_true(inh_tbl$risk)
  expect_lt(ind_tbl$aucr, 0.8)
  expect_true(ind_tbl$risk)
})


test_that("mech_stat_cyp_risk rejects wrong input classes", {
  perp <- make_risk_perp()
  inh <- inhibitor(
    tibble::tribble(
      ~object, ~ki, ~source,
      "CYP3A4",  10, "study"
    ),
    precipitant = "testdrug"
  )
  ind <- inducer(
    tibble::tribble(
      ~object, ~emax, ~ec50, ~max_c, ~source,
      "CYP3A4",    2,     1,    100, "study"
    ),
    precipitant = "testdrug"
  )

  expect_error(
    mech_stat_cyp_risk("not-a-perp", inh, ind),
    "perp must be a perpetrator object"
  )
  expect_error(
    mech_stat_cyp_risk(perp, "not-an-inhibitor", ind),
    "cyp_inh must be an inhibitor object"
  )
  expect_error(
    mech_stat_cyp_risk(perp, inh, "not-an-inducer"),
    "cyp_ind must be an inducer object"
  )
  expect_error(
    mech_stat_cyp_risk(perp, inh, ind, cyp_tdi = "not-an-inhibitor"),
    "cyp_tdi must be an inhibitor object"
  )
  expect_error(
    mech_stat_cyp_risk(perp, inh, NULL),
    "cyp_ind must be an inducer object"
  )
})


test_that("examplinib parent CYP inhibition risk matches independent R values", {
  Iu <- imaxssu(examplinib_parent, molar = TRUE)
  Igut <- igut(examplinib_parent, molar = TRUE)
  tbl <- basic_cyp_inhibition_risk(
    examplinib_parent, examplinib_cyp_inh_parent
  )@table
  cyp3a4 <- tbl[tbl$object == "CYP3A4", ]
  cyp2c19 <- tbl[tbl$object == "CYP2C19", ]

  expect_equal(cyp3a4$r, round(Iu / 12.5, 4))
  expect_equal(cyp3a4$r_gut, round(Igut / 12.5, 4))
  expect_true(cyp3a4$risk_intest)
  expect_equal(cyp2c19$r, round(Iu / 0.25, 4))
  expect_true(cyp2c19$risk_hep)
  expect_true(is.na(cyp2c19$r_gut))
})


test_that("examplinib parent UGT, TDI, induction and transporter assessments run", {
  expect_no_error(
    ugt <- basic_ugt_inhibition_risk(
      examplinib_parent, examplinib_ugt_inh_parent
    )
  )
  expect_no_error(
    tdi <- basic_cyp_tdi_risk(
      examplinib_parent, examplinib_cyp_tdi_parent
    )
  )
  expect_no_error(
    static <- static_cyp_induction_risk(
      examplinib_parent, examplinib_cyp_ind_parent
    )
  )
  expect_no_error(
    kinetic <- kinetic_cyp_induction_risk(
      examplinib_parent, examplinib_cyp_ind_parent
    )
  )
  expect_no_error(
    transp <- transporter_inh_risk(
      examplinib_parent, examplinib_transporter_inh_parent
    )
  )
  expect_no_error(
    msm <- mech_stat_cyp_risk(
      examplinib_parent,
      examplinib_cyp_inh_parent,
      examplinib_cyp_ind_parent,
      examplinib_cyp_tdi_parent
    )
  )

  expect_true(all(c("UGT1A1", "UGT1A9") %in% ugt@table$object))
  expect_equal(tdi@table$object, "CYP3A4")
  expect_true(static@table$risk[static@table$object == "CYP3A4"])
  expect_true(kinetic@table$risk[kinetic@table$object == "CYP3A4"])
  expect_true(all(c("Pgp_int", "Pgp_sys", "OATP1B1") %in% transp@table$object))
  expect_true("CYP3A4" %in% msm@table$object)
  expect_false(anyNA(msm@table$aucr[msm@table$object == "CYP3A4"]))
})


test_that("examplinib metabolite uses IV concentration branches in risk tables", {
  expect_equal(igut(examplinib_metabolite, molar = TRUE), 0)
  cyp <- basic_cyp_inhibition_risk(
    examplinib_metabolite, examplinib_cyp_inh_metabolite
  )@table
  cyp3a4 <- cyp[cyp$object == "CYP3A4", ]
  transp <- transporter_inh_risk(
    examplinib_parent, examplinib_transporter_inh_parent
  )@table

  if (nrow(cyp3a4) > 0 && !all(is.na(cyp3a4$ki))) {
    expect_equal(cyp3a4$r_gut, 0)
  }
  expect_s4_class(
    basic_ugt_inhibition_risk(
      examplinib_metabolite, examplinib_ugt_inh_metabolite
    ),
    "ddi_risk"
  )
  expect_true("Pgp_int" %in% transp$object)
})


test_that("ddi_risk objects returned by the assessment functions can be shown and printed", {
  perp <- make_risk_perp()
  inh <- inhibitor(
    tibble::tribble(
      ~object, ~ki, ~source,
      "CYP3A4",  10, "study"
    ),
    precipitant = "testdrug"
  )
  res <- basic_cyp_inhibition_risk(perp, inh)
  show_text <- paste(capture.output(show(res)), collapse = "\n")
  print_text <- paste(capture.output(print(res)), collapse = "\n")

  expect_match(show_text, "Direct CYP inhibition risk for testdrug")
  expect_match(print_text, "Direct CYP inhibition risk for testdrug")
  expect_s3_class(print(res), "knitr_kable")
})
