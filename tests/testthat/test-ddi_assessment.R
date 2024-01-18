test_that("integrated ddi assessment for single drug runs without error", {
  expect_no_error({
    compounds <- read_perpetrators(
      test_path("fixtures", "examplinib_compounds_single.csv"))

    cyp_inhibition_data <- read_cyp_inhibitor_data(
      test_path("fixtures", "examplinib_cyp_inhibition.csv"))

    cyp_tdi_data <- read_tdi_data(
      test_path("fixtures", "examplinib_cyp_tdi.csv"))

    cyp_induction_data <- read_inducer_data(
      test_path("fixtures", "examplinib_cyp_induction.csv"))

    ugt_inhibition_data <- read_ugt_inhibitor_data(
      test_path("fixtures", "examplinib_ugt_inhibition.csv"))

    transporter_inhibition_data <- read_transporter_inhibitor_data(
      test_path("fixtures", "examplinib_transporter_inhibition.csv"))

    cat(compound_names_string(compounds))

    p <- compounds

    property_table(p)
    conc_table(p)

    basic_cyp_inhibition_risk_table(p, cyp_inhibition_data)
    basic_cyp_tdi_risk_table(p, cyp_tdi_data, na.rm=T)
    static_cyp_induction_risk_table(p, cyp_induction_data, na.rm=T)
    kinetic_cyp_induction_risk_table(p, cyp_induction_data, na.rm=T)
    mech_stat_cyp_risk_table(p, cyp_inhibition_data, cyp_induction_data,
                             cyp_tdi_data, include_induction=F)
    mech_stat_cyp_risk_table(p, cyp_inhibition_data, cyp_induction_data,
                             cyp_tdi_data, include_induction=T, na.rm=T)
    basic_ugt_inhibition_risk_table(p, ugt_inhibition_data)
    transporter_inhibition_risk_table(p, transporter_inhibition_data)
  })
})


test_that("integrated ddi assessment for drug with metabolite runs without error", {
  expect_no_error({
    compounds <- read_perpetrators(
      test_path("fixtures", "examplinib_compounds.csv"))

    cyp_inhibition_data <- read_cyp_inhibitor_data(
      test_path("fixtures", "examplinib_cyp_inhibition.csv"))

    cyp_tdi_data <- read_tdi_data(
      test_path("fixtures", "examplinib_cyp_tdi.csv"))

    cyp_induction_data <- read_inducer_data(
      test_path("fixtures", "examplinib_cyp_induction.csv"))

    ugt_inhibition_data <- read_ugt_inhibitor_data(
      test_path("fixtures", "examplinib_ugt_inhibition.csv"))

    transporter_inhibition_data <- read_transporter_inhibitor_data(
      test_path("fixtures", "examplinib_transporter_inhibition.csv"))

    cat(compound_names_string(compounds))

    p <- compounds[[1]]
    m <- compounds[[2]]

    property_table(p)
    property_table(m)
    conc_table(p)
    conc_table(m)

    basic_cyp_inhibition_risk_table(p, cyp_inhibition_data)
    basic_cyp_inhibition_risk_table(m, cyp_inhibition_data)

    basic_cyp_tdi_risk_table(p, cyp_tdi_data, na.rm=T)
    basic_cyp_tdi_risk_table(m, cyp_tdi_data, na.rm=T)

    static_cyp_induction_risk_table(p, cyp_induction_data, na.rm=T)
    static_cyp_induction_risk_table(m, cyp_induction_data, na.rm=T)

    kinetic_cyp_induction_risk_table(p, cyp_induction_data, na.rm=T)
    kinetic_cyp_induction_risk_table(m, cyp_induction_data, na.rm=T)

    mech_stat_cyp_risk_table(p, cyp_inhibition_data, cyp_induction_data,
                             cyp_tdi_data, include_induction=F)
    mech_stat_cyp_risk_table(m, cyp_inhibition_data, cyp_induction_data,
                             cyp_tdi_data, include_induction=F)

    mech_stat_cyp_risk_table(p, cyp_inhibition_data, cyp_induction_data,
                             cyp_tdi_data, include_induction=T, na.rm=T)
    mech_stat_cyp_risk_table(m, cyp_inhibition_data, cyp_induction_data,
                             cyp_tdi_data, include_induction=T, na.rm=T)

    basic_ugt_inhibition_risk_table(p, ugt_inhibition_data)
    basic_ugt_inhibition_risk_table(m, ugt_inhibition_data)

    transporter_inhibition_risk_table(p, transporter_inhibition_data)
    transporter_inhibition_risk_table(m, transporter_inhibition_data)
  })
})
