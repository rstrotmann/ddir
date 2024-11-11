test_that("CYP inhibition summary", {
  expect_no_error(
    invisible(capture_output({
      cyp_inhibition_risk_summary.perpetrator(
        examplinib_parent, examplinib_cyp_inhibition_data)

      cyp_inhibition_risk_summary.list(
        examplinib_compounds, examplinib_cyp_inhibition_data)

      cyp_inhibition_risk_summary(
        examplinib_compounds, examplinib_cyp_inhibition_data)

      static_cyp_induction_risk(
        examplinib_parent, examplinib_cyp_induction_data)
    }))
  )
})

