### ----- DIRECT CYP INHIBITION -----

#' Summary statement for CYP inhibition
#'
#' @inheritParams basic_cyp_inhibition_risk
#'
#' @return Text output.
#' @export
#'
#' @examples
#' cyp_inhibition_risk_summary(examplinib_parent, examplinib_cyp_inhibition_data)
#' cyp_inhibition_risk_summary(examplinib_compounds, examplinib_cyp_inhibition_data)
cyp_inhibition_risk_summary <- function(perp, cyp_inh) {
  UseMethod("cyp_inhibition_risk_summary")
}


#' Summary statement for CYP inhibition
#'
#' @inheritParams cyp_inhibition_risk_summary
#'
#' @return Text output.
#' @export
#'
#' @examples
#' cyp_inhibition_risk_summary.perpetrator(examplinib_parent, examplinib_cyp_inhibition_data)
#' cyp_inhibition_risk_summary.perpetrator(examplinib_metabolite, examplinib_cyp_inhibition_data)
cyp_inhibition_risk_summary.perpetrator <- function(perp, cyp_inh) {
  temp <- basic_cyp_inhibition_risk(perp, cyp_inh) %>%
    filter(risk_hep == TRUE)

  if(nrow(temp) > 0){
    out <- paste0("* ", name(perp),
                  " has a clinical risk for direct inhibition of ",
                  nice_enumeration(temp$cyp), "\n")
  } else {
    out <- ""
  }
  return(out)
}


#' Summary statement for CYP inhibition for a list of perpetrators
#'
#' @param perp Perpetrator object.
#' @param ... Further parameters.
#'
#' @return Text output.
#' @export
#'
#' @examples
#' cyp_inhibition_risk_summary.list(examplinib_compounds, examplinib_cyp_inhibition_data)
cyp_inhibition_risk_summary.list <- function(perp, ...) {
  out <- ""
  for(i in perp) {
    temp <- cyp_inhibition_risk_summary.perpetrator(i, ...)
    out <- paste0(out, temp)
  }
  return(out)
}


### ----- TIME-DEPENDENT CYP INHIBITION -----

#' Summary statement for CYP TDI
#'
#' @inheritParams basic_cyp_tdi_risk
#' @return Text.
#' @export
#'
#' @examples
#' cyp_tdi_risk_summary(examplinib_parent, examplinib_cyp_tdi_data)
#' cyp_tdi_risk_summary(examplinib_compounds, examplinib_cyp_tdi_data)
cyp_tdi_risk_summary <- function(perp, cyp_tdi, cyp_kdeg=cyp_turnover) {
  UseMethod("cyp_tdi_risk_summary")
}


#' Summary statement for CYP TDI
#'
#' @inheritParams basic_cyp_inhibition_risk
#'
#' @return Text output.
#' @export
#'
#' @examples
#' cyp_tdi_risk_summary.perpetrator(examplinib_parent, examplinib_cyp_tdi_data)
cyp_tdi_risk_summary.perpetrator <- function(perp, cyp_tdi, cyp_kdeg=cyp_turnover) {
  temp <- basic_cyp_tdi_risk(perp, cyp_tdi, cyp_kdeg) %>%
    filter(risk == TRUE)

  if(nrow(temp) > 0){
    out <- paste0("* ", name(perp),
                  " has a clinical risk for time-dependent inhibition of ",
                  nice_enumeration(temp$cyp), "\n")
  } else {
    out <- ""
  }
  return(out)
}


#' Summary statement for CYP TDI
#'
#' @inheritParams cyp_tdi_risk_summary
#'
#' @return Text output.
#' @export
#' @examples
#' cyp_tdi_risk_summary(examplinib_compounds, examplinib_cyp_tdi_data)
cyp_tdi_risk_summary.list <- function(perp, ...) {
  out <- ""
  for(i in perp) {
    temp <- cyp_tdi_risk_summary.perpetrator(i, ...)
    out <- paste0(out, temp)
  }
  return(out)
}
