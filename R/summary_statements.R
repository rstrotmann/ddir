
#' Title
#'
#' @inheritParams cyp_inhibition_risk_summary
#' @inheritParams cyp_tdi_risk_summary
#' @inheritParams cyp_induction_risk_summary
#' @inheritParams ugt_inhibition_risk_summary
#'
#' @return Text output.
#' @export
#'
#' @examples
#' ddi_risk_summary(examplinib_compounds, examplinib_cyp_inhibition_data,
#'   examplinib_cyp_tdi_data, examplinib_cyp_induction_data,
#'   examplinib_ugt_inhibition_data, examplinib_transporter_inhibition_data)
#' ddi_risk_summary(examplinib_parent, examplinib_cyp_inhibition_data)
ddi_risk_summary <- function(
    perp,
    cyp_inh = NULL,
    cyp_tdi = NULL,
    cyp_ind = NULL,
    ugt_inh = NULL,
    transporter_inh = NULL,
    d = 1,
    include_induction = TRUE,
    substr = cyp_reference_substrates,
    cyp_kdeg = cyp_turnover) {
  out <- c(
    cyp_inhibition_risk_summary(perp, cyp_inh),
    cyp_tdi_risk_summary(perp, cyp_tdi),
    cyp_induction_risk_summary(perp, cyp_ind, d),
    mech_stat_cyp_risk_summary(perp, cyp_inh, cyp_ind, cyp_tdi, d,
                               include_induction, substr, cyp_kdeg),
    ugt_inhibition_risk_summary(perp, ugt_inh),
    transporter_inhibition_risk_summary(perp, transporter_inh)
  )
  out <- paste(out, collapse = "\n")
  return(out)
}


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
cyp_inhibition_risk_summary.perpetrator <- function(perp, cyp_inh) {
  temp <- basic_cyp_inhibition_risk(perp, cyp_inh) %>%
    filter(risk_hep == TRUE)

  if(nrow(temp) > 0){
    out <- paste0(
      "* using basic modeling, ", name(perp),
      " has a clinical risk for direct inhibition of ",
      nice_enumeration(temp$cyp))
  } else {
    out <- NULL
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
cyp_inhibition_risk_summary.list <- function(perp, ...) {
  out <- NULL
  for(i in perp) {
    temp <- cyp_inhibition_risk_summary.perpetrator(i, ...)
    out <- c(out, temp)
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
#' @inheritParams basic_cyp_tdi_risk
#'
#' @return Text output.
#' @export
cyp_tdi_risk_summary.perpetrator <- function(perp, cyp_tdi, cyp_kdeg=cyp_turnover) {
  out <- NULL
  if(!is.null(cyp_tdi)){
    temp <- basic_cyp_tdi_risk(perp, cyp_tdi, cyp_kdeg) %>%
      filter(risk == TRUE)

    if(nrow(temp) > 0){
      out <- paste0("* ", name(perp),
                    " has a clinical risk for time-dependent inhibition of ",
                    nice_enumeration(temp$cyp))
    }
  }
  return(out)
}


#' Summary statement for CYP TDI
#'
#' @param perp A list of perpetrator objects.
#' @param ... Further parameters.
#'
#' @return Text output.
#' @export
#' @examples
#' cyp_tdi_risk_summary(examplinib_compounds, examplinib_cyp_tdi_data)
cyp_tdi_risk_summary.list <- function(perp, ...) {
  out <- NULL
  for(i in perp) {
    temp <- cyp_tdi_risk_summary.perpetrator(i, ...)
    out <- c(out, temp)
  }
  return(out)
}


### ----- CYP INDUCTION -----

#' Title
#'
#' @inheritParams static_cyp_induction_risk
#' @inheritParams kinetic_cyp_induction_risk
#'
#' @return Text output.
#' @export
cyp_induction_risk_summary <- function(perp, cyp_ind, d = 1) {
  UseMethod("cyp_induction_risk_summary")
}


#' Title
#'
#' @inheritParams cyp_induction_risk_summary
#' @return Text output.
#' @export
cyp_induction_risk_summary.perpetrator <- function(perp, cyp_ind, d = 1) {
  out <- NULL
  if(!is.null(cyp_ind)){
    temp <- static_cyp_induction_risk(perp, cyp_ind) %>%
      filter(risk == TRUE)
    if(nrow(temp) > 0){
      out <- paste0(
        "* ", name(perp),
        " has a clinical risk for induction of ", nice_enumeration(temp$cyp),
        " (fold-change method)")
    }

    temp <- kinetic_cyp_induction_risk(perp, cyp_ind, d = d) %>%
      filter(risk == TRUE)
    if(nrow(temp) > 0){
      out <- c(
        out, paste0(
        "* ", name(perp),
        " has a clinical risk for induction of ", nice_enumeration(temp$cyp),
        " (basic kinetic method)"))
    }
  }

  return(out)
}


#' Title
#'
#' @param perp A list of perpetrator objects.
#' @param ... Further arguments.
#'
#' @return Text output.
#' @export
cyp_induction_risk_summary.list <- function(perp, ...) {
  out <- NULL
  for(i in perp) {
    temp <- cyp_induction_risk_summary.perpetrator(i, ...)
    out <- c(out, temp)
  }
  return(out)
}


### ----- MECHANISTIC STATIC MODELING -----


#' Title
#'
#' @inheritParams mech_stat_cyp_risk
#' @return Text ouput.
#' @export
#'
#' @examples
#' mech_stat_cyp_risk_summary(examplinib_parent,
#'   examplinib_cyp_inhibition_data, examplinib_cyp_induction_data)
#' mech_stat_cyp_risk_summary(examplinib_metabolite,
#'   examplinib_cyp_inhibition_data, examplinib_cyp_induction_data)
#' mech_stat_cyp_risk_summary(examplinib_compounds,
#'   examplinib_cyp_inhibition_data, examplinib_cyp_induction_data)
mech_stat_cyp_risk_summary <- function(
    perp,
    cyp_inh,
    cyp_ind,
    cyp_tdi = NULL,
    d = 1,
    include_induction = TRUE,
    substr = cyp_reference_substrates,
    cyp_kdeg = cyp_turnover){
  UseMethod("mech_stat_cyp_risk_summary")
}


#' Title
#'
#' @inheritParams mech_stat_cyp_risk
#'
#' @return Text output.
#' @export
mech_stat_cyp_risk_summary.perpetrator <- function(
    perp,
    cyp_inh,
    cyp_ind = NULL,
    cyp_tdi = NULL,
    d = 1,
    include_induction = TRUE,
    substr = cyp_reference_substrates,
    cyp_kdeg = cyp_turnover) {
  out <- NULL
  if(!is.null(cyp_inh)){
    temp <- mech_stat_cyp_risk(
      perp,
      cyp_inh,
      cyp_ind,
      cyp_tdi = NULL,
      d = 1,
      include_induction = TRUE,
      substr = cyp_reference_substrates,
      cyp_kdeg = cyp_turnover) %>%
      filter(risk == TRUE)

    inhib <- filter(temp, aucr > 1)
    induct <- filter(temp, aucr < 1)

    if(nrow(inhib) > 0){
      out <- paste0(
        "* ", "based on mechanistic-static modeling (",
        nice_enumeration(inhib$substrate), "), ",
        name(perp),
        " has a clinical risk for inhibition of ", nice_enumeration(inhib$cyp))
    }
    if(nrow(induct) > 0){
      out <- c(
        out,
        paste0(
          "* ", "based on mechanistic-static modeling (",
          nice_enumeration(induct$substrate), "), ",
          name(perp),
          " has a clinical risk for induction of ", nice_enumeration(induct$cyp))
      )
    }
  }
  return(out)
}


#' Title
#'
#' @param perp A perpetrator object.
#' @param ... Further parameters.
#'
#' @return Text output.
#' @export
mech_stat_cyp_risk_summary.list <- function(perp, ...) {
  out <- NULL
  for(i in perp) {
    temp <- mech_stat_cyp_risk_summary.perpetrator(i, ...)
    out <- c(out, temp)
  }
  return(out)
}


### ----- UGT INHIBITION -----

#' Title
#'
#' @inheritParams basic_ugt_inhibition_risk
#'
#' @return Text output.
#' @export
#' @examples
#' ugt_inhibition_risk_summary(
#'   examplinib_parent, examplinib_ugt_inhibition_data)
#' ugt_inhibition_risk_summary(
#'   examplinib_compounds, examplinib_ugt_inhibition_data)
ugt_inhibition_risk_summary <- function(perp, ugt_inh){
  UseMethod("ugt_inhibition_risk_summary")
}


#' Title
#'
#' @inheritParams ugt_inhibition_risk_summary
#' @return Text output.
#' @export
ugt_inhibition_risk_summary.perpetrator <- function(perp, ugt_inh) {
  out <- NULL
  if(!is.null(ugt_inh)) {
    temp <- basic_ugt_inhibition_risk(perp, ugt_inh) %>%
      filter(risk == TRUE)
    if(nrow(temp) > 0){
      out <- paste0(
        "* ", name(perp),
        " has a clinical risk for inhibition of ", nice_enumeration(temp$ugt))
    }
  }

  return(out)
}


#' Title
#'
#' @param perp A perpetrator object.
#' @param ... Further parameters.
#'
#' @return Text output.
#' @export
ugt_inhibition_risk_summary.list <- function(perp, ...) {
  out <- NULL
  for(i in perp) {
    temp <- ugt_inhibition_risk_summary.perpetrator(i, ...)
    out <- c(out, temp)
  }
  return(out)
}


### ----- TRANSPORTER INHIBITION -----

#' Title
#'
#' @inheritParams transporter_inhibition_risk
#'
#' @return Text output.
#' @export
#' @examples
#' transporter_inhibition_risk_summary(examplinib_parent,
#'   examplinib_transporter_inhibition_data)
#' transporter_inhibition_risk_summary(examplinib_metabolite,
#'   examplinib_transporter_inhibition_data)
#' transporter_inhibition_risk_summary(examplinib_compounds,
#'   examplinib_transporter_inhibition_data)
transporter_inhibition_risk_summary <- function
(perp, transporter_inh, transporter_ref = transporter_reference_data){
  UseMethod("transporter_inhibition_risk_summary")
}


#' Title
#'
#' @inheritParams transporter_inhibition_risk_summary
#'
#' @return Text output
#' @export
transporter_inhibition_risk_summary.perpetrator <- function(
    perp, transporter_inh, transporter_ref = transporter_reference_data) {
  out <- NULL
  if(!is.null(transporter_inh)) {
    temp <- transporter_inhibition_risk(perp, transporter_inh, transporter_ref) %>%
      mutate(transporter = case_match(transporter,
        "Pgp_int" ~ "P-gp (intestinal)", "Pgp_sys" ~ "P-gp (hepatic)",
        "BCRP_int" ~ "BCRP (intestinal)", "BCRP_sys" ~ "BCRP (hepatic)",
        .default = transporter)) %>%
      filter(risk == TRUE)
    if(nrow(temp) > 0){
      out <- paste0(
        "* ", name(perp),
        " has a clinical risk for inhibition of ", nice_enumeration(temp$transporter))
    }
  }

  return(out)
}


#' Title
#'
#' @param perp A perpetrator object.
#' @param ... Further parameters.
#'
#' @return Text output.
#' @export
transporter_inhibition_risk_summary.list <- function(perp, ...) {
  out <- NULL
  for(i in perp) {
    temp <- transporter_inhibition_risk_summary.perpetrator(i, ...)
    out <- c(out, temp)
  }
  return(out)
}


