#' Drug transporter inhibition risk
#'
#' @details
#' The metric for the assessment of transporter interactions is
#' \eqn{R=[I]/IC_{50}}.
#'
#' The relevant perpetrator concentrations \eqn{[I]} and regulatory thresholds
#' of concern are:
#'
#' | I                     | transporter                                                                    | threshold |
#' | ---                   | ---                                                                            | ---       |
#' | \eqn{I_{gut}}         | P-gp and BRCR when drugs are orally administered                               | 10        |
#' | \eqn{C_{max,ss,u}}    | P-gp and BRCR when drugs are administered parenterally or for drug metabolites | 0.02      |
#' | \eqn{I_{max,inlet,u}} | hepatic basolateral transporters OCT1, OATP1B1 and OATP1B3                     | 0.1       |
#' | \eqn{C_{max,ss,u}}    | renal basolateral transporters OAT1, OAT3 and OCT2                             | 0.1       |
#' | \eqn{C_{max,ss,u}}    | apical transporters MATE1 and MATE2-K                                          | 0.02      |
#'
#' @param perp The perpetrator object.
#' @param transporter_inh Transporter inhibition data as data frame. The
#' following fields are expected:
#' * 'name' The perpetrator compound name.
#' * 'cyp' The UGT enzyme as (upper case) character.
#' * 'ic50' The \eqn{IC_{50}} of the inhibition in μM.
#' * 'source' Optional source information as character.
#' @param transporter_ref Transporter reference data, see
#' [transporter_reference_data] for details.
#' @return A data frame.
#' @seealso [transporter_inhibition_risk_table()]
#' @seealso [read_transporter_inhibitor_data()]
#' @export
#' @examples
#' transporter_inhibition_risk(examplinib_parent,
#'   examplinib_transporter_inhibition_data)
transporter_inhibition_risk <- function(
    perp,
    transporter_inh,
    transporter_ref = transporter_reference_data) {
  ic50 <- transporter_inh %>%
    filter(name == name(perp)) %>%
    select(-name) %>%
    mutate(ic50 = as.num(ic50))

  temp <-  key_concentrations(perp, molar=TRUE)
  i <- data.frame(
    i=names(temp),
    conc=temp)

  # duplicate rows Pgp and BCRP, if applicable, and assign intestinal and
  #   systemic scope
  out <- ic50 %>%
    bind_rows(filter(ic50, transporter %in% c("Pgp", "BCRP")) %>%
                mutate(transporter=paste0(transporter, "_sys"))) %>%
    bind_rows(filter(ic50, transporter %in% c("Pgp", "BCRP")) %>%
                mutate(transporter=paste0(transporter, "_int"))) %>%
    filter(!transporter %in% c("Pgp", "BCRP")) %>%
    left_join(transporter_ref %>% mutate(rank = row_number()),
              by="transporter") %>%
    left_join(i, by="i") %>%
    mutate(r=case_when(is.na(ic50) ~ NA, .default=conc/ic50)) %>%
    mutate(risk=r>threshold) %>%
    arrange(rank) %>%
    select(transporter, ic50, source, i, r, threshold, risk)
  return(out)
}


#' Table of drug transporter inhibition risks
#'
#' @inheritParams transporter_inhibition_risk
#' @param na.rm Exclude rows with lacking \eqn{IC_{50}} data from the
#' output, as logical.
#' @param show_dose Show_dose Show perpetrator dose in table title, as logical.
#' Defaults to `FALSE.`
#' @seealso [transporter_inhibition_risk()]
#' @return A markdown-formatted table.
#' @export
#' @examples
#' transporter_inhibition_risk_table(examplinib_parent,
#'   examplinib_transporter_inhibition_data)
#' transporter_inhibition_risk_table(examplinib_compounds,
#'   examplinib_transporter_inhibition_data)
#' transporter_inhibition_risk_table(examplinib_compounds,
#'   examplinib_transporter_inhibition_data, show_dose = TRUE)
transporter_inhibition_risk_table <- function(
    perp,
    transporter_inh,
    transporter_ref=transporter_reference_data,
    na.rm=F,
    show_dose = FALSE) {
  UseMethod("transporter_inhibition_risk_table")
}


#' Table of drug transporter inhibition risks
#' @inheritParams transporter_inhibition_risk_table
#' @export
#' @noRd
transporter_inhibition_risk_table.perpetrator <- function(
    perp,
    transporter_inh,
    transporter_ref=transporter_reference_data,
    na.rm=F,
    show_dose = FALSE) {
  i_names <- tribble(
    ~i,            ~i_label,
    "igut",       "$I_{gut}$",
    "imaxssu",    "$C_{max,ss,u}$",
    "imaxinletu", "$I_{max,inlet,u}$",
    "imaxintest", "$I_{max,intest,u$"
  )
  temp <- transporter_inhibition_risk(
    perp, transporter_inh,
    transporter_ref=transporter_reference_data) %>%
    mutate(risk = case_match(
      as.character(risk),
      "TRUE" ~ "Yes", "FALSE" ~ "No", .default = "")) %>%
    mutate(r=round(r, 2)) %>%
    # left_join(i_names, by = "i") %>%
    # select(transporter, ic50, source, i_label, r, threshold, risk)
    select(transporter, ic50, source, r, threshold, risk)


  if(na.rm==TRUE) {
    temp <- temp %>%
      filter(!is.na(ic50))
  }

  # labels <- c("transporter", "$IC_{50}$", "source", "$I$", "$R$", "threshold",
  #             "risk")
  labels <- c("transporter", "$IC_{50}$", "source", "$R$", "threshold", "risk")
  if(nrow(temp)!=0) {
    caption <- paste0("Risk for drug transporter inhibition by ",
                     name(perp), conditional_dose_string(perp, show_dose))
    out <- knitr::kable(temp, caption = caption, col.names=labels, signif=2)
    return(out)
  }
}


#' Table of drug transporter inhibition risks
#' @inheritParams transporter_inhibition_risk_table
#' @export
#' @noRd
transporter_inhibition_risk_table.list <- function(perp, ...) {
  for(i in perp) {
    temp <- transporter_inhibition_risk_table.perpetrator( i, ...)
    if(!is.null(temp)) {
      print(temp)
    }
  }
}
