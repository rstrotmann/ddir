
#' Drug transporter inhibition risk
#'
#' @details
#' The relevant metric for the assessment of transporter interactions is
#' \eqn{R=[I]/K_i}. In in vitro transporter inhibition studies, \eqn{IC_{50}}
#' values are experimentally determined. Since the transporter substrate
#' concentration is usually kept very low in relation to \eqn{K_m} to
#' minimize passive permeation, \eqn{K_i = IC_{50}} can be assumed.
#'
#' The relevant perpetrator concentrations \eqn{[I]} are:
#' * \eqn{I_{gut}} for intestinal P-gp and BRCR
#' * \eqn{I_{max,inlet,u}} for the hepatic basolateral transporters OCT1,
#' OATP1B1 and OATP1B3
#' * \eqn{I_{max,ss,u}} for the renal basolateral transporters OAT1, OAT3 and
#' OCT2, as well as the apical transporters outside the intestinal mucosa, i.e.,
#' hepatic P-gp and BCRP, and MATE1, MATE2-k.
#'
#' Note that the FDA and EMA guidelines differ in their definitions for the
#' thresholds for assumed clinically relevant effects.
#'
#' @param perp The perpetrator object.
#' @param transporter_inh Transporter inhibition data as data frame.
#' @param transporter_ref Transporter reference data.
#'
#' @return A data frame.
#' @export
#' @examples
#' transporter_inhibition_risk(examplinib_parent, examplinib_transporter_inhibition_data)
transporter_inhibition_risk <- function(
    perp,
    transporter_inh,
    transporter_ref = transporter_reference_data) {
  ic50 <- transporter_inh %>%
    filter(name==name(perp)) %>%
    select(-name) %>%
    mutate(ic50=as.num(ic50))

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
    left_join(transporter_ref, by="transporter") %>%
    left_join(i, by="i") %>%
    mutate(r=case_when(is.na(ic50) ~ NA, .default=conc/ic50)) %>%
    mutate(fda_risk=r>fda_thld) %>%
    mutate(ema_risk=r>ema_thld) %>%
    arrange(rank) %>%
    select(transporter, ic50, source, r, fda_thld, fda_risk, ema_thld, ema_risk)
  return(out)
}


#' Table of drug transporter inhibition risks
#'
#' @param perp The perpetrator object.
#' @param transporter_ic50 Transporter inhibition IC50 as data frame.
#' @param na.rm Boolean to indicate whether rows with lacking ki data should be
#' removed from the output (i.e., where ki == NA).
#' @param transporter_ref Transporter reference data.
#' @seealso [transporter_inhibition_risk()]
#' @return A markdown-formatted table.
#' @export
#' @examples
#' transporter_inhibition_risk_table(examplinib_parent, examplinib_transporter_inhibition_data)
transporter_inhibition_risk_table <- function(
    perp,
    transporter_ic50,
    transporter_ref=transporter_reference_data,
    na.rm=F) {
  temp <- transporter_inhibition_risk(
      perp, transporter_ic50,
      transporter_ref=transporter_reference_data) %>%
    mutate(r=round(r, 3))

if(na.rm==TRUE) {
    temp <- temp %>%
      filter(!is.na(ic50))
  }

  labels <- c("transporter", "$IC_{50}$", "source", "$R$", "thld FDA", "risk FDA",
              "thld EMA", "risk EMA")
  if(nrow(temp)!=0) {
    out <- knitr::kable(
      temp,
      caption=paste("Risk for drug transporter inhibition by",
                    name(perp)), col.names=labels,
      signif=2)
    return(out)
  }
}
