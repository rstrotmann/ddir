#### BASIC MODELING

#' Basic CYP inhibition risk
#'
#' This function evaluates the clinical risk for direct (reversible) CYP
#' inhibition according to the basic model defined in the relevant regulatory
#' guidelines.
#'
#' @details
#' For the basic modeling of direct (reversible) CYP enzyme inhibition, the
#' ratio of the relevant inhibitor concentration to the \eqn{K_i} of the
#' respective CYP enzyme is
#' considered, i.e., \eqn{R_1} for hepatic enzymes and \eqn{R_{1,gut}} for
#' intestinal enzymes (refer to fig. 1 of  the the
#' [FDA guidance](https://www.fda.gov/media/134582/download)).
#'
#' ## Liver
#'
#' \deqn{R_1=1+\frac{I_{max,ss,u}}{K_{i,u}}}
#'
#' ## Gut wall
#'
#' \deqn{R_{1,gut}=1+\frac{I_{gut}}{K_{i,u}}}
#'
#' \eqn{R_1} values > 1.02, i.e., maximal unbound perpetrator concentrations
#' 50-fold over \eqn{K_i} are considered to indicate a clinical risk using this
#' method.
#'
#' Note that section 5.3.3.1 of the
#' [EMA DDI guidance](https://www.ema.europa.eu/en/documents/scientific-guideline/guideline-investigation-drug-interactions-revision-1_en.pdf)
#' defines the relevant ratio as \eqn{R_1=\frac{I_{max,ss,u}}{K_{i,u}}} and the
#' threshold as 0.02.
#'
#' In the output, the columns `risk_hep` and `risk_intest` indicate whether the
#' regulatory threshold is reached for the respecive enzyme.
#'
#' Refer to the documentation to the [key_concentrations()] function for details
#' on the calculation of \eqn{I_{max,ss,u}} and \eqn{I_{gut}}.
#'
#' @param perp The perpetrator object.
#' @param cyp_inh CYP inhibition data as data frame.
#'
#' @return A data frame.
#' @seealso [basic_cyp_inhibition_risk_table()]
#' @seealso [key_concentrations()]
#' @seealso [mech_stat_cyp_risk()]
#' @export
#' @examples
#' basic_cyp_inhibition_risk(examplinib_parent, examplinib_cyp_inhibition_data)
#'
basic_cyp_inhibition_risk <- function(perp, cyp_inh) {
  ki <- cyp_inh %>%
    filter(name==name(perp)) #%>%
    # filter(param!="name")

  i <- key_concentrations(perp, molar=TRUE)
  fumic <- as.num(perp["fumic", "value"])

  out <- ki %>%
    mutate(cyp=item) %>%
    mutate(ki=as.num(ki)) %>%
    mutate(kiu=ki*fumic) %>%
    mutate(r1=1 + (i["imaxssu"]/kiu)) %>%
    mutate(r1gut=case_when(cyp=="CYP3A4" ~ 1+(i["igut"]/ki), .default=NA)) %>%
    select(-name, item) %>%
    mutate(risk_hep=r1>1.02) %>%
    mutate(risk_intest=r1gut>11) %>%
    mutate(r1=format(r1, digits=4)) %>%
    mutate(r1gut=format(r1gut, digits=4)) %>%
    select(cyp, ki, kiu, r1, risk_hep, r1gut, risk_intest)

  return(out)
}


#' Basic CYP inhibition risk table
#'
#' This function generates a markdown-formatted table of the direct (reversible)
#' CYP inhibition risk assessment. See [basic_cyp_inhibition_risk()] for details
#' on the calculation of the risk.
#'
#' @param perp The perpetrator object.
#' @param cyp_inh CYP inhibition data as data frame.
#' @param na.rm Boolean to define whether rows with lacking \eqn{K_i} data are
#' removed from the output (i.e., where `ki == NA`).
#' @return A markdown-formatted table, or an empty string (if no CYP inhibition
#' data are available).
#' @export
#' @seealso [basic_cyp_inhibition_risk()]
#' @examples
#' basic_cyp_inhibition_risk_table(examplinib_parent, examplinib_cyp_inhibition_data)
#' basic_cyp_inhibition_risk_table(examplinib_parent, examplinib_cyp_inhibition_data, na.rm = TRUE)
basic_cyp_inhibition_risk_table <- function(perp, cyp_inh, na.rm=F) {
  temp <- basic_cyp_inhibition_risk(perp, cyp_inh)

  if(na.rm==TRUE) {
    temp <- temp %>%
      filter(!is.na(ki))
  }

  labels <- c("CYP", "$K_{i}$", "$K_{i,u}$", "$R_1$",
              "risk (hepatic)", "$R_{1,gut}$", "risk (intestinal)")
  if(nrow(temp)!=0) {
    out <- knitr::kable(temp,
      caption=paste("Risk for direct CYP inhibition by", name(perp), "(basic model)"),
      col.names=labels)
    return(out)
  }
}


#' Basic modeling of the CYP time-dependent inhibition risk
#'
#' In the present version, only the risk for hepatic TDI for CYP enzymes is
#' calculated.
#' @param perp The perpetrator object.
#' @param tdi The CYP TDI data as data frame.
#' @param cyp_kdeg The CYP turnover data as data frame. Defaults to the
#' built-in reference data.
#' @seealso [cyp_turnover]
#' @return A data frame.
#' @export
#'
#' @examples
#' basic_cyp_tdi_risk(examplinib_parent, examplinib_cyp_tdi_data)
basic_cyp_tdi_risk <- function(perp, tdi, cyp_kdeg=cyp_turnover) {
  tdi <- tdi %>%
    filter(name==name(perp))

  i <- key_concentrations(perp, molar=TRUE)
  imaxssu <- i["imaxssu"]
  fumic <- as.num(perp["fumic", "value"])
  fu <- as.num(perp["fu", "value"])

  tdi %>%
    mutate(kobs=kinact*50*imaxssu/(ki * fu + 50 * imaxssu)) %>%
    mutate(fu=fu) %>%
    left_join(cyp_kdeg, by="cyp") %>%
    mutate(kdeg=kdeg_hepatic) %>%
    mutate(r2=(kobs + kdeg)/kdeg) %>%
    mutate(risk=(r2>1.25)) %>%
    select(cyp, ki, fu, kinact, kdeg, source, r2, risk)
}


#' Basic CYP time-dependent inhibition risk table
#'
#' @param perp The perpetrator object.
#' @param tdi The CYP TDI data as data frame.
#' @param cyp_kdeg The CYP turnover data as data frame. Defaults to the
#' built-in reference data.
#' @param na.rm Boolean to define whether rows with lacking \eqn{K_I} or
#' \eqn{k_{deg}} data are removed from the output.
#' @return A markdown-formatted table.
#' @export
#'
#' @examples
#' basic_cyp_tdi_risk_table(examplinib_parent, examplinib_cyp_tdi_data)
basic_cyp_tdi_risk_table <- function(perp, tdi, cyp_kdeg=cyp_turnover,
                                     na.rm = TRUE) {
  temp <- basic_cyp_tdi_risk(perp, tdi, cyp_kdeg)
  if(na.rm==TRUE) {
    temp <- temp %>%
      filter(!is.na(ki) & !is.na(kdeg))
  }

  labels <- c("CYP", "$K_{I}$", "$f_u$", "$k_{inact}$", "$k_{deg}$", "source",
              "$R_2$", "risk")
  if(nrow(temp)!=0) {
    out <- knitr::kable(temp,
      caption=paste("Risk for CYP TDI by", name(perp), "(basic model)"),
      digits = 2,
      col.names=labels)
    return(out)
  }
}


#### CYP INDUCTION

#' Basic static CYP induction risk
#'
#' @details
#' The basic (EMA) or fold-change (FDA) methods evaluate whether the maximal
#' fold-change in mRNA expression is > 2-fold at the expected unbound hepatic
#' concentration of the drug.
#'
#' Regarding the relevant drug concentrations, the FDA guidance suggests
#' considering \eqn{30*I_{max,ss,u}} while the EMA guidance considers
#' \eqn{50*I_{max,ss,u}} for hepatic and \eqn{0.1*I_{gut}} for intestinal
#' CYP enzyme induction. It is expected that the concentrations in the
#' respective in vitro assays cover these concentrations.
#'
#'
#' @param perp The perpetrator object or a list thereof.
#' @param cyp_ind The CYP induction data as data frame.
#'
#' @return A data frame.
#' @export
#' @examples
#' static_cyp_induction_risk(examplinib_parent, examplinib_cyp_induction_data)
static_cyp_induction_risk <- function(perp, cyp_ind)  {
  i <- key_concentrations(perp, molar=TRUE)

  cyp_ind %>%
    filter(name==name(perp)) %>%
    mutate(maxc_imaxssu=round(maxc/i["imaxssu"], 1)) %>%
    mutate(risk=emax>2) %>%
    mutate(note=case_when(
      maxc_imaxssu<30 ~ "Maximal tested concentration is below EMA/FDA expectations",
      (maxc_imaxssu>30 & maxc_imaxssu<50)~"Maximal tested concentration is below FDA expectations",
      .default="")) %>%
    select(-c(name, ec50))
}


#' Table of the basic static CYP induction risk
#'
#' @param perp The perpetrator object.
#' @param cyp_ind The CYP induction data as data frame.
#' @param na.rm Remove rows with lacking ki data (i.e., where ki == NA).
#'
#' @return A markdown-formatted table.
#' @export
#' @seealso [static_cyp_induction_risk()]
#' @examples
#' static_cyp_induction_risk_table(examplinib_parent, examplinib_cyp_induction_data)
static_cyp_induction_risk_table <- function(perp, cyp_ind, na.rm=F) {
  temp <- static_cyp_induction_risk(perp, cyp_ind)

  if(na.rm==TRUE) {
    temp <- temp %>%
      filter(!is.na(emax))
  }

  labels <- c("CYP", "$E_{max}$", "$max c$", "source", "$max c/I_{max,ss,u}$",
              "risk", "notes")

  if(nrow(temp)!=0) {
    out <- knitr::kable(
      temp, caption=paste("Risk for hepatic CYP induction by", name(perp),
        "(basic static model)"),
      col.names=labels)
    return(out)
  }
}


#' Basic kinetic CYP induction risk
#'
#' @details
#' Basic kinetic modeling of the CYP induction risk considers \eqn{R_3} (refer
#' to [FDA 2020](https://www.fda.gov/media/134582/download)):
#'
#' \deqn{R_3 = \frac {1}{1+d* \frac {E_{max}*10*I_{max,u}}{EC_{50} + 10*I_{max,u}}}}
#'
#' For the risk assessment, a threshold of 0.8 applies to \eqn{R_3}.
#'
#' @param perp The perpetrator object.
#' @param cyp_ind The CYP induction data as data frame.
#'
#' @return A data frame.
#' @export
#' @examples
#' kinetic_cyp_induction_risk(examplinib_parent, examplinib_cyp_induction_data)
kinetic_cyp_induction_risk <- function(perp, cyp_ind) {
  i <- key_concentrations(perp, molar=TRUE)

  out <- cyp_ind %>%
    filter(name==name(perp)) %>%
    mutate(r3=(1/(1 + emax*10*i["imaxssu"]/(ec50 + 10*i["imaxssu"])))) %>%
    mutate(risk=r3<=0.8) %>%
    select(-name)
  return(out)
}


#' Table of the basic kinetc CYP induction risk
#'
#' @param perp The perpetrator object.
#' @param cyp_ind The CYP induction data as data frame.
#' @param na.rm Remove rows with lacking ki data (i.e., where ki == NA).
#'
#' @return A markdown-formatted table.
#' @export
#' @seealso [kinetic_cyp_induction_risk()]
#' @examples
#' kinetic_cyp_induction_risk_table(examplinib_parent, examplinib_cyp_induction_data)
kinetic_cyp_induction_risk_table <- function(perp, cyp_ind, na.rm=F) {
  temp <- kinetic_cyp_induction_risk(perp, cyp_ind) %>%
    mutate(r3=format(r3, digit=3))

  if(na.rm==TRUE) {
    temp <- temp %>%
      filter(!is.na(emax))
  }

  labels <- c("CYP", "$E_{max}$", "$EC_{50}$", "$max c$",
              "source", "$R_3$", "risk")
  if(nrow(temp)!=0) {
    out <- knitr::kable(temp,
                        caption=paste("Risk for CYP induction by", name(perp),
                                      "(basic kinetic model)"),
                        col.names=labels)
    return(out)
  }
}


#### MECHANISTIC STATIC MODELING

#' CYP perpetration risk as per mechanistic-static modeling
#'
#' @details
#' In this approach, AUC ratios for probe substrates are calculated based on
#' their known intestinal and hepatic metabolism. Both direct (competitive) and
#' time-dependent inhibition are considered.
#'
#' The below formula given by the FDA guideline (refer to
#' fig. 7 of [FDA 2020](https://www.fda.gov/media/134582/download))
#' also includes intestinal and hepatic enzyme induction terms (\eqn{C_g} and
#' \eqn{C_h}, respectively). At the same time, the guideline states that both
#' inhibition and induction should be considered separately.
#'
#' \deqn{AUCR = \frac{1}{A_g*B_g*C_g* \left(1-F_g \right)+F_g} * \frac{1}{A_h*B_h*C_h*f_m+\left(1-f_m\right)}}
#'
#' The individual terms are:
#'
#' ## Reversible inhibition
#'
#' \deqn{A_g = \frac{1}{1+\frac{I_g}{K_i}}}
#'
#' \deqn{A_h = \frac{1}{1+\frac{I_h}{K_i}}}
#'
#' ## Time-dependent inhibition
#'
#' \deqn{B_g = \frac{k_{deg,g}}{k_{deg,g} + \frac{I_g*k_{inact}}{I_g+K_I}}}
#'
#' \deqn{B_h = \frac{k_{deg,h}}{k_{deg,h} + \frac{I_h*k_{inact}}{I_h+K_I}}}
#'
#' ## Induction
#'
#' \deqn{C_g = 1 + \frac{d*E_{max}*I_g}{I_g+EC_{50}}}
#'
#' \deqn{C_h = 1 + \frac{d*E_{max}*I_h}{I_h+EC_{50}}}
#'
#' with the hepatic inlet concentration \eqn{I_h=I_{max,inlet,u}} and the
#' intestinal concentration \eqn{I_g=I_{enteric,u}}, see
#' [`key_concentrations()`].
#'
#' \eqn{d} is an induction scaling factor (assumed to
#' be 1, but can be adjusted).
#'
#' @param perp The perpetrator object.
#' @param cyp_inh CYP inhibition data as data frame.
#' @param cyp_ind CYP induction data as data frame.
#' @param cyp_tdi CYP TDI data as data frame.
#' @param include_induction Boolean value to define whether induction should be
#'   included in the calculation (C-terms as per the FDA guideline)
#' @param substr The CYP reference substrates to be used as data frame.
#' @param kdeg The CYP turnover data as data frame. Defaults to the
#' built-in reference data.
#' @return A data frame.
#' @export
#' @examples
#' mech_stat_cyp_risk(examplinib_parent, examplinib_cyp_inhibition_data,
#'   examplinib_cyp_induction_data)
#' mech_stat_cyp_risk(examplinib_parent, examplinib_cyp_inhibition_data,
#'   examplinib_cyp_induction_data, examplinib_cyp_tdi_data)
#' mech_stat_cyp_risk(examplinib_parent, examplinib_cyp_inhibition_data,
#'   examplinib_cyp_induction_data, examplinib_cyp_tdi_data,
#'   include_induction = FALSE)
mech_stat_cyp_risk <- function(
    perp,
    cyp_inh,
    cyp_ind,
    cyp_tdi = NULL,
    include_induction = TRUE,
    substr=cyp_reference_substrates,
    kdeg=cyp_turnover) {
  fumic <- as.num(perp["fumic", "value"])

  i <- key_concentrations(perp, molar = TRUE)
  Ig <- i["imaxintestu"]
  Ih <- i["imaxinletu"]

  if(is.null(cyp_tdi)) {
    cyp_tdi = data.frame(
      name="",
      cyp = "",
      ki = NA,
      kinact = NA,
      source = "")}

  out <- cyp_inh %>%
    filter(name == name(perp)) %>%
    mutate(cyp = item) %>%
    select(-source) %>%

    # direct inhibition
    mutate(ki = as.num(ki)) %>%
    mutate(kiu = ki * fumic) %>%
    mutate(Ag = case_when(!is.na(ki) ~ 1 / (1 + (Ig / kiu)),
                          .default = 1)) %>%
    mutate(Ah = case_when(!is.na(ki) ~ 1 / (1 + (Ih / kiu)),
                          .default = 1)) %>%

    # TDI
    left_join(
      cyp_tdi %>%
        filter(name == name(perp)) %>%
        mutate(ki_tdi = ki) %>%
        select(-ki, -name, -source),
      by = "cyp") %>%
    left_join(kdeg, by = "cyp") %>%
    mutate(Bg = case_when(
      !is.na(ki_tdi) ~ kdeg_intestinal / (kdeg_intestinal +
                                            (Ig * kinact /(Ig + ki_tdi))),
      .default = 1)) %>%
    mutate(Bh = case_when(
      !is.na(ki_tdi) ~ kdeg_hepatic / (kdeg_hepatic +
                                         (Ih * kinact / (Ih + ki_tdi))),
      .default = 1)) %>%

    # induction
    left_join(cyp_ind %>%
                filter(name == name(perp)) %>%
                select(-source, -name),
              by=c("cyp")) %>%
    # mutate(Cg = 1, Ch = 1) %>% # temporary
    # mutate(Cg = 1 + (emax * Ig / (Ig + ec50))) %>%
    mutate(Cg = case_when((is.na(ec50) | include_induction == FALSE) ~ 1,
                          .default = 1 + (emax * Ig / (Ig + ec50)))) %>%
    # mutate(Ch = 1 + (emax * Ih / (Ih + ec50))) %>%
    mutate(Ch = case_when((is.na(ec50) | include_induction == FALSE) ~ 1,
                        .default = 1 + (emax * Ih / (Ih + ec50))))  %>%
    select(-name, -item) %>%

    # substrate
    left_join(substr, by="cyp") %>%
    mutate(aucr = 1 / (Ag * Bg * Cg * (1 - fgut) + fgut) *
             1 / (Ah * Bh * Ch * fm * fmcyp + (1 - fm * fmcyp))) %>%
    mutate(risk = aucr>1.25 | aucr < 0.8) %>%
    select(cyp, substrate, kiu, fgut, fm, fmcyp, Ag, Ah, Bg, Bh, Cg, Ch, aucr,
           risk) %>%
    as.data.frame()

  return(out)
}


#' CYP perpetration risk table as per mechanistic-static modeling
#'
#' @param perp The perpetrator object.
#' @param cyp_inh CYP inhibiton data as data frame.
#' @param cyp_ind CYP induction data as data frame.
#' @param include_induction Boolean value to define whether induction should be
#'   included in the calculation (C-terms as per the FDA guideline)
#' @param substr The CYP reference substrates to be used as data frame.
#' @param na.rm Remove rows with lacking ki data (i.e., where ki == NA).
#' @param cyp_tdi CYP TDI data as data frame.
#' @param kdeg The CYP turnover data as data frame. Defaults to the
#' built-in reference data.
#' @return A markdown-formatted table.
#' @seealso [mech_stat_cyp_risk()]
#' @export
#' @examples
#'
#' mech_stat_cyp_risk_table(
#'   examplinib_parent,
#'   examplinib_cyp_inhibition_data,
#'   examplinib_cyp_induction_data)
#' mech_stat_cyp_risk_table(
#'   examplinib_parent,
#'   examplinib_cyp_inhibition_data,
#'   examplinib_cyp_induction_data,
#'   examplinib_cyp_tdi_data)
mech_stat_cyp_risk_table <- function(
    perp,
    cyp_inh,
    cyp_ind,
    cyp_tdi = NULL,
    include_induction = TRUE,
    substr = cyp_reference_substrates,
    kdeg = cyp_turnover,
    na.rm = FALSE) {
  temp <- mech_stat_cyp_risk(perp, cyp_inh, cyp_ind, cyp_tdi,
                             include_induction=include_induction,
                             substr=substr, kdeg=kdeg) %>%
    select(cyp, kiu, substrate, fgut, fm, fmcyp, Ag, Ah, Bg, Bh, Cg, Ch, aucr,
           risk) %>%
    # mutate(across(Ag:Ch, ~ format(.x, signif=2))) %>%
    mutate(across(Ag:Ch, ~ signif(., digits=2))) %>%
    mutate(aucr=format(aucr, digits=3))

  if(na.rm==TRUE) {
    temp <- temp %>%
      filter(!is.na(risk))
  }

  labels <- c("CYP", "$K_{i,u}$", "substrate", "$F_{gut}$", "$f_m$",
              "$f_{m,CYP}$", "$A_g$", "$A_h$", "$B_g$", "$B_h$", "$C_g$",
              "$C_h$", "AUCR", "risk")

  if(nrow(temp)!=0) {
    out <- knitr::kable(
      temp,
      caption=paste("Mechanistic static modeling of CYP inhibition risk for ",
                    name(perp)),
      col.names=labels)
    return(out)
  }
}













