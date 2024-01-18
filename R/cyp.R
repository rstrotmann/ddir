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
#' @param cyp_inh CYP inhibition data as data frame. The following fields are
#' expected:
#' * 'name' The name of the perpetrator compound.
#' * 'cyp' The CYP enzyme as (upper case) character.
#' * 'ki' The \eqn{k_i} in µM as numeric.
#' * 'source' Optional source information as character.
#' @return A data frame.
#' @seealso [basic_cyp_inhibition_risk_table()]
#' @seealso [key_concentrations()]
#' @seealso [mech_stat_cyp_risk()]
#' @export
#' @examples
#' basic_cyp_inhibition_risk(examplinib_parent, examplinib_cyp_inhibition_data)
basic_cyp_inhibition_risk <- function(perp, cyp_inh) {
  ki <- cyp_inh %>%
    filter(name==name(perp)) #%>%

  i <- key_concentrations(perp, molar=TRUE)
  fumic <- as.num(perp["fumic", "value"])

  out <- ki %>%
    mutate(ki=as.num(ki)) %>%
    mutate(kiu=ki*fumic) %>%
    mutate(r1=1 + (i["imaxssu"]/kiu)) %>%
    mutate(r1gut=case_when(cyp=="CYP3A4" ~ 1+(i["igut"]/ki), .default=NA)) %>%
    select(-name) %>%
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
#' @inheritParams basic_cyp_inhibition_risk
#' @param na.rm Switch to define whether rows with lacking \eqn{K_i} data are
#' removed from the output (i.e., where `ki == NA`). Defaults to `FALSE`.
#' @return A markdown-formatted table, or an empty string.
#' @export
#' @seealso [basic_cyp_inhibition_risk()]
#' @examples
#' basic_cyp_inhibition_risk_table(examplinib_parent, examplinib_cyp_inhibition_data)
#' basic_cyp_inhibition_risk_table(examplinib_parent, examplinib_cyp_inhibition_data, na.rm = TRUE)
basic_cyp_inhibition_risk_table <- function(perp, cyp_inh, na.rm = FALSE) {
  temp <- basic_cyp_inhibition_risk(perp, cyp_inh)

  if(na.rm==TRUE) {
    temp <- temp %>%
      filter(!is.na(ki))
  }

  labels <- c("CYP", "$K_{i}$ (µM)", "$K_{i,u}$ (µM)", "$R_1$",
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
#' This function calculates the risk for time-dependent inhibition of CYP
#' enzymes.
#' @details
#' The relevant metric is \eqn{R_2} as defined in fig. 2 of the
#' [FDA guidance](https://www.fda.gov/media/134582/download):
#'
#' \deqn{R_2=\frac {k_{obs} + k_{deg}}{k_{deg}}}
#'
#' where
#'
#' \deqn{k_{obs}=\frac {50*k_{inact}*I_{max,u}}{K_{I,u} + 50*I_{max,u}}}
#'
#' Values of \eqn{R_2 > 1.25} suggest a relevant TDI potential.
#'
#' The CYP degradation rates, \eqn{k_{deg}} are physiological constants that
#' should be derived from the scientific literature. This package provides
#' standard values for \eqn{k_{deg}} in [cyp_turnover] that are commonly used.
#'
#' In the present version, only the risk for hepatic TDI for CYP enzymes is
#' calculated.
#' @param perp The perpetrator object.
#' @param cyp_tdi The CYP TDI data as data frame. The following fields are expected:
#' * 'name' The perpetrator compound name as character.
#' * 'cyp' The CYP enzyme as character.
#' * 'ki' The \eqn{K_I} in µM as numeric.
#' * 'kinact' The \eqn{k_{inact}} in 1/h as numeric.
#' * 'source' Optional source information as character,
#' @param cyp_kdeg The CYP turnover data as data frame. Defaults to the
#' built-in reference data, [cyp_turnover].
#' @seealso [cyp_turnover]
#' @seealso [read_tdi_data()]
#' @return A data frame.
#' @export
#' @examples
#' basic_cyp_tdi_risk(examplinib_parent, examplinib_cyp_tdi_data)
basic_cyp_tdi_risk <- function(perp, cyp_tdi, cyp_kdeg=cyp_turnover) {
  cyp_tdi <- cyp_tdi %>%
    filter(name==name(perp))

  i <- key_concentrations(perp, molar=TRUE)
  imaxssu <- i["imaxssu"]
  fumic <- as.num(perp["fumic", "value"])
  fu <- as.num(perp["fu", "value"])

  cyp_tdi %>%
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
#' @inheritParams basic_cyp_tdi_risk
#' @param na.rm Switch to define whether rows with lacking \eqn{K_i} data are
#' removed from the output. Defaults to `FALSE`.
#' @return A markdown-formatted table.
#' @seealso [basic_cyp_tdi_risk()]
#' @export
#' @examples
#' basic_cyp_tdi_risk_table(examplinib_parent, examplinib_cyp_tdi_data)
basic_cyp_tdi_risk_table <- function(perp, cyp_tdi, cyp_kdeg=cyp_turnover,
                                     na.rm = TRUE) {
  temp <- basic_cyp_tdi_risk(perp, cyp_tdi, cyp_kdeg)
  if(na.rm==TRUE) {
    temp <- temp %>%
      filter(!is.na(ki) & !is.na(kdeg))
  }

  labels <- c("CYP", "$K_{I} (µM)$", "$f_u$", "$k_{inact}$ (1/h)",
              "$k_{deg}$ (1/h)", "source", "$R_2$", "risk")
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
#' @param perp The perpetrator object or a list thereof.
#' @param cyp_ind The CYP induction data as data frame. The following fields
#' are expected:
#' * 'name' The name of the perpetrator compound as character.
#' * 'cyp' The CYP enzyme as (upper case) character.
#' * 'emax' The \eqn{E_{max}}, i.e., the maximum induction effect determined in
#' vitro as numeric.
#' * 'ec50' The \eqn{EC_{50}} in µM as numeric.
#' * 'maxc' The maximal concentration in µM tested in the in vitro assay as
#' numeric.
#' * 'source' Optional source information as character.
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
#' @inheritParams static_cyp_induction_risk
#' @param na.rm Switch to define whether rows with lacking \eqn{E_{max}} data
#' are removed from the output. Defaults to `FALSE`.
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

  labels <- c("CYP", "$E_{max}$", "$max c$ (µM)", "source", "$max c/I_{max,ss,u}$",
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
#' to fig. 4 of the [FDA guideline](https://www.fda.gov/media/134582/download)):
#'
#' \deqn{R_3 = \frac {1}{1+d* \frac {E_{max}*10*I_{max,u}}{EC_{50} + 10*I_{max,u}}}}
#'
#' \eqn{d} is a scaling factor with a standard value of 1. A different value can
#' be used if warranted by prior experience with the experimental setup.
#'
#' \eqn{R_3 \le 0.8} suggests a relevant in vivo induction potential.
#'
#' @inheritParams static_cyp_induction_risk
#' @param d Scaling factor, defaults to 1.
#' @return A data frame.
#' @export
#' @examples
#' kinetic_cyp_induction_risk(examplinib_parent, examplinib_cyp_induction_data)
kinetic_cyp_induction_risk <- function(perp, cyp_ind, d = 1) {
  i <- key_concentrations(perp, molar=TRUE)

  out <- cyp_ind %>%
    filter(name==name(perp)) %>%
    mutate(r3=(1/(1 + d * emax * 10 * i["imaxssu"] /
                    (ec50 + 10 * i["imaxssu"])))) %>%
    mutate(risk = r3 <= 0.8) %>%
    select(-name)
  return(out)
}


#' Table of the basic kinetc CYP induction risk
#'
#' @inheritParams kinetic_cyp_induction_risk
#' @param na.rm Switch to define whether rows with lacking \eqn{E_{max}} data
#' are removed from the output. Defaults to `FALSE`.
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

  labels <- c("CYP", "$E_{max}$", "$EC_{50}$ (µM)", "$max c$ (µM)",
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
#' fig. 7 of [FDA, 2020](https://www.fda.gov/media/134582/download))
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
#' \eqn{d} is a scaling factor for the CYP induction term with a standard value
#' of 1. A different value can be used if warranted by prior experience with the
#' experimental setup.
#'
#' @inheritParams basic_cyp_inhibition_risk
#' @inheritParams basic_cyp_tdi_risk
#' @inheritParams kinetic_cyp_induction_risk
#' @param include_induction Switch to define whether induction effects should be
#' included in the calculation (C-terms as per the FDA guideline)
#' @param substr The CYP probe substrates to be used as data frame, defaults to
#' [cyp_reference_substrates]. The data frame is expected to have the following
#' fields:
#' * 'cyp' The CYP enzyme as (upper case) character.
#' * 'substrate' The substrate name as character.
#' * 'fgut' The fraction of the drug escaping gut metabolism.
#' * 'fm' The fraction of the drug that undergoes hepatic metabolism.
#' * 'fmcyp' The fraction metabolized by the respective CYP enzyme.
#' @return A data frame.
#' @export
#' @seealso [mech_stat_cyp_risk_table()]
#' @seealso [cyp_reference_substrates]
#' @seealso [cyp_turnover]
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
    d = 1,
    include_induction = TRUE,
    substr = cyp_reference_substrates,
    cyp_kdeg = cyp_turnover) {
  fumic <- as.num(perp["fumic", "value"])

  i <- key_concentrations(perp, molar = TRUE)
  Ig <- i["imaxintest"]
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
    left_join(cyp_kdeg, by = "cyp") %>%
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
    mutate(Cg = case_when((is.na(ec50) | include_induction == FALSE) ~ 1,
                          .default = 1 + (d * emax * Ig / (Ig + ec50)))) %>%
    mutate(Ch = case_when((is.na(ec50) | include_induction == FALSE) ~ 1,
                        .default = 1 + (d * emax * Ih / (Ih + ec50))))  %>%
    select(-name) %>%

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
#' @inheritParams mech_stat_cyp_risk
#' @param na.rm Switch to define whether rows with lacking \eqn{K_I} or
#' \eqn{k_{deg}} data are removed from the output.
#' @return A markdown-formatted table.
#' @seealso [mech_stat_cyp_risk()]
#' @seealso [cyp_reference_substrates]
#' @export
#' @examples
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
    cyp_kdeg = cyp_turnover,
    na.rm = FALSE) {
  temp <- mech_stat_cyp_risk(perp, cyp_inh, cyp_ind, cyp_tdi,
                             include_induction=include_induction,
                             substr=substr, cyp_kdeg=cyp_kdeg) %>%
    select(cyp, kiu, substrate, fgut, fm, fmcyp, Ag, Ah, Bg, Bh, Cg, Ch, aucr,
           risk) %>%
    # mutate(across(Ag:Ch, ~ signif(., digits=2))) %>%
    mutate(across(Ag:Ch, ~ format(., digits=3, nsmall=3))) %>%
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













