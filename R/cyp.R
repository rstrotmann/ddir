#### BASIC MODELING

#' Basic CYP inhibition risk
#'
#' This function evaluates the clinical risk for direct (reversible) CYP
#' inhibition according to the basic model defined in the relevant
#' [regulatory guideline](https://www.ema.europa.eu/en/documents/scientific-guideline/ich-m12-guideline-drug-interaction-studies-step-5_en.pdf).
#'
#' @details For the basic modeling of direct (reversible) CYP enzyme inhibition,
#'   the ratio of the relevant inhibitor concentration to the \eqn{K_i} of the
#'   respective CYP enzyme is considered, i.e., \eqn{R} for hepatic enzymes and
#'   \eqn{R_{gut}} for intestinal enzymes (refer to Section 2.1.2.1 of the [ICH
#'   M12 guidance
#'   document](https://www.ema.europa.eu/en/documents/scientific-guideline/ich-m12-guideline-drug-interaction-studies-step-5_en.pdf)).
#'
#'   ## Liver
#'
#'   \deqn{R=\frac{C_{max,ss,u}}{K_{i,u}}}
#'
#'   \eqn{R} values > 0.02, i.e., maximal unbound perpetrator concentrations
#'   50-fold over \eqn{K_i} are considered to indicate a potential clinical CYP
#'   inhibition risk using this method.
#'
#'   ## Intestine
#'
#'   \deqn{R_{gut}=\frac{I_{gut}}{K_{i,u}}}
#'
#'   where
#'
#'   \deqn{I_{gut}=\frac{Dose}{250\ mg}}
#'
#'   \eqn{R} values > 10 are considered to indicate a clinical risk for
#'   intestinal CYP3A inhibition.
#'
#'   In the output, the columns `risk_hep` and `risk_intest` indicate whether
#'   the regulatory threshold is reached for the respective enzyme.
#'
#' @param perp The perpetrator object.
#' @param cyp_inh CYP inhibition data as data frame. The following fields are
#'   expected:
#' * 'name' The name of the perpetrator compound.
#' * 'cyp' The CYP enzyme as (upper case) character.
#' * 'ki' The \eqn{K_i} in \eqn{\mu M} as numeric.
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
    filter(name==name(perp))

  i <- key_concentrations(perp, molar=TRUE)
  fumic <- as.num(perp["fumic", "value"])

  out <- ki %>%
    mutate(ki=as.num(ki)) %>%
    mutate(kiu=ki*fumic) %>%
    mutate(r=i["imaxssu"]/kiu) %>%
    mutate(r_gut=case_when(cyp=="CYP3A4" ~ i["igut"]/kiu, .default=NA)) %>%
    select(-name) %>%
    mutate(risk_hep=r>0.02) %>%
    mutate(risk_intest=r_gut>10) %>%
    # mutate(r=format(r, digits=4)) %>%
    mutate(r=round(r, digits=4)) %>%
    # mutate(r_gut=format(r_gut, digits=4)) %>%
    mutate(r_gut=round(r_gut, digits=4)) %>%
    select(cyp, ki, kiu, r, risk_hep, r_gut, risk_intest)

  return(out)
}


#' Basic CYP inhibition risk table
#'
#' This function generates a markdown-formatted table of the direct (reversible)
#' CYP inhibition risk assessment. See [basic_cyp_inhibition_risk()] for details
#' on the calculation of the risk.
#' @inheritParams basic_cyp_inhibition_risk
#' @param na.rm Switch to define whether rows with lacking \eqn{K_i} data are
#' removed from the output (i.e., where `ki == NA`). Defaults to `FALSE`.
#' @param show_dose Show perpetrator dose in table title, defaults to `FALSE.`
#' @return A markdown-formatted table, or an empty string.
#' @export
#' @seealso [basic_cyp_inhibition_risk()]
#' @examples
#' basic_cyp_inhibition_risk_table(examplinib_parent,
#'   examplinib_cyp_inhibition_data)
#' basic_cyp_inhibition_risk_table(examplinib_compounds,
#'   examplinib_cyp_inhibition_data, na.rm = TRUE)
#' basic_cyp_inhibition_risk_table(examplinib_compounds,
#'   examplinib_cyp_inhibition_data, na.rm = TRUE, show_dose = TRUE)
basic_cyp_inhibition_risk_table <- function(perp, cyp_inh, na.rm = FALSE,
                                            show_dose = FALSE) {
  UseMethod("basic_cyp_inhibition_risk_table")
}


#' Basic CYP inhibition risk table
#' @inheritParams basic_cyp_inhibition_risk_table
#' @export
#' @noRd
basic_cyp_inhibition_risk_table.perpetrator <- function(
    perp, cyp_inh, na.rm = FALSE, show_dose = FALSE) {
  temp <- basic_cyp_inhibition_risk(perp, cyp_inh) %>%
    mutate(r = round(r, digits = 3)) %>%
    mutate(r_gut = round(r_gut, 1)) %>%
    mutate(risk_hep = case_match(as.character(risk_hep),
      "TRUE" ~ "Yes", "FALSE" ~ "No", .default = "")) %>%
    mutate(risk_intest = case_match(as.character(risk_intest),
      "TRUE" ~ "Yes", "FALSE" ~ "No", .default = "")) %>%
    mutate(r_gut = case_when(
      !is.na(r_gut) ~ as.character(r_gut), .default = ""))
    # mutate(risk_intest = case_when(
    #   !is.na(risk_intest) ~ as.character(risk_intest), .default = ""))

  if(na.rm==TRUE) {
    temp <- temp %>%
      filter(!is.na(ki))}
  labels <- c("CYP", "$K_{i}$ ($\\mu M$)", "$K_{i,u}$ ($\\mu M$)", "$R$",
              "risk (hepatic)", "$R_{gut}$", "risk (intestinal)")
  if(nrow(temp)!=0) {
    caption <- paste0("Risk for direct CYP inhibition by ", name(perp),
                      conditional_dose_string(perp, show_dose),
                      ", basic model")
    out <- knitr::kable(temp,
      caption = caption,
      col.names=labels)
    return(out)
  }
}


#' Basic CYP inhibition risk table
#' @inheritParams basic_cyp_inhibition_risk_table
#' @export
#' @noRd
basic_cyp_inhibition_risk_table.list <- function(perp, ...) {
  for(i in perp) {
    temp = basic_cyp_inhibition_risk_table.perpetrator(i, ...)
    if(!is.null(temp)){
      print(temp)
    }
  }
}


#' Basic modeling of the CYP time-dependent inhibition risk
#'
#' This function calculates the risk for time-dependent inhibition of CYP
#' enzymes.
#'
#' @details
#' The risk assessment is based on:
#'
#' \deqn{R=\frac {k_{obs} + k_{deg}}{k_{deg}}}
#'
#' where
#'
#' \deqn{k_{obs}=\frac {5*k_{inact}*C_{max,u}}{K_{I,u} + 5 * C_{max,u}}}
#'
#' Values of \eqn{R > 1.25} suggest a relevant TDI potential.
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
#' * 'ki' The \eqn{K_I} in \eqn{\mu M} as numeric.
#' * 'kinact' The \eqn{k_{inact}} in 1/h as numeric.
#' * 'source' Optional source information as character,
#'
#' @param cyp_kdeg The CYP turnover data as data frame. Defaults to the
#' built-in reference data, [cyp_turnover].
#' @seealso [cyp_turnover]
#' @seealso [read_tdi_data()]
#' @return A data frame.
#' @export
#' @examples
#' basic_cyp_tdi_risk(examplinib_parent, examplinib_cyp_tdi_data)
basic_cyp_tdi_risk<- function(perp, cyp_tdi, cyp_kdeg=cyp_turnover) {
  cyp_tdi <- cyp_tdi %>%
    filter(name==name(perp))

  i <- key_concentrations(perp, molar=TRUE)
  imaxssu <- i["imaxssu"]
  fumic <- as.num(perp["fumic", "value"])
  fu <- as.num(perp["fu", "value"])

  cyp_tdi %>%
    mutate(kobs=kinact*5*imaxssu/(ki * fu + 5 * imaxssu)) %>%
    mutate(fu=fu) %>%
    left_join(cyp_kdeg, by="cyp") %>%
    mutate(kdeg=kdeg_hepatic) %>%
    mutate(r=(kobs + kdeg)/kdeg) %>%
    mutate(risk=(r>1.25)) %>%
    select(cyp, ki, fu, kinact, kdeg, source, r, risk)
}


#' Basic CYP time-dependent inhibition risk table
#'
#' @inheritParams basic_cyp_tdi_risk
#' @param na.rm Switch to define whether rows with lacking \eqn{K_i} data are
#' removed from the output. Defaults to `FALSE`.
#' @param show_dose Show perpetrator dose in table title, defaults to `FALSE.`
#' @return A markdown-formatted table.
#' @seealso [basic_cyp_tdi_risk()]
#' @export
#' @examples
#' basic_cyp_tdi_risk_table(examplinib_parent, examplinib_cyp_tdi_data)
#' basic_cyp_tdi_risk_table(examplinib_compounds, examplinib_cyp_tdi_data)
#' basic_cyp_tdi_risk_table(examplinib_compounds, examplinib_cyp_tdi_data,
#'   show_dose = TRUE)
basic_cyp_tdi_risk_table <- function(perp, cyp_tdi, cyp_kdeg = cyp_turnover,
                                     na.rm = TRUE, show_dose = FALSE) {
  UseMethod("basic_cyp_tdi_risk_table")
}


#' Basic CYP time-dependent inhibition risk table
#' @inheritParams basic_cyp_tdi_risk_table
#' @export
#' @noRd
basic_cyp_tdi_risk_table.perpetrator <- function(
    perp, cyp_tdi, cyp_kdeg=cyp_turnover, na.rm = TRUE, show_dose = FALSE) {
  temp <- basic_cyp_tdi_risk(perp, cyp_tdi, cyp_kdeg) %>%
    mutate(risk = case_match(
      as.character(risk),
      "TRUE" ~ "Yes", "FALSE" ~ "No", .default = ""))

  if(na.rm==TRUE) {
    temp <- temp %>%
      filter(!is.na(ki) & !is.na(kdeg))
  }

  labels <- c("CYP", "$K_{I}$ ($\\mu M$)", "$f_u$", "$k_{inact}$ (1/h)",
              "$k_{deg}$ (1/h)", "source", "$R$", "risk")
  if(nrow(temp)!=0) {
    caption <- paste0("Risk for CYP TDI by ", name(perp),
                      conditional_dose_string(perp, show_dose),
                      ", basic model")
    out <- knitr::kable(temp,
                        caption = caption,
                        digits = 2,
                        col.names=labels)
    return(out)
  }
}


#' Basic CYP time-dependent inhibition risk table
#' @inheritParams basic_cyp_tdi_risk
#' @export
#' @noRd
basic_cyp_tdi_risk_table.list <- function(perp, ...) {
  for(i in perp) {
    temp <- basic_cyp_tdi_risk_table.perpetrator(i, ...)
    if(!is.null(temp)) {
      print(temp)
    }
  }
}



#### CYP INDUCTION

#' Basic static CYP induction risk
#'
#' @details
#' The basic fold-change method evaluates whether the maximal fold-change in
#' mRNA expression for the CYP enzyme, relative to the vehicle control, is >
#' 2-fold at concentrations up to 50-fold over the unbound Cmax,ss - provided
#' that concentration-dependent increases in mRNA were observed in vitro.
#'
#' If in the in vitro experiments concentrations of 50-fold the unbound Cmax,ss
#' were not tested, this is noted in the output.
#'
#' Also note that the ICH guideline further specifies that the positive control
#' used in the in vitro studies should result in at least 6-fold mRNA increases,
#' otherwise an induction risk cannot be ruled out if the mRNA increase for the
#' CYP enzyme is > 20% of the positive control. This, as well as the condition
#' of concentration-dependent increase in mRNA increase (see above), must be
#' confirmed manually by the user.
#'
#' For details, refer to section 2.1.4.1 of the
#' [ICH M1s2 guideline](https://www.ema.europa.eu/en/documents/scientific-guideline/ich-m12-guideline-drug-interaction-studies-step-5_en.pdf)
#'
#'
#' @param perp The perpetrator object or a list thereof.
#' @param cyp_ind The CYP induction data as data frame. The following fields
#' are expected:
#' * 'name' The name of the perpetrator compound as character.
#' * 'cyp' The CYP enzyme as (upper case) character.
#' * 'emax' The \eqn{E_{max}}, i.e., the maximum induction effect determined in
#' vitro as numeric.
#' * 'ec50' The \eqn{EC_{50}} in \eqn{\mu M} as numeric.
#' * 'maxc' The maximal concentration in \eqn{\mu M} tested in the in vitro assay as
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
      maxc_imaxssu<50 ~ "Low maxc",
      .default="")) %>%
    select(-c(name, ec50))
}


#' Table of the basic static CYP induction risk
#'
#' @inheritParams static_cyp_induction_risk
#' @param na.rm Switch to define whether rows with lacking \eqn{E_{max}} data
#' are removed from the output. Defaults to `FALSE`.
#' @param show_dose Show perpetrator dose in table title, defaults to `FALSE.`
#' @return A markdown-formatted table.
#' @export
#' @seealso [static_cyp_induction_risk()]
#' @examples
#' static_cyp_induction_risk_table(examplinib_parent,
#'   examplinib_cyp_induction_data)
#' static_cyp_induction_risk_table(examplinib_compounds,
#'   examplinib_cyp_induction_data)
#' static_cyp_induction_risk_table(examplinib_compounds,
#'   examplinib_cyp_induction_data, show_dose = TRUE)
static_cyp_induction_risk_table <- function(
    perp, cyp_ind, na.rm = F, show_dose = FALSE) {
  UseMethod("static_cyp_induction_risk_table")
}


#' Table of the basic static CYP induction risk
#' @inheritParams static_cyp_induction_risk_table
#' @export
#' @noRd
static_cyp_induction_risk_table.perpetrator <- function(
    perp, cyp_ind, na.rm = F, show_dose = FALSE) {
  temp <- static_cyp_induction_risk(perp, cyp_ind) %>%
    mutate(risk = case_match(
      as.character(risk),
      "TRUE" ~ "Yes", "FALSE" ~ "No", .default = ""))

  if(na.rm==TRUE) {
    temp <- temp %>%
      filter(!is.na(emax))}
  labels <- c("CYP", "$E_{max}$", "$max c$ ($\\mu M$)", "source",
              "$max c/C_{max,ss,u}$", "risk", "notes")
  if(nrow(temp)!=0) {
    caption <- paste0("Risk for hepatic CYP induction by ", name(perp),
                      conditional_dose_string(perp, show_dose),
                      ", basic static model")
    out <- knitr::kable(
      temp, caption = caption,
      col.names=labels)
    return(out)
  }
}


#' Table of the basic static CYP induction risk
#' @inheritParams static_cyp_induction_risk_table
#' @export
#' @noRd
static_cyp_induction_risk_table.list <- function(perp, ...) {
  for(i in perp) {
    temp <- static_cyp_induction_risk_table.perpetrator(i, ...)
    if(!is.null(temp)) {
      print(temp)
    }
  }
}


#' Basic kinetic CYP induction risk
#'
#' @details
#' Basic kinetic modeling of the CYP induction risk considers:
#'
#' \deqn{R = \frac {1}{1 + d * \frac {E_{max}*10*C_{max,ss,u}}{EC_{50,u} + 10 * C_{max,ss,u}}}}
#'
#' \eqn{d} is a scaling factor with a standard value of 1. A different value can
#' be used if warranted by prior experience with the experimental conditions.
#'
#' \eqn{R \le 0.8} suggest a relevant in vivo CYP induction potential.
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
    # mutate(d = d) %>%
    mutate(r=(1/(1 + d * emax * 10 * i["imaxssu"] /
                    (ec50 + 10 * i["imaxssu"])))) %>%
    mutate(risk = r <= 0.8) %>%
    select(-name)
  return(out)
}


#' Table of the basic kinetic CYP induction risk
#'
#' @inheritParams kinetic_cyp_induction_risk
#' @param na.rm Switch to define whether rows with lacking \eqn{E_{max}} data
#' are removed from the output. Defaults to `FALSE`.
#' @param show_dose Show perpetrator dose in table title, defaults to `FALSE.`
#' @return A markdown-formatted table.
#' @export
#' @seealso [kinetic_cyp_induction_risk()]
#' @examples
#' kinetic_cyp_induction_risk_table(examplinib_parent,
#'   examplinib_cyp_induction_data)
#' kinetic_cyp_induction_risk_table(examplinib_compounds,
#'   examplinib_cyp_induction_data)
#' kinetic_cyp_induction_risk_table(examplinib_compounds,
#'   examplinib_cyp_induction_data, show_dose = TRUE)
kinetic_cyp_induction_risk_table <- function(
    perp, cyp_ind, na.rm = F, show_dose = FALSE) {
  UseMethod("kinetic_cyp_induction_risk_table")
}


#' Table of the basic kinetic CYP induction risk
#' @inheritParams kinetic_cyp_induction_risk_table
#' @export
#' @noRd
kinetic_cyp_induction_risk_table.perpetrator <- function(
    perp, cyp_ind, na.rm = F, show_dose = FALSE) {
  temp <- kinetic_cyp_induction_risk(perp, cyp_ind) %>%
    mutate(risk = case_match(
      as.character(risk),
      "TRUE" ~ "Yes", "FALSE" ~ "No", .default = "")) %>%
    mutate(r=format(r, digit=2)) %>%
    select(-maxc)

  if(na.rm==TRUE) {
    temp <- temp %>%
      filter(!is.na(emax))
  }

  labels <- c("CYP", "$E_{max}$", "$EC_{50}$ ($\\mu M$)",
              "source", "$R$", "risk")
  if(nrow(temp)!=0) {
    caption <- paste0("Risk for CYP induction by ", name(perp),
          conditional_dose_string(perp, show_dose), ", basic kinetic model")
    out <- knitr::kable(temp, caption = caption, col.names = labels)
    return(out)
  }
}


#' Table of the basic kinetic CYP induction risk
#' @inheritParams kinetic_cyp_induction_risk_table
#' @export
#' @noRd
kinetic_cyp_induction_risk_table.list <- function(perp, ...) {
  for(i in perp) {
    temp <- kinetic_cyp_induction_risk_table.perpetrator(i, ...)
    if(!is.null(temp)) {
      print(temp)
    }
  }
}



#### MECHANISTIC STATIC MODELING

#' CYP perpetration risk as per mechanistic-static modeling
#'
#' @details
#' In this approach, AUC ratios for CYP-specific probe substrates are calculated
#' based on their known intestinal and hepatic metabolism. Direct (competitive)
#' and time-dependent inhibition terms, as well as enzyme induction terms are
#' considered. For details, refer to section 7.5.1.2 of the [ICH M12
#' guideline](https://www.ema.europa.eu/en/documents/scientific-guideline/ich-m12-guideline-drug-interaction-studies-step-5_en.pdf).
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
#' @param show_dose Show perpetrator dose in table title, defaults to `FALSE.`
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
#'   examplinib_compounds,
#'   examplinib_cyp_inhibition_data,
#'   examplinib_cyp_induction_data,
#'   examplinib_cyp_tdi_data)
#' mech_stat_cyp_risk_table(
#'   examplinib_compounds,
#'   examplinib_cyp_inhibition_data,
#'   examplinib_cyp_induction_data,
#'   examplinib_cyp_tdi_data, show_dose = TRUE)
mech_stat_cyp_risk_table <- function(
    perp,
    cyp_inh,
    cyp_ind,
    cyp_tdi = NULL,
    include_induction = TRUE,
    substr = cyp_reference_substrates,
    cyp_kdeg = cyp_turnover,
    na.rm = FALSE,
    show_dose = FALSE) {
  UseMethod("mech_stat_cyp_risk_table")
}


#' CYP perpetration risk table as per mechanistic-static modeling
#' @inheritParams mech_stat_cyp_risk_table
#' @export
#' @noRd
mech_stat_cyp_risk_table.perpetrator <- function(
    perp,
    cyp_inh,
    cyp_ind,
    cyp_tdi = NULL,
    include_induction = TRUE,
    substr = cyp_reference_substrates,
    cyp_kdeg = cyp_turnover,
    na.rm = FALSE,
    show_dose = FALSE) {
  temp <- mech_stat_cyp_risk(perp, cyp_inh, cyp_ind, cyp_tdi,
                             include_induction=include_induction,
                             substr=substr, cyp_kdeg=cyp_kdeg) %>%
    mutate(risk = case_match(
      as.character(risk),
      "TRUE" ~ "Yes", "FALSE" ~ "No", .default = "")) %>%
    select(cyp, substrate, fgut, fm, fmcyp, Ag, Ah, Bg, Bh, Cg, Ch, aucr,
           risk) %>%
    # mutate(across(Ag:Ch, ~ signif(., digits=2))) %>%
    # mutate(across(Ag:Ch, ~ format(., digits=2, nsmall=3))) %>%
    mutate(across(Ag:Ch, ~ round(., digits=2))) %>%
    mutate(aucr=format(aucr, digits=3))

  if(na.rm==TRUE) {
    temp <- temp %>%
      filter(!is.na(risk))
  }

  labels <- c("CYP", "substrate", "$F_{gut}$", "$f_m$",
              "$f_{m,CYP}$", "$A_g$", "$A_h$", "$B_g$", "$B_h$", "$C_g$",
              "$C_h$", "AUCR", "risk")

  if(nrow(temp)!=0) {
    caption <- paste0("Mechanistic static modeling of the CYP inhibition risk for ",
      name(perp), conditional_dose_string(perp, show_dose))
    out <- knitr::kable(
      temp,
      caption = caption,
      col.names = labels)
    return(out)
  }
}


#' CYP perpetration risk table as per mechanistic-static modeling
#' @inheritParams mech_stat_cyp_risk_table
#' @export
#' @noRd
mech_stat_cyp_risk_table.list <- function(perp, ...) {
  for(i in perp) {
    temp <- mech_stat_cyp_risk_table.perpetrator(i, ...)
    if(!is.null(temp)) {
      print(temp)
    }
  }
}











