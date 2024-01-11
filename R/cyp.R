#### BASIC MODELING

#' Basic CYP inhibition risk
#'
#' This function evaluates the clinical risk for direct (reversible) CYP
#' inhibition according to the basic model defined in the relevant regulatory
#' guideines.
#'
#' @details
#' For the basic modeling of direct (reversible) CYP enzyme inhibition, the
#' ratios of the relevant inhibitor concentration to the \eqn{K_i} are
#' considered, i.e., \eqn{R_1} for hepatic enzymes and \eqn{R_{1,gut}} for
#' intestinal enzymes (refer to fig. 1 of  the the
#' [FDA guidance](https://www.fda.gov/media/134582/download)). A threshold
#' of 1.02 applies.
#'
#' Note that section 5.3.3.1 of the
#' [EMA DDI guidance](https://www.ema.europa.eu/en/documents/scientific-guideline/guideline-investigation-drug-interactions-revision-1_en.pdf)
#' defines the ratio as \eqn{R_1=\frac{I_{max,ss,u}}{K_{i,u}}} and the
#' threshold as 0.02.
#'
#' ## Liver
#'
#' \deqn{R_1=1+\frac{I_{max,ss,u}}{K_{i,u}}}
#'
#' ## Gut wall
#'
#' \deqn{R_{1,gut}=1+\frac{I_{gut}}{K_{i,u}}}
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
#' @export
#' @examples
#' basic_cyp_inhibition_risk(examplinib_parent, examplinib_cyp_inhibition_data)
#'
basic_cyp_inhibition_risk <- function(perp, cyp_inh) {
  ki <- cyp_inh %>%
    filter(name==name(perp)) %>%
    filter(param!="name")

  i <- key_concentrations(perp, molar=TRUE)
  fumic <- as.num(perp["fumic", "value"])

  out <- ki %>%
    mutate(ki=as.num(value)) %>%
    mutate(kiu=as.num(value)*fumic) %>%
    mutate(r1=1 + (i["imaxssu"]/kiu)) %>%
    mutate(r1gut=case_when(param=="CYP3A4" ~ 1+(i["igut"]/ki), .default=NA)) %>%
    select(-name, -value) %>%
    mutate(risk_hep=r1>1.02) %>%
    mutate(risk_intest=r1gut>11) %>%
    mutate(r1=format(r1, digits=4)) %>%
    mutate(r1gut=format(r1gut, digits=4)) %>%
    select(param, ki, kiu, r1, risk_hep, r1gut, risk_intest)

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
#' @param na.rm Boolean to define whether rows with lacking ki data are removed
#' (i.e., where ki == NA).
#' @return A markdown-formatted table, or an empty string (if no CYP inhibition
#' data are vavailable).
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
#' Basic kinetic modeling of the CYP induction risk considers \eqn{R_3}:
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
#' The below formula given by the FDA guideline (refer to FDA 2020, Fig. 7) also
#' includes intestinal and hepatic enzyme induction terms (\eqn{C_g} and
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
#' @param cyp_inh CYP inhibiton data as data frame.
#' @param cyp_ind CYP induction data as data frame.
#' @param substr The CYP reference substrates to be used as data frame.
#' @param include_induction Boolean value to define whether induction should be
#'   included in the calculation (C-terms as per the FDA guideline)
#' @return A data frame.
#' @export
#' @examples
#' mech_stat_cyp_risk(examplinib_parent, examplinib_cyp_inhibition_data, examplinib_cyp_induction_data)
#'
mech_stat_cyp_risk <- function(
    perp,
    cyp_inh,
    cyp_ind,
    include_induction=T,
    substr=cyp_reference_substrates) {
  fumic <- as.num(perp["fumic", "value"])

  i <- key_concentrations(perp, molar=TRUE)
  Ig <- i["imaxintestu"]
  Ih <- i["imaxinletu"]

  out <- cyp_inh %>%
    filter(name==name(perp)) %>%
    filter(param!="name") %>%
    mutate(cyp=param) %>%
    mutate(ki=as.num(value)) %>%
    left_join(cyp_ind %>% filter(name==name(perp)), by=c("cyp", "name")) %>%
    mutate(kiu=ki*fumic) %>%
    mutate(Ag=1/(1+(Ig/kiu))) %>%
    mutate(Ah=1/(1+(Ih/kiu))) %>%
    mutate(Cg=1, Ch=1) %>% # temporary
    mutate(Cg=1+(emax*Ig/(Ig+ec50))) %>%
    mutate(Cg=case_when((is.na(Cg) | include_induction==F)~1, .default=Cg)) %>%
    mutate(Ch=1+(emax*Ih/(Ih+ec50))) %>%
    mutate(Ch=case_when((is.na(Ch) | include_induction==F)~1, .default=Ch))  %>%
    select(-name, -value, -starts_with("source"), -param) %>%
    left_join(substr, by="cyp") %>%
    mutate(aucr=1/(Ag*Cg*(1-fgut)+fgut) * 1/(Ah*Ch*fm*fmcyp+(1-fm*fmcyp))) %>%
    mutate(risk=aucr>1.25 | aucr < 0.8) %>%
    select(cyp, substrate, kiu, fgut, fm, fmcyp, Ag, Ah, Cg, Ch, aucr, risk) %>%
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
#'
#' @return A markdown-formatted table.
#' @seealso [mech_stat_cyp_risk()]
#' @export
#' @examples
#' mech_stat_cyp_risk_table(examplinib_parent, examplinib_cyp_inhibition_data, examplinib_cyp_induction_data)
#'
mech_stat_cyp_risk_table <- function(
    perp,
    cyp_inh,
    cyp_ind,
    include_induction=T,
    substr=cyp_reference_substrates,
    na.rm=FALSE) {
  temp <- mech_stat_cyp_risk(perp, cyp_inh, cyp_ind,
                             include_induction=include_induction,
                             substr=substr) %>%
    select(cyp, kiu, substrate, fgut, fm, fmcyp, Ag, Ah, Cg, Ch, aucr, risk) %>%
    mutate(across(Ag:Ch, ~ format(.x, digits=3))) %>%
    mutate(aucr=format(aucr, digits=3))

  if(na.rm==TRUE) {
    temp <- temp %>%
      filter(!is.na(kiu))
  }

  labels <- c("CYP", "$K_{i,u}$", "substrate", "$F_{gut}$", "$f_m$",
              "$f_{m,CYP}$", "$A_g$", "$A_h$", "$C_g$", "$C_h$", "AUCR", "risk")

  if(nrow(temp)!=0) {
    out <- knitr::kable(
      temp,
      caption=paste("Mechanistic static modeling of CYP inhibition risk for ",
                    name(perp)),
      col.names=labels)
    return(out)
  }
}













