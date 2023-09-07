#' Load DMPK data (CYP inhibition, transporter inhibition, etc.)
#'
#' @param filename The filename.
#'
#' @return A list of data frames.
#' @import dplyr
#' @export
load_dmpk_data <- function(filename) {
  raw <- as.data.frame(read.csv(filename,
                                col.names=c("name", "param", "value", "source"),
                                header = F,
                                comment.char = '#')) %>%
    dplyr::mutate(across(everything(), trimws)) %>%
    dplyr::group_by(name) %>%
    dplyr::group_modify(~ tibble::add_row(param="name", value=.y$name, source="", .x, , .before=1)) %>%
    dplyr::ungroup() %>%
    as.data.frame()
  return(raw)
}

#' Load CYP inhibitor data
#'
#' @param filename The filename as string.
#' @return The data as data frame.
#' @export
load_cyp_inhibitor_data <- function(filename) {
  return(load_dmpk_data(filename))
}



#' Read csv-formatted DMPK data from string
#'
#' @param x The string containing the data in csv format.
#' @return The data as data frame.
#' @export
read_string_dmpk_data <- function(x) {
  raw <- as.data.frame(read.csv(text=x,
                                col.names=c("name", "param", "value", "source"),
                                header = F,
                                comment.char = '#')) %>%
    dplyr::mutate(across(everything(), trimws)) %>%
    dplyr::group_by(name) %>%
    dplyr::group_modify(~ tibble::add_row(param="name", value=.y$name, source="", .x, , .before=1)) %>%
    dplyr::ungroup() %>%
    as.data.frame()
  return(raw)
}


#' Load CYP inducer data from file
#'
#' @param filename The filename.
#' @return A list of data frames.
#' @export
load_cyp_inducer_data <- function(filename) {
  raw <- as.data.frame(read.csv(filename,
                                col.names=c("name", "cyp", "emax", "ec50",
                                            "maxc", "source"),
                                header = F,
                                comment.char = '#')) %>%
    dplyr::mutate(across(everything(), trimws)) %>%
    mutate(across(3:5, as.num)) %>%
    as.data.frame()
  return(raw)
}

#' Read csv-formatted CYP inducer data from string
#'
#' @param x The string containing the data in csv format.
#' @return The data as data frame.
#' @export
read_string_cyp_inducer_data <- function(x) {
  raw <- as.data.frame(read.csv(text=x,
                                col.names=c("name", "cyp", "emax", "ec50",
                                            "maxc", "source"),
                                header = F,
                                comment.char = '#')) %>%
    dplyr::mutate(across(everything(), trimws)) %>%
    mutate(across(3:5, as.num)) %>%
    as.data.frame()
  return(raw)
}

#### BASIC MODELING


#' Basic CYP inhibition parameters
#'
#' @param perp The perpetrator object.
#' @param cyp_inh CYP inhibiton data as data frame.
#'
#' @return Data frame.
#' @export
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


#' Table of basic CYP inhibition risk
#'
#' @param perp The perpetrator object.
#' @param cyp_inh CYP inhibition data as data frame.
#' @param na.rm Remove rows with lacking ki data (i.e., where ki == NA).
#'
#' @return Basic evaluation of CYP inhibition risk as markdown-formatted table,
#'   or an empty string if no CYP inhibition data available.
#' @export
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

#' Basic kinetc CYP induction risk
#'
#' @param perp The perpetrator object.
#' @param cyp_ind The CYP induction data as data frame.
#'
#' @return A data frame.
#' @export
kinetic_cyp_induction_risk <- function(perp, cyp_ind) {
  i <- key_concentrations(perp, molar=TRUE)

  out <- cyp_ind %>%
    filter(name==name(perp)) %>%
    mutate(r3=(1/(1 + emax*10*i["imaxssu"]/(ec50 + 10*i["imaxssu"])))) %>%
    mutate(risk=r3<=0.8) %>%
    select(-name)
}


#' Table of the basic kinetc CYP induction risk
#'
#' @param perp The perpetrator object.
#' @param cyp_ind The CYP induction data as data frame.
#' @param na.rm Remove rows with lacking ki data (i.e., where ki == NA).
#'
#' @return A markdown-formatted table.
#' @export
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
                        caption=paste("Risk for CYP inductio by", name(perp),
                                      "(basic kinetic model)"),
                        col.names=labels)
    return(out)
  }
}


#### MECHANISTIC STATIC MODELING

# cyp_reference_substrates <- data.frame(
#   cyp=c("CYP1A2", "CYP2C8", "CYP2C9", "CYP2C19", "CYP3A4"),
#   substrate=c("tizanidine", "repaglinide", "S-warfarin", "omeprazole", "midazolam"),
#   fgut=c(1, 1, 1, 1, 0.57),
#   fm=c(0.95, 1, 1, 1, 0.96),
#   fmcyp=c(0.98, 0.61, 0.91, 0.87, 1)
# )


#' CYP perpetration risk as per mechanistic-static modeling
#'
#' @param perp The perpetrator object.
#' @param cyp_inh CYP inhibiton data as data frame.
#' @param cyp_ind CYP induction data as data frame.
#' @param substr The CYP reference substrates to be used as data frame.
#' @param include_induction Boolean value to define whether induction should be
#'   included in the calculation (C-terms as per the FDA guideline)
#' @return A data frame.
#' @export
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
#' @export
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













