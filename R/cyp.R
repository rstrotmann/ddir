#### BASIC MODELING


#' Generic function to evaluate the basic CYP inhibition risk
#'
#' @param perp The perpetrator as perpetrator object or list thereof.
#' @param cyp_inh The CYP inhibitor data as data frame.
#'
#' @return A data frame.
#' @export
basic_cyp_inhibition_risk <- function(perp, cyp_inh) {
  UseMethod("basic_cyp_inhibition_risk")
}


#' Basic CYP inhibition risk evaluation for a perpetrator object
#'
#' @param perp The perpetrator object.
#' @param cyp_inh CYP inhibiton data as data frame.
#'
#' @return A data frame.
#' @export
#' @examples
#' basic_cyp_inhibition_risk(examplinib_parent, examplinib_cyp_inhibition_data)
#'
basic_cyp_inhibition_risk.perpetrator <- function(perp, cyp_inh) {
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


#' Basic CYP inhibition risk evaluation for a list of perpetrator objects
#'
#' @param perp Perpetrator objects as list.
#' @param cyp_inh CYP inhibition data as data frame.
#'
#' @return A list of data frames.
#' @export
#' @examples
#' basic_cyp_inhibition_risk(examplinib_compounds, examplinib_cyp_inhibition_data)

basic_cyp_inhibition_risk.list <- function(perp, cyp_inh) {
  lapply(perp, basic_cyp_inhibition_risk, cyp_inh=cyp_inh)
}


#' Generic function for basic CYP inhibition risk evaluation as table
#'
#' @param perp The perpetrator as perpetrator object or list thereof.
#' @param cyp_inh CYP inhibition data as data frame.
#' @param na.rm Remove rows with lacking ki data (i.e., where ki == NA).
#' @export
#' @examples
#' basic_cyp_inhibition_risk_table(examplinib_parent, examplinib_cyp_inhibition_data)
basic_cyp_inhibition_risk_table <- function(perp, cyp_inh, na.rm=F) {
  UseMethod("basic_cyp_inhibition_risk_table")
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
#' @examples
#' basic_cyp_inhibition_risk_table(examplinib_parent, examplinib_cyp_inhibition_data)
basic_cyp_inhibition_risk_table.perpetrator <- function(perp, cyp_inh, na.rm=F) {
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


#' Table of basic CYP inhibition risk
#'
#' @param perp A list of perpetrator objects.
#' @param cyp_inh CYP inhibition data as data frame.
#' @param na.rm Remove rows with lacking ki data (i.e., where ki == NA).
#'
#' @return Basic evaluation of CYP inhibition risk as markdown-formatted tables,
#'   or an empty string if no CYP inhibition data available.
#' @export
#' @examples
#' basic_cyp_inhibition_risk_table(examplinib_compounds, examplinib_cyp_inhibition_data)
basic_cyp_inhibition_risk_table.list <- function(perp, cyp_inh, na.rm=F) {
  for(i in perp) {
    print(basic_cyp_inhibition_risk_table(i, cyp_inh, na.rm=na.rm))
  }
}



#### CYP INDUCTION

#' Generic function for basic static CYP induction risk
#'
#' @param perp The perpetrator object or a list thereof.
#' @param cyp_ind The CYP induction data as data frame.
#'
#' @return A data frame.
#' @export
#' @examples
#' static_cyp_induction_risk(examplinib_parent, examplinib_cyp_induction_data)
static_cyp_induction_risk <- function(perp, cyp_ind) {
  UseMethod("static_cyp_induction_risk")
}


#' Basic static CYP induction risk
#'
#' @param perp The perpetrator object or a list thereof.
#' @param cyp_ind The CYP induction data as data frame.
#'
#' @return A data frame.
#' @export
#' @examples
#' static_cyp_induction_risk(examplinib_parent, examplinib_cyp_induction_data)
static_cyp_induction_risk.perpetrator <- function(perp, cyp_ind)  {
  i <- key_concentrations(perp, molar=TRUE)

  cyp_ind %>%
    filter(name==name(perp)) %>%
    mutate(maxc_imaxssu=round(maxc/i["imaxssu"], 1)) %>%
    mutate(risk=emax>2) %>%
    mutate(note=case_when(
      maxc_imaxssu<30~"Maximal tested concentration is below EMA/FDA expectations",
      (maxc_imaxssu>30 & maxc_imaxssu<50)~"Maximal tested concentration is below FDA expectations",
      .default="")) %>%
    select(-c(name, ec50))
}


#' Basic static CYP induction risk
#'
#' @param perp The perpetrator object or a list thereof.
#' @param cyp_ind The CYP induction data as data frame.
#'
#' @return A data frame.
#' @export
#' @examples
#' static_cyp_induction_risk(examplinib_compounds, examplinib_cyp_induction_data)
static_cyp_induction_risk.list <- function(perp, cyp_ind) {
  lapply(perp, static_cyp_induction_risk, cyp_ind=cyp_ind)
}


#' Generic function: Table of the basic static CYP induction risk
#'
#' @param perp The perpetrator object or list thereof
#' @param cyp_ind The CYP induction data as data frame.
#' @param na.rm Remove rows with lacking ki data (i.e., where ki == NA).
#' @export
#' @examples
#' static_cyp_induction_risk_table(examplinib_parent, examplinib_cyp_induction_data)
static_cyp_induction_risk_table <- function(perp, cyp_ind, na.rm=F) {
  UseMethod("static_cyp_induction_risk_table")
}


#' Table of the basic static CYP induction risk
#'
#' @param perp The perpetrator object.
#' @param cyp_ind The CYP induction data as data frame.
#' @param na.rm Remove rows with lacking ki data (i.e., where ki == NA).
#'
#' @return A markdown-formatted table.
#' @export
#' @examples
#' static_cyp_induction_risk_table(examplinib_parent, examplinib_cyp_induction_data)
static_cyp_induction_risk_table.perpetrator <- function(perp, cyp_ind, na.rm=F) {
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


#' Tables of the basic static CYP induction risk
#'
#' @param perp The list of perpetrator objects.
#' @param cyp_ind The CYP induction data as data frame.
#' @param na.rm Remove rows with lacking ki data (i.e., where ki == NA).
#'
#' @return A list of markdown-formatted tables.
#' @export
#' @examples
#' static_cyp_induction_risk_table(examplinib_compounds, examplinib_cyp_induction_data)
static_cyp_induction_risk_table.list <- function(perp, cyp_ind, na.rm=F) {
  for(i in perp) {
    print(static_cyp_induction_risk_table(i, cyp_ind, na.rm=na.rm))
  }
}





#' Generic function for basic kinetic CYP induction risk
#'
#' @param perp The perpetrator object or a list thereof.
#' @param cyp_ind The CYP induction data as data frame.
#' @export
#' @examples
#' kinetic_cyp_induction_risk(examplinib_parent, examplinib_cyp_induction_data)
kinetic_cyp_induction_risk <- function(perp, cyp_ind) {
  UseMethod("kinetic_cyp_induction_risk")
}


#' Basic kinetic CYP induction risk
#'
#' @param perp The perpetrator object.
#' @param cyp_ind The CYP induction data as data frame.
#'
#' @return A data frame.
#' @export
#' @examples
#' kinetic_cyp_induction_risk(examplinib_parent, examplinib_cyp_induction_data)
kinetic_cyp_induction_risk.perpetrator <- function(perp, cyp_ind) {
  i <- key_concentrations(perp, molar=TRUE)

  out <- cyp_ind %>%
    filter(name==name(perp)) %>%
    mutate(r3=(1/(1 + emax*10*i["imaxssu"]/(ec50 + 10*i["imaxssu"])))) %>%
    mutate(risk=r3<=0.8) %>%
    select(-name)
  return(out)
}


#' Basic kinetic CYP induction risk
#'
#' @param perp A list of perpetrator objects.
#' @param cyp_ind The CYP induction data as data frame.
#'
#' @return A list of data frames.
#' @export
#' @examples
#' kinetic_cyp_induction_risk(examplinib_compounds, examplinib_cyp_induction_data)
kinetic_cyp_induction_risk.list <- function(perp, cyp_ind) {
  lapply(perp, kinetic_cyp_induction_risk, cyp_ind=cyp_ind)
}


#' Generic function: Table of the basic kinetc CYP induction risk
#'
#' @param perp The perpetrator object or list thereof
#' @param cyp_ind The CYP induction data as data frame.
#' @param na.rm Remove rows with lacking ki data (i.e., where ki == NA).
#' @export
#' @examples
#' kinetic_cyp_induction_risk_table(examplinib_parent, examplinib_cyp_induction_data)
kinetic_cyp_induction_risk_table <- function(perp, cyp_ind, na.rm=F) {
  UseMethod("kinetic_cyp_induction_risk_table")
}


#' Table of the basic kinetc CYP induction risk
#'
#' @param perp The perpetrator object.
#' @param cyp_ind The CYP induction data as data frame.
#' @param na.rm Remove rows with lacking ki data (i.e., where ki == NA).
#'
#' @return A markdown-formatted table.
#' @export
#' @examples
#' kinetic_cyp_induction_risk_table(examplinib_parent, examplinib_cyp_induction_data)
kinetic_cyp_induction_risk_table.perpetrator <- function(perp, cyp_ind, na.rm=F) {
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


#' Table of the basic kinetc CYP induction risk
#'
#' @param perp A list of perpetrator objects.
#' @param cyp_ind The CYP induction data as data frame.
#' @param na.rm Remove rows with lacking ki data (i.e., where ki == NA).
#'
#' @return A markdown-formatted table.
#' @export
kinetic_cyp_induction_risk_table.list <- function(perp, cyp_ind, na.rm=F) {
  for(i in perp) {
    print(kinetic_cyp_induction_risk_table(i, cyp_ind, na.rm=na.rm))
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
#' @export
#' @examples
#' mech_stat_cyp_risk_table(examplinib_parent, examplinib_cyp_inhibition_data, examplinib_cyp_induction_data)
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













