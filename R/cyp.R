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
#' @param ra.rm Remove rows with lacking ki data (i.e., where ki == NA).
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








