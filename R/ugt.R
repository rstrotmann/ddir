

#' UGT inhibition risk
#'
#' This function assumes that the UGT inhibition data is provided as ICC50.
#'   According to Cheng-Prusoff, Ki can be assumed to be IC50/2 at the
#'   usual experimental conditions where substrate concentrations close to Km
#'   are used.
#'
#' @param perp The perpetrator object.
#' @param ugt_inh UGT inhibiton data as data frame, with value repersenting the
#'   IC50 (!).
#'
#' @return A markdown-formatted table.
#' @export
basic_ugt_inhibition_risk <- function(perp, ugt_inh) {
  ki <- ugt_inh %>%
    filter(name==name(perp)) %>%
    filter(param!="name")

  i <- key_concentrations(perp, molar=TRUE)
  fumic <- as.num(perp["fumic", "value"])

  out <- ki %>%
    mutate(ki=as.num(value)/2) %>%
    mutate(kiu=ki*fumic) %>%
    mutate(r1=1 + (i["imaxssu"]/kiu)) %>%
    select(-name, -value) %>%
    mutate(risk=r1>1.02) %>%
    select(param, kiu, r1, risk)
  return(out)
}


#' UGT inhibition risk table
#'
#' @param perp The perpetrator object.
#' @param ugt_inh UGT inhibiton data as data frame, with value repersenting the
#'   IC50 (!).
#' @param na.rm Remove rows with lacking ki data (i.e., where ki == NA).
#'
#' @return A markdown-formatted table.
#' @export
basic_ugt_inhibition_risk_table <- function(perp, ugt_inh, na.rm=F) {
  temp <- basic_ugt_inhibition_risk(perp, ugt_inh)

  if(na.rm==TRUE) {
    temp <- temp %>%
      filter(!is.na(ki))
  }

  labels <- c("UGT", "$K_{i,u}$", "$R_1$", "risk")
  if(nrow(temp)!=0) {
    out <- knitr::kable(temp,
                        caption=paste("Risk for UGT inhibition by", name(perp), "(basic model)"),
                        col.names=labels)
    return(out)
  }
}
