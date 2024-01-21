

#' UGT inhibition risk
#'
#' This function evaluates the clinical risk for reversible inhibition of UGT
#' enzymes according to the relevant
#' [FDA](https://www.fda.gov/media/134582/download) and
#' [EMA](https://www.ema.europa.eu/en/documents/scientific-guideline/guideline-investigation-drug-interactions-revision-1_en.pdf)
#' guidelines.
#'
#' @details
#' This function assumes that the UGT inhibition data is provided as
#' \eqn{IC_{50}}. According to
#' [Cheng, Prusoff 1973](https://doi.org/10.1016/0006-2952(73)90196-2)),
#' \eqn{K_i} can be assumed to be \eqn{IC_{50}/2} at the experimental conditions
#' commonly used in the in vitro ihibition studies where substrate
#' concentrations are close to \eqn{K_M}.
#'
#' @details
#' The relevant metric for basic modeling of the UGT inhibition risk is
#' \eqn{R_1=I_{max,ss,u}/K_{i,u}} (refer to fig. 1 of  the the
#' [FDA guidance](https://www.fda.gov/media/134582/download)).
#'
#' \eqn{R_1>1.02} are considered to indicate a potential UGT inhibition risk.
#'
#' Refer to the documentation to the [key_concentrations()] function for details
#' on the calculation of \eqn{I_{max,ss,u}}.
#'
#' @param perp The perpetrator object.
#' @param ugt_inh UGT inhibition data as data frame, The following fields are
#' expected:
#' * 'name' The name of the perpetrator compound.
#' * 'ugt' The UGT enzyme as (upper case) character.
#' * 'ic50' The \eqn{IC_{50}} in \eqn{\mu M} as numeric.
#' * 'source' Optional source information as character.
#' @seealso [key_concentrations()]
#' @seealso [basic_ugt_inhibition_risk_table()]
#' @seealso [read_ugt_inhibitor_data()]
#' @return A markdown-formatted table.
#' @export
#' @examples
#' basic_ugt_inhibition_risk(examplinib_parent, examplinib_ugt_inhibition_data)
basic_ugt_inhibition_risk <- function(perp, ugt_inh) {
  temp <- ugt_inh %>%
    filter(name==name(perp))

  i <- key_concentrations(perp, molar=TRUE)
  fumic <- as.num(perp["fumic", "value"])

  out <- temp %>%
    mutate(ki=as.num(ic50)/2) %>%
    mutate(kiu=ki*fumic) %>%
    mutate(r1=1 + (i["imaxssu"]/kiu)) %>%
    select(-name) %>%
    mutate(risk=r1>1.02) %>%
    select(ugt, kiu, r1, risk)
  return(out)
}


#' UGT inhibition risk table
#'
#' This function generates a markdown-formatted table of the reversible
#' UGT inhibition risk assessment. See [basic_ugt_inhibition_risk()] for details
#' on the calculation of the risk.
#'
#' @inheritParams basic_ugt_inhibition_risk
#' @param na.rm Boolean to define whether rows with lacking \eqn{K_i} data are
#' removed from the output (i.e., where `ki == NA`).
#' @seealso [basic_ugt_inhibition_risk()]
#' @return A markdown-formatted table.
#' @export
#' @examples
#' basic_ugt_inhibition_risk_table(examplinib_parent, examplinib_ugt_inhibition_data)
#' basic_ugt_inhibition_risk_table(examplinib_compounds, examplinib_ugt_inhibition_data)
basic_ugt_inhibition_risk_table <- function(perp, ugt_inh, na.rm=F) {
  UseMethod("basic_ugt_inhibition_risk_table")
}


#' UGT inhibition risk table
#'
#' This function generates a markdown-formatted table of the reversible
#' UGT inhibition risk assessment. See [basic_ugt_inhibition_risk()] for details
#' on the calculation of the risk.
#'
#' @inheritParams basic_ugt_inhibition_risk
#' @param na.rm Boolean to define whether rows with lacking \eqn{K_i} data are
#' removed from the output (i.e., where `ki == NA`).
#' @seealso [basic_ugt_inhibition_risk()]
#' @return A markdown-formatted table.
#' @export
#' @noRd
#' @examples
#' basic_ugt_inhibition_risk_table(examplinib_parent, examplinib_ugt_inhibition_data)
basic_ugt_inhibition_risk_table.perpetrator <- function(perp, ugt_inh, na.rm=F) {
  temp <- basic_ugt_inhibition_risk(perp, ugt_inh)

  if(na.rm==TRUE) {
    temp <- temp %>%
      filter(!is.na(ki))
  }

  labels <- c("UGT", "$K_{i,u}$", "$R_1$", "risk")
  if(nrow(temp)!=0) {
    out <- temp %>%
      mutate(r1=round(r1, 3)) %>%
      knitr::kable(caption=paste("Risk for UGT inhibition by",
                                 name(perp), "(basic model)"),
                   col.names=labels)
    return(out)
  }
}


#' UGT inhibition risk table
#'
#' This function generates a markdown-formatted table of the reversible
#' UGT inhibition risk assessment. See [basic_ugt_inhibition_risk()] for details
#' on the calculation of the risk.
#'
#' @inheritParams basic_ugt_inhibition_risk
#' @param na.rm Boolean to define whether rows with lacking \eqn{K_i} data are
#' removed from the output (i.e., where `ki == NA`).
#' @seealso [basic_ugt_inhibition_risk()]
#' @return A markdown-formatted table.
#' @export
#' @noRd
#' @examples
#' basic_ugt_inhibition_risk_table(examplinib_parent, examplinib_ugt_inhibition_data)
basic_ugt_inhibition_risk_table.list <- function(perp, ugt_inh, na.rm=F) {
  for(i in perp) {
    temp <- basic_ugt_inhibition_risk_table.perpetrator(i, ugt_inh=ugt_inh,
           na.rm=na.rm)
    if(!is.null(temp)) {
      print(temp)
    }
  }
}
