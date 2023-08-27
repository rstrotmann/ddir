
transporter_ref <- data.frame(
  param=c("Pgp_int", "Pgp_sys", "BCRP_int", "BCRP_sys", "OCT1", "OATP1B1",
          "OATP1B3", "OAT1", "OAT3", "BSEP", "OCT2", "MATE1", "MATE2k"),
  rank=seq(1, 13),
  fda_thld=c(10, 0.1, 10, 0.1, NA, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1),
  ema_thld=c(10, 0.02, 10, 0.02, 0.04, 0.04, 0.04, 0.04, 0.04, 0.02, 0.02, 0.02, 0.02),
  i=c("igut", "imaxssu", "igut", "imaxssu", "imaxinletu", "imaxinletu",
      "imaxinletu", "imaxssu", "imaxssu", "imaxssu", "imaxssu", "imaxssu", "imaxssu")
)



#' Title
#'
#' @param perp
#' @param transporter_ic50
#'
#' @return
#' @export
#'
#' @examples
transporter_inhibition_risk <- function(perp, transporter_ic50) {
  ic50 <- transporter_ic50 %>%
    filter(name==name(perp)) %>%
    filter(param!="name") %>%
    select(-name) %>%
    mutate(value=as.num(value))

  temp <-  key_concentrations(perp, molar=TRUE)
  i <- data.frame(
    i=names(temp),
    conc=temp)

  # duplicate rows Pgp and BCRP, if applicable, and assign intestinal and
  #   systemic scope
  out <- ic50 %>%
    bind_rows(filter(ic50, param %in% c("Pgp", "BCRP")) %>%
                mutate(param=paste0(param, "_sys"))) %>%
    bind_rows(filter(ic50, param %in% c("Pgp", "BCRP")) %>%
                mutate(param=paste0(param, "_int"))) %>%
    filter(!param %in% c("Pgp", "BCRP")) %>%
    left_join(transporter_ref, by="param") %>%
    left_join(i, by="i") %>%
    mutate(r=case_when(is.na(value) ~ NA, .default=conc/value)) %>%
    mutate(fda_risk=r>fda_thld) %>%
    mutate(ema_risk=r>ema_thld) %>%
    mutate(ic50=value) %>%
    arrange(rank) %>%
    select(param, ic50, source, r, fda_thld, fda_risk, ema_thld, ema_risk)
  return(out)
}

transporter_inhibition_risk_table <- function(perp, transporter_ic50, na.rm=F) {
  temp <- transporter_inhibition_risk(perp, transporter_ic50)

  if(na.rm==TRUE) {
    temp <- temp %>%
      filter(!is.na(ic50))
  }

  labels <- c("transporter", "$IC_{50}$", "source", "$R$", "thld FDA", "risk FDA",
              "thld EMA", "risk EMA")
  if(nrow(temp)!=0) {
    out <- knitr::kable(
      temp, caption=paste("Risk for drug transporter inhibition by",
      name(perp)), col.names=labels)
    return(out)
  }
}
