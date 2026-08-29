#' Basic CYP inhibition risk
#'
#' Evaluate the clinical risk for direct (reversible) CYP inhibition according
#' to the basic model defined in the [ICH M12
#' guideline](https://www.ema.europa.eu/en/documents/scientific-guideline/ich-m12-guideline-drug-interaction-studies-step-5_en.pdf).
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
#'   the regulatory threshold is reached.
#'
#' @param perp The perpetrator object.
#' @param cyp_inh CYP inhibition data object.
#'
#' @return DDI risk object.
#' @importFrom methods new
#' @export
#' @examples
#' basic_cyp_inhibition_risk(examplinib, examplinib_cyp_inhibition)
#'
basic_cyp_inhibition_risk <- function(perp, cyp_inh) {
  # input validation
  validate_perpetrator(perp)
  validate_inhibition_data(cyp_inh)
  if (perp$name != attr(cyp_inh, "precipitant"))
    warning("Perpetrator name and precipitant do not match!")

  allowed_object <- c("CYP1A2", "CYP2B6", "CYP2C8", "CYP2C9", "CYP2C19",
                      "CYP2D6", "CYP3A4")
  ki <- filter(cyp_inh, .data$object %in% allowed_object)
  if (nrow(ki) == 0)
    stop("No inhibition data for known CYP enzymes found")
  excluded_object <- setdiff(unique(cyp_inh$object), allowed_object)
  if (length(excluded_object) > 0)
    warning(paste0(
      "Non-CYP data were excluded (",
      nice_enumeration(excluded_object),
      ")"
    ))

  # business logic
  out <- ki |>
    mutate(kiu = ki * perp$fumic) |>
    mutate(r = imaxssu(perp, molar = TRUE) / kiu) |>
    mutate(r_gut = case_when(
      .data$object == "CYP3A4" ~ igut(perp, molar = TRUE) / kiu,
      .default = NA)) |>
    mutate(risk_hep = r > 0.02) |>
    mutate(risk_intest = r_gut > 10)  |>
    mutate(r = round(r, digits = 4)) |>
    mutate(r_gut = round(r_gut, digits = 4)) |>
    select(c("object", "ki", "kiu", "source", "r", "risk_hep", "r_gut", "risk_intest"))

  risk(
    out,
    precipitant = perp$name,
    title = paste0("Direct CYP inhibition risk for ", perp$name)
  )
}


#' UGT inhibition risk
#'
#' Evaluate the clinical risk for reversible inhibition of UGT
#' enzymes according to the [ICH M12 guideline](https://www.ema.europa.eu/en/documents/scientific-guideline/ich-m12-guideline-drug-interaction-studies-step-5_en.pdf).
#'
#' @details
#' This function assumes that the UGT inhibition data is provided as
#' \eqn{IC_{50}}. According to
#' [Cheng, Prusoff 1973](https://doi.org/10.1016/0006-2952(73)90196-2)),
#' \eqn{K_i} can be assumed to be \eqn{IC_{50}/2} at the experimental conditions
#' commonly used in the in vitro inhibition studies where substrate
#' concentrations are close to \eqn{K_M}.
#'
#' @details
#' The relevant metric for basic modeling of the UGT inhibition risk is
#' \eqn{R=C_{max,ss,u}/K_{i,u}}
#'
#' Refer to Section 2.1.2.1 of the [ICH M12 guidance document](https://www.ema.europa.eu/en/documents/scientific-guideline/ich-m12-guideline-drug-interaction-studies-step-5_en.pdf))
#' for details.
#'
#' \eqn{R>0.02} are considered to indicate a potential UGT inhibition risk.
#'
#' @param perp The perpetrator object.
#' @param ugt_inh UGT inhibition data object.
#' @return DDI risk object.
#' @importFrom methods new
#' @export
#' @examples
#' basic_ugt_inhibition_risk(examplinib, examplinib_ugt_inhibition)
#'
basic_ugt_inhibition_risk <- function(perp, ugt_inh) {
  # input validation
  validate_perpetrator(perp)
  validate_inhibition_data(ugt_inh)
  if (perp$name != attr(ugt_inh, "precipitant"))
    warning("Perpetrator name and precipitant do not match!")

  allowed_objects <- c("UGT1A1", "UGT1A3", "UGT1A4", "UGT1A6", "UGT1A9",
                      "UGT2B7", "UGT2B15", "UGT2B17")
  ki <- filter(ugt_inh, .data$object %in% allowed_objects)
  if (nrow(ki) == 0)
    stop("No inhibition data for known UGT enzymes found")
  excluded_objects <- setdiff(unique(ugt_inh$object), allowed_objects)
  if (length(excluded_objects) > 0)
    warning(paste0(
      "Non-UGT data were excluded (",
      nice_enumeration(excluded_objects),
      ")"
    ))

  out <- ki |>
    # mutate(ki = ic50/2) |>
    mutate(kiu = .data$ki * perp$fumic) |>
    mutate(r = imaxssu(perp) / .data$kiu) |>
    mutate(risk = .data$r > 0.02) |>
    mutate(r = round(.data$r, digits = 4)) |>
    select(c("object", "ki", "kiu", "source", "r", "risk"))

  risk(
    out,
    precipitant = perp$name,
    title = paste0("UGT inhibition risk for ", perp$name)
  )
}


#' Drug transporter inhibition risk
#'
#' @details
#' The metric for the assessment of transporter interactions is
#' \eqn{R=[I]/IC_{50}}.
#'
#' The relevant perpetrator concentrations \eqn{[I]} and regulatory thresholds
#' of concern are:
#'
#' | I                     | transporter                                                                    | threshold |
#' | ---                   | ---                                                                            | ---       |
#' | \eqn{I_{gut}}         | P-gp and BRCR when drugs are orally administered                               | 10        |
#' | \eqn{C_{max,ss,u}}    | P-gp and BRCR when drugs are administered parenterally or for drug metabolites | 0.02      |
#' | \eqn{I_{max,inlet,u}} | hepatic basolateral transporters OCT1, OATP1B1 and OATP1B3                     | 0.1       |
#' | \eqn{C_{max,ss,u}}    | renal basolateral transporters OAT1, OAT3 and OCT2                             | 0.1       |
#' | \eqn{C_{max,ss,u}}    | apical transporters MATE1 and MATE2-K                                          | 0.02      |
#'
#' @param perp Perpetrator object.
#' @param transporter_inh Inhibitor object.
#' @param transporter_ref Data frame.
#' @param qh Hepatic blood flow in l/min, defaults to 1.616 l/min.
#'
#' @returns DDI risk object.
#' @export
#' @examples
#' transporter_inhibition_risk(examplinib, examplinib_transporter_inhibition)
#'
transporter_inhibition_risk <- function(
    perp,
    transporter_inh,
    transporter_ref = transporter_reference_data,
    qh = 1.616
){
  # input validation
  validate_perpetrator(perp)
  validate_inhibition_data(transporter_inh)

  if (perp$name != attr(transporter_inh, "precipitant"))
    warning("Perpetrator name and precipitant do not match!")

  allowed_objects <- c("Pgp", "BCRP", "OATP1B1", "OATP1B3", "OAT1", "OAT3",
                      "BSEP", "OCT1", "OCT2", "MATE1", "MATE2k")

  in_vitro <- filter(transporter_inh, .data$object %in% allowed_objects)
  if (nrow(in_vitro) == 0)
    stop("No inhibition data for known transporters found")
  excluded_objects <- setdiff(unique(transporter_inh$object), allowed_objects)
  if (length(excluded_objects) > 0)
    warning(paste0(
      "Non-transporter data were excluded (",
      nice_enumeration(excluded_objects),
      ")"
    ))

  perp_conc <- data.frame(
    i = c("igut", "imaxssu", "imaxinletu"),
    conc = c(igut(perp, molar = TRUE), imaxssu(perp, molar = TRUE),
             imaxinletu(perp, qh = qh, molar = TRUE))
  )

  out <- in_vitro |>
    bind_rows(filter(in_vitro, .data$object %in% c("Pgp", "BCRP")) |>
                mutate(object = paste0(.data$object, "_sys"))) |>
    bind_rows(filter(in_vitro, .data$object %in% c("Pgp", "BCRP")) |>
                mutate(object = paste0(.data$object, "_int"))) |>
    filter(!.data$object %in% c("Pgp", "BCRP")) |>
    left_join(
      transporter_ref |>
        mutate(rank = row_number()),
      by = "object") |>
    left_join(perp_conc, by = "i") |>
    mutate(r = case_when(
      is.na(.data$ic50) ~ NA,
      .default = .data$conc / .data$ic50)) |>
    mutate(risk = .data$r >= threshold) |>
    arrange(.data$rank) |>
    select(all_of(c("object", "ic50", "source", "i", "r", "threshold", "risk")))

  risk(out, precipitant = perp$name, title = paste0(
    "Transporter inhibition risk for ", perp$name
  ))
}


#' CYP time-dependent inhibition risk
#'
#' @param perp Perpetrator object.
#' @param cyp_tdi inhibition_data object.
#' @param cyp_kdeg CYP turnover rates as data frame, defaults to
#' [ddir::cyp_turnover].
#'
#' @returns DDI risk object.
#' @export
#' @examples
#' basic_cyp_tdi_risk(examplinib, examplinib_cyp_tdi)
basic_cyp_tdi_risk <- function(
    perp,
    cyp_tdi,
    cyp_kdeg = cyp_turnover
) {
  # input validation
  validate_perpetrator(perp)
  validate_inhibition_data(cyp_tdi)
  expected_columns <- c("object", "ki", "kinact", "source")
  missing_columns <- setdiff(expected_columns, names(cyp_tdi))
  if (length(missing_columns) > 0)
    stop(paste0(
      "Missing columns in cyp_tdi: ", nice_enumeration(missing_columns)))

  allowed_object <- c("CYP1A2", "CYP2B6", "CYP2C8", "CYP2C9", "CYP2C19",
                      "CYP2D6", "CYP3A4")

  in_vitro <- filter(cyp_tdi, .data$object %in% allowed_object)
  if (nrow(in_vitro) == 0)
    stop("No TDI data for known CYP enzymes found")
  excluded_objects <- setdiff(unique(cyp_tdi$object), allowed_object)
  if (length(excluded_objects) > 0)
    warning(paste0(
      "Non-CYP data were excluded (",
      nice_enumeration(excluded_objects),
      ")"
    ))

  # business logic
  out <- in_vitro |>
    mutate(
      kobs = .data$kinact * 5 * imaxssu(perp) /
        (.data$ki * perp$fumic + 5 * imaxssu(perp))
    ) |>
    mutate(fu = perp$fu) |>
    left_join(cyp_kdeg, by = "object") |>
    mutate(kdeg = .data$kdeg_hepatic) |>
    mutate(r = (.data$kobs + .data$kdeg) / .data$kdeg) |>
    mutate(risk = (.data$r >= 1.25)) |>
    select(c("object", "ki", "fu", "kinact", "kdeg", "source", "r", "risk"))

  risk(
    out,
    precipitant = perp$name,
    title = paste0("Time-dependent CYP inhibition risk for ", perp$name)
  )
}


#' Static risk assessment for CYP induction
#'
#' @param perp Perpetrator object.
#' @param cyp_ind Induction object.
#'
#' @returns DDI risk object
#' @export
#' @examples
#' static_cyp_induction_risk(examplinib, examplinib_cyp_induction)
static_cyp_induction_risk <- function(perp, cyp_ind)  {
  # input validation
  validate_perpetrator(perp)
  validate_induction_data(cyp_ind)

  allowed_object <- c("CYP1A2", "CYP2B6", "CYP2C8", "CYP2C9", "CYP2C19",
                      "CYP2D6", "CYP3A4")

  in_vitro <- filter(cyp_ind, .data$object %in% allowed_object)

  if (nrow(in_vitro) == 0)
    stop("No CYP induction data for known CYP enzymes found")

  excluded_object <- setdiff(unique(cyp_ind$object), allowed_object)
  if (length(excluded_object) > 0)
    warning(paste0(
      "Non-CYP data were excluded (",
      nice_enumeration(excluded_object),
      ")"
    ))

  # assess risk
  out <- in_vitro |>
    mutate(maxc_imaxssu = round(.data$max_c / imaxssu(perp), 1)) |>
    mutate(risk = .data$emax >= 2)|>
    mutate(note = case_when(
      .data$maxc_imaxssu < 50 ~ "Not tested up to 50-fold Cmax,u",
      .default = "")) |>
    select(-c("ec50", "maxc_imaxssu"))

  risk(
    out,
    precipitant = perp$name,
    paste0("Static CYP induction risk for ", perp$name)
  )
}


#' Kinetic assessment of CYP induction risk
#'
#' @param perp Perpetrator object
#' @param cyp_ind induction_data object.
#' @param d Scaling factor, defaults to 1.
#'
#' @returns DDI risk object
#' @export
#' @examples
#' kinetic_cyp_induction_risk(examplinib, examplinib_cyp_induction)
#'
kinetic_cyp_induction_risk <- function(perp, cyp_ind, d = 1) {
  # input validation
  validate_perpetrator(perp)
  validate_induction_data(cyp_ind)

  allowed_object <- c("CYP1A2", "CYP2B6", "CYP2C8", "CYP2C9", "CYP2C19",
                      "CYP2D6", "CYP3A4")

  in_vitro <- filter(cyp_ind, .data$object %in% allowed_object)

  if (nrow(in_vitro) == 0)
    stop("No CYP induction data for known CYP enzymes found")

  excluded_object <- setdiff(unique(cyp_ind$object), allowed_object)
  if (length(excluded_object) > 0)
    warning(paste0(
      "Non-CYP data were excluded (",
      nice_enumeration(excluded_object),
      ")"
    ))

  # business logic
  out <- cyp_ind |>
    mutate(r = 1 / (1 + d * .data$emax * 10 * imaxssu(perp) /
                      (.data$ec50 + 10 * imaxssu(perp)))
    ) |>
    mutate(risk = r <= 0.8)

  risk(
    out,
    precipitant = perp$name,
    paste0("Kinetic CYP induction risk for ", perp$name)
  )
}


#' Mechanistic-static risk assessment
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
#' @param perp Perpetrator object.
#' @param cyp_inh Inhibitor object.
#' @param cyp_ind Inducer object.
#' @param cyp_tdi Inhibitor object.
#' @param d Hepatocyte batch scaling factor, defaults to 1.
#' @param include_induction Logical.
#' @param substr The CYP probe substrates to be used as data frame, defaults to
#' `ddir::cyp_reference_substrates`.
#' @param cyp_kdeg The CYP turnover data as data frame. Defaults to the built-in
#' reference data `ddir::cyp_turnover`.
#' @param qh Hepatic blood flow in l/min, defaults to 1.616 l/min.
#' @param qent Enteric blood flow in l/min, defaults to 0.3 l/min (18 l/h).
#'
#' @returns DDI risk object
#' @export
#' @examples
#' mech_stat_cyp_risk(
#'   examplinib,
#'   examplinib_cyp_inhibition,
#'   examplinib_cyp_induction,
#'   examplinib_cyp_tdi
#'  )
#'
mech_stat_cyp_risk <- function(
    perp,
    cyp_inh,
    cyp_ind,
    cyp_tdi = NULL,
    d = 1,
    include_induction = TRUE,
    substr = cyp_reference_substrates,
    cyp_kdeg = cyp_turnover,
    qh = 1.616,
    qent = 18/60
) {
  # input validation
  validate_perpetrator(perp)
  validate_inhibition_data(cyp_inh)
  validate_inhibition_data(cyp_tdi, allow_null = TRUE)
  validate_induction_data(cyp_ind, allow_null = TRUE)

  # risk assessment
  fumic <- perp$fumic
  Ig <- imaxintest(perp, qent = qent)
  Ih <- imaxinletu(perp, qh = qh)

  if(is.null(cyp_ind)) {
    cyp_ind <- induction_data(NULL)
  }

  if(is.null(cyp_tdi)) {
    cyp_tdi <- inhibition_data(NULL) |>
      mutate(kinact = numeric())
  }

  out <- cyp_inh |>
    as.data.frame() |>

    # direct inhibition
    mutate(kiu = .data$ki * fumic) |>
    mutate(Ag = case_when(
      !is.na(.data$ki) ~ 1 / (1 + (Ig / .data$kiu)),
      .default = 1)) |>
    mutate(Ah = case_when(
      !is.na(.data$ki) ~ 1 / (1 + (Ih / .data$kiu)),
      .default = 1)) |>

    # TDI
    left_join(
      cyp_tdi |>
        as.data.frame() |>
        mutate(ki_tdi = .data$ki) |>
        select(-c("ki", "source")),
      by = "object") |>
    left_join(cyp_kdeg, by = "object") |>
    mutate(Bg = case_when(
      !is.na(.data$ki_tdi) ~ .data$kdeg_intestinal /
        (.data$kdeg_intestinal + (Ig * .data$kinact /(Ig + .data$ki_tdi))),
      .default = 1)) |>
    mutate(Bh = case_when(
      !is.na(.data$ki_tdi) ~ .data$kdeg_hepatic /
        (.data$kdeg_hepatic + (Ih * .data$kinact / (Ih + .data$ki_tdi))),
      .default = 1)) |>

    # induction
    left_join(
      cyp_ind |>
        as.data.frame() |>
        select(-"source"),
      by = c("object")) |>
    mutate(Cg = case_when(
      (is.na(.data$ec50) | include_induction == FALSE) ~ 1,
      .default = 1 + (d * .data$emax * Ig / (Ig + .data$ec50)))) |>
    mutate(Ch = case_when(
      (is.na(.data$ec50) | include_induction == FALSE) ~ 1,
      .default = 1 + (d * .data$emax * Ih / (Ih + .data$ec50))))  |>

    # substrate
    left_join(substr, by = "object") |>
    mutate(aucr = 1 / (Ag * Bg * Cg * (1 - .data$fgut) + .data$fgut) *
             1 / (Ah * Bh * Ch * .data$fm * .data$fmcyp + (1 - .data$fm * .data$fmcyp))) |>
    mutate(risk = .data$aucr > 1.25 | .data$aucr < 0.8) |>
    select(c("object", "substrate", "kiu", "fgut", "fm", "fmcyp", "Ag", "Ah",
             "Bg", "Bh", "Cg", "Ch", "aucr", "risk"))

  risk(
    out,
    precipitant = perp$name,
    paste0("Mechanistic-static risk assessment for ", perp$name)
  )
}
