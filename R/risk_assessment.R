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
#' @return DDI risk object.
#' @export
basic_cyp_inhibition_risk <- function(perp, cyp_inh) {
  # Validate inputs
  if (!inherits(perp, "perpetrator")) {
    stop("perp must be a perpetrator object")
  }
  if (!inherits(cyp_inh, "inhibitor")) {
    stop("cyp_inh must be an inhibitor object")
  }

  allowed_target <- c("CYP1A2", "CYP2B6", "CYP2C8", "CYP2C9", "CYP2C19",
                      "CYP2D6", "CYP3A4")
  ki <- filter(cyp_inh@data, target %in% allowed_target)
  if (nrow(ki) == 0)
    stop("No inhibition data for known CYP enzymes found")
  excluded_target <- setdiff(unique(cyp_inh@data$target), allowed_target)
  if (length(excluded_target) > 0)
    warning(paste0(
      "Non-CYP data were excluded (",
      nice_enumeration(excluded_target),
      ")"
    ))

  # risk assessmen
  out <- ki |>
    mutate(kiu = ki * perp@fumic) |>
    mutate(r = imaxssu(perp, molar = TRUE) / kiu) |>
    mutate(r_gut = case_when(
      target == "CYP3A4" ~ igut(perp, molar = TRUE) / kiu,
      .default = NA)) |>
    mutate(risk_hep = r > 0.02) |>
    mutate(risk_intest = r_gut > 10)  |>
    mutate(r = round(r, digits = 4)) |>
    mutate(r_gut = round(r_gut, digits = 4)) |>
    select(target, ki, kiu, source, r, risk_hep, r_gut, risk_intest)

  new(
    "ddi_risk",
    table = out,
    object = perp,
    title = paste0(
      "Direct CYP inhibition risk for ", perp@name
    ))
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
#' @export
basic_ugt_inhibition_risk <- function(perp, ugt_inh) {
  # Validate inputs
  if (!inherits(perp, "perpetrator")) {
    stop("perp must be a perpetrator object")
  }
  if (!inherits(ugt_inh, "inhibitor")) {
    stop("ugt_inh must be an inhibitor object")
  }

  allowed_target <- c("UGT1A1", "UGT1A3", "UGT1A4", "UGT1A6", "UGT1A9",
                      "UGT2B7", "UGT2B15", "UGT2B17")
  ki <- filter(ugt_inh@data, target %in% allowed_target)
  if (nrow(ki) == 0)
    stop("No inhibition data for known UGT enzymes found")
  excluded_target <- setdiff(unique(ugt_inh@data$target), allowed_target)
  if (length(excluded_target) > 0)
    warning(paste0(
      "Non-UGT data were excluded (",
      nice_enumeration(excluded_target),
      ")"
    ))

  out <- ki |>
    mutate(kiu = ki * perp@fumic) |>
    mutate(r = imaxssu(perp) / kiu) |>
    mutate(risk = r > 0.02) |>
    mutate(r = round(r, digits = 4)) |>
    select(target, ki, kiu, source, r, risk)

  new(
    "ddi_risk",
    table = out,
    object = perp,
    title = paste0(
      "UGT inhibition risk for ", perp@name
    ))
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
#'
#' @returns DDI risk object.
#' @export
transporter_inh_risk <- function(
    perp,
    transporter_inh,
    transporter_ref = transporter_reference_data,
    qh = 1.616
){
  # input validation
  if (!inherits(perp, "perpetrator")) {
    stop("perp must be a perpetrator object")
  }
  if (!inherits(transporter_inh, "inhibitor")) {
    stop("transporter_inh must be an inhibitor object")
  }

  allowed_target <- c("Pgp", "BCRP", "OATP1B1", "OATP1B3", "OAT1", "OAT3",
                      "BSEP", "OCT1", "OCT2", "MATE1", "MATE2k")

  in_vitro <- filter(transporter_inh@data, target %in% allowed_target)
  if (nrow(in_vitro) == 0)
    stop("No inhibition data for known transporters found")
  excluded_target <- setdiff(unique(transporter_inh@data$target), allowed_target)
  if (length(excluded_target) > 0)
    warning(paste0(
      "Non-transporter data were excluded (",
      nice_enumeration(excluded_target),
      ")"
    ))

  perp_conc <- data.frame(
    i = c("igut", "imaxssu", "imaxinletu"),
    conc = c(igut(perp, molar = TRUE), imaxssu(perp, molar = TRUE),
             imaxinletu(perp, qh = qh, molar = TRUE))
  )

  transporter_ref <- rename(transporter_ref, target = transporter)

  out <- in_vitro %>%
    bind_rows(filter(in_vitro, target %in% c("Pgp", "BCRP")) %>%
                mutate(target = paste0(target, "_sys"))) %>%
    bind_rows(filter(in_vitro, target %in% c("Pgp", "BCRP")) %>%
                mutate(target = paste0(target, "_int"))) %>%
    filter(!target %in% c("Pgp", "BCRP")) %>%
    left_join(
      transporter_ref |>
        mutate(rank = row_number()),
      by = "target") |>
    left_join(perp_conc, by = "i") |>
    mutate(r = case_when(
      is.na(ic50) ~ NA,
      .default = conc / ic50)) |>
    mutate(risk = r >= threshold) |>
    arrange(rank) |>
    select(target, ic50, source, i, r, threshold, risk)

  new(
    "ddi_risk",
    table = out,
    object = perp,
    title = paste0(
      "Transporter inhibition risk for ", perp@name
    ))
}


#' Title
#'
#' @param perp Perpetrator object.
#' @param cyp_tdi Inhibitor object.
#' @param cyp_kdeg Data frame
#'
#' @returns DDI risk object
#' @export
basic_cyp_tdi_risk <- function(
    perp,
    cyp_tdi,
    cyp_kdeg = cyp_turnover
) {
  # input validation
  if (!inherits(perp, "perpetrator")) {
    stop("perp must be a perpetrator object")
  }
  if (!inherits(cyp_tdi, "inhibitor")) {
    stop("cyp_tdi must be an inhibitor object")
  }
  expected_columns <- c("target", "ki", "kinact", "source")
  missing_columns <- setdiff(expected_columns, names(cyp_tdi@data))
  if (length(missing_columns) > 0)
    stop(paste0(
      "Missing columns in cyp_tdi: ", nice_enumeration(missing_columns)))

  allowed_target <- c("CYP1A2", "CYP2B6", "CYP2C8", "CYP2C9", "CYP2C19",
                      "CYP2D6", "CYP3A4")

  in_vitro <- filter(cyp_tdi@data, target %in% allowed_target)
  if (nrow(in_vitro) == 0)
    stop("No TDI data for known CYP enzymes found")
  excluded_target <- setdiff(unique(cyp_tdi@data$target), allowed_target)
  if (length(excluded_target) > 0)
    warning(paste0(
      "Non-CYP data were excluded (",
      nice_enumeration(excluded_target),
      ")"
    ))

  cyp_kdeg <- cyp_kdeg |>
    rename(target = cyp)

  out <- cyp_tdi@data |>
    mutate(kobs = .data$kinact * 5 * imaxssu(perp) / (.data$ki * perp@fumic + 5 * imaxssu(perp))) |>
    mutate(fu = perp@fu) |>
    left_join(cyp_kdeg, by = "target") |>
    mutate(kdeg = kdeg_hepatic) |>
    mutate(r = (kobs + kdeg) / kdeg) |>
    mutate(risk = (r > 1.25)) |>
    select(target, ki, fu, kinact, kdeg, source, r, risk)

  new(
    "ddi_risk",
    table = out,
    object = perp,
    title = paste0(
      "Time-dependent CYP inhibition risk for ", perp@name
    ))
}


#' Title
#'
#' @param perp Perpetrator object
#' @param cyp_ind Induction object.
#'
#' @returns DDI risk object
#' @export
static_cyp_induction_risk <- function(perp, cyp_ind)  {
  # input validation
  if (!inherits(perp, "perpetrator")) {
    stop("perp must be a perpetrator object")
  }
  if (!inherits(cyp_ind, "inducer")) {
    stop("cyp_ind must be an induction object")
  }
  expected_columns <- c("target", "emax", "ec50", "max_c", "source")
  missing_columns <- setdiff(expected_columns, names(cyp_ind@data))
  if (length(missing_columns) > 0)
    stop(paste0(
      "Missing columns in cyp_ind: ", nice_enumeration(missing_columns)))

  allowed_target <- c("CYP1A2", "CYP2B6", "CYP2C8", "CYP2C9", "CYP2C19",
                      "CYP2D6", "CYP3A4")

  in_vitro <- filter(cyp_ind@data, target %in% allowed_target)
  if (nrow(in_vitro) == 0)
    stop("No CYP induction data for known CYP enzymes found")
  excluded_target <- setdiff(unique(cyp_ind@data$target), allowed_target)
  if (length(excluded_target) > 0)
    warning(paste0(
      "Non-CYP data were excluded (",
      nice_enumeration(excluded_target),
      ")"
    ))

  # assess risk
  out <- cyp_ind@data |>
    mutate(maxc_imaxssu = round(max_c / imaxssu(perp), 1)) |>
    mutate(risk = emax >= 2)|>
    mutate(note = case_when(
      maxc_imaxssu < 50 ~ "Not tested up to 50-fold Cmax,u",
      .default = "")) |>
    select(-c(ec50, maxc_imaxssu))

  new(
    "ddi_risk",
    table = out,
    object = perp,
    title = paste0(
      "Static CYP induction risk for ", perp@name
    ))
}


#' Title
#'
#' @param perp Perpetrator object
#' @param cyp_ind Induction object
#' @param d Numeric
#'
#' @returns DDI risk object
#' @export
kinetic_cyp_induction_risk <- function(perp, cyp_ind, d=1) {
  # input validation
  if (!inherits(perp, "perpetrator")) {
    stop("perp must be a perpetrator object")
  }
  if (!inherits(cyp_ind, "inducer")) {
    stop("cyp_ind must be an inducer object")
  }
  expected_columns <- c("target", "emax", "ec50", "max_c", "source")
  missing_columns <- setdiff(expected_columns, names(cyp_ind@data))
  if (length(missing_columns) > 0)
    stop(paste0(
      "Missing columns in cyp_ind: ", nice_enumeration(missing_columns)))

  allowed_target <- c("CYP1A2", "CYP2B6", "CYP2C8", "CYP2C9", "CYP2C19",
                      "CYP2D6", "CYP3A4")

  in_vitro <- filter(cyp_ind@data, target %in% allowed_target)
  if (nrow(in_vitro) == 0)
    stop("No CYP induction data for known CYP enzymes found")
  excluded_target <- setdiff(unique(cyp_ind@data$target), allowed_target)
  if (length(excluded_target) > 0)
    warning(paste0(
      "Non-CYP data were excluded (",
      nice_enumeration(excluded_target),
      ")"
    ))

  out <- cyp_ind@data |>
    mutate(r = (1 / (1 + d * emax * 10 * imaxssu(perp) / (ec50 + 10 * imaxssu(perp))))) |>
    mutate(risk = r <= 0.8)

  new(
    "ddi_risk",
    table = out,
    object = perp,
    title = paste0(
      "Kinetic CYP induction risk for ", perp@name
    ))
}


#' Title
#'
#' @param perp Perpetrator object
#' @param cyp_inh Inhibitor object
#' @param cyp_ind Inducer object
#' @param cyp_tdi Inhibitor object
#' @param d Numeric
#' @param include_induction Logical
#' @param substr Data frame
#' @param cyp_kdeg Data frame
#'
#' @returns DDI risk object
#' @export
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
  if (!inherits(perp, "perpetrator")) {
    stop("perp must be a perpetrator object")
  }
  if (!inherits(cyp_inh, "inhibitor")) {
    stop("cyp_inh must be an inhibitor object")
  }
  if (!inherits(cyp_ind, "inducer")) {
    stop("cyp_ind must be an inducer object")
  }
  if (!inherits(cyp_tdi, "inhibitor")) {
    stop("cyp_tdi must be an inhibitor object")
  }

  # allowed_target <- c("CYP1A2", "CYP2B6", "CYP2C8", "CYP2C9", "CYP2C19",
  #                     "CYP2D6", "CYP3A4")
  # ki <- filter(cyp_inh@data, target %in% allowed_target)
  # if (nrow(ki) == 0)
  #   stop("No inhibition data for known CYP enzymes found")
  #
  # excluded_target <- setdiff(unique(cyp_inh@data$target), allowed_target)
  # if (length(excluded_target) > 0)
  #   warning(paste0(
  #     "Non-CYP data were excluded (",
  #     nice_enumeration(excluded_target),
  #     ")"
  #   ))

  # risk assessment
  fumic <- perp@fumic
  Ig <- imaxintest(perp, qent = qent)
  Ih <- imaxinletu(perp, qh = qh)

  if(is.null(cyp_ind)) {
    cyp_ind <- inducer(NULL)
  }

  if(is.null(cyp_tdi)) {
    cyp_tdi <- inhibitor(NULL)
  }

  cyp_kdeg <- rename(cyp_kdeg, target = cyp)
  substr <- rename(substr, target = cyp)

  out <- cyp_inh@data |>
    select(-source) |>

    # direct inhibition
    mutate(kiu = ki * fumic) |>
    mutate(Ag = case_when(
      !is.na(ki) ~ 1 / (1 + (Ig / kiu)),
      .default = 1)) |>
    mutate(Ah = case_when(
      !is.na(ki) ~ 1 / (1 + (Ih / kiu)),
      .default = 1)) |>

    # TDI
    left_join(
      cyp_tdi@data %>%
        mutate(ki_tdi = ki) |>
        select(-c(ki, source)),
      by = "target") |>
    left_join(cyp_kdeg, by = "target") |>
    mutate(Bg = case_when(
      !is.na(ki_tdi) ~ kdeg_intestinal / (kdeg_intestinal +
                                            (Ig * kinact /(Ig + ki_tdi))),
      .default = 1)) |>
    mutate(Bh = case_when(
      !is.na(ki_tdi) ~ kdeg_hepatic / (kdeg_hepatic +
                                         (Ih * kinact / (Ih + ki_tdi))),
      .default = 1)) |>

    # induction
    left_join(
      cyp_ind@data |>
        select(-source),
      by=c("target")) |>
    mutate(Cg = case_when(
      (is.na(ec50) | include_induction == FALSE) ~ 1,
      .default = 1 + (d * emax * Ig / (Ig + ec50)))) |>
    mutate(Ch = case_when(
      (is.na(ec50) | include_induction == FALSE) ~ 1,
      .default = 1 + (d * emax * Ih / (Ih + ec50))))  |>

    # substrate
    left_join(substr, by="target") |>
    mutate(aucr = 1 / (Ag * Bg * Cg * (1 - fgut) + fgut) *
             1 / (Ah * Bh * Ch * fm * fmcyp + (1 - fm * fmcyp))) |>
    mutate(risk = aucr > 1.25 | aucr < 0.8) %>%
    select(c("target", "substrate", "kiu", "fgut", "fm", "fmcyp", "Ag", "Ah",
             "Bg", "Bh", "Cg", "Ch", "aucr", "risk"))

  new(
    "ddi_risk",
    table = out,
    object = perp,
    title = paste0(
      "Mechanistic-static CYP inhibition risk for ", perp@name
    ))
}
