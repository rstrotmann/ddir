#' Plot in vitro CYP induction data
#'
#' @param x The in vitro data set as data frame. The following fields are
#' expected:
#' * DONOR The hepatocyte donor batch as character.
#' * PRECIPITANT The DDI precipitant as character.
#' * CONC The precipitant concentration.
#' * OBJECT The DDI target object as character.
#' * FOLD The mRNA or enzyme activity change at the respective precipitant
#' concentration.
#' * REL The mRNA or enzyme activity induction, as percent of positive control
#' at the respective precipitant concentration.
#'
#' FOLD and REL are optional but a column that corresponds to the 'type'
#' argument must be present.
#' @param type The DDI metric, e.g. FOLD or REL, as character. A different type
#' can be given but a column of that name must be present in the data set.
#' @param log Plot logarithmic x axis, as logical.
#'
#' @returns A ggplot object.
#' @import ggplot2
#' @export
#'
#' @examples
#' induction_plot(examplinib_in_vitro_ind)
#' induction_plot(examplinib_in_vitro_ind, type = "FOLD")
induction_plot <- function(
    x,
    type = "REL",
    log = TRUE
) {
  # input validation
  validate_argument(type, values = c("FOLD", "REL"))
  validate_df_argument(
    x, expected_fields = c("DONOR", "SAMPLE", "CONC", "OBJECT", type))
  validate_argument(log, "logical")

  # business logic
  p <- x |>
    filter(.data$SAMPLE == "test") |>
    filter(!is.na(.data[[type]])) |>
    ggplot(aes(
      x = .data$CONC,
      y = .data[[type]],
      group = .data$DONOR,
      color = .data$DONOR)
    ) +
    geom_point() +
    geom_line() +
    facet_wrap(~OBJECT, scales = "free") +
    theme_bw()

  if (isTRUE(log))
    p <- p + scale_x_log10()

  return(p)
}


#' Flag concentration-dependent CYP mRNA down-regulation
#'
#' @description
#' `r lifecycle::badge("experimental")`
#'
#' For each CYP isoform (and, by default, each hepatocyte donor), decide
#' whether vehicle-normalized mRNA fold-change falls as precipitant
#' concentration rises, rather than being induced.
#'
#' @details
#' Only `SAMPLE == "test"` rows with `CONC > 0` and non-missing `FOLD` are
#' used. A curve is flagged when both of the following hold:
#'
#' * fold-change at the highest tested concentration is below
#'   `fold_threshold` (default 0.5, i.e. ≤ 50% of vehicle);
#' * Spearman's rank correlation of `CONC` with `FOLD` is negative, so
#'   fold-change tends to decrease as concentration increases.
#'
#' `FOLD` is fold vs vehicle (1 = no change, &lt; 1 = suppression). Percent
#' of positive control (`REL`) is not used. Isolated `FOLD < 1` points at
#' low concentration, which occur in ordinary induction curves, do not
#' trigger the flag.
#'
#' The rank correlation uses order only, not a linear fit on the µM scale.
#' No p-value is required; typical induction experiments have too few
#' concentrations for a Spearman test to be informative.
#'
#' @param x An `induction_experiment` object. A `FOLD` column is required.
#' @param fold_threshold Maximum fold-change at the highest concentration
#'   that still counts as suppression, as numeric. Defaults to 0.5.
#' @param by_donor If `TRUE` (default), one row per isoform and donor. If
#'   `FALSE`, donors are pooled and there is one row per isoform.
#'
#' @returns A data frame with one row per isoform, or per isoform and
#'   donor, and the columns:
#' * `OBJECT` CYP isoform.
#' * `DONOR` Hepatocyte donor (omitted if `by_donor = FALSE`).
#' * `n` Number of concentrations used.
#' * `min_fold`, `max_fold` Range of `FOLD`.
#' * `fold_maxc` `FOLD` at the highest `CONC`.
#' * `rho` Spearman correlation of `CONC` and `FOLD`.
#' * `downregulation` `TRUE` if the curve meets both criteria above.
#'
#' @importFrom stats cor
#' @export
#'
#' @examples
#' induction_downregulation(examplinib_in_vitro_ind)
#' induction_downregulation(examplinib_in_vitro_ind1)
induction_downregulation <- function(x, fold_threshold = 0.5, by_donor = TRUE) {
  if (!inherits(x, "induction_experiment")) {
    stop("x must be an induction_experiment object")
  }
  validate_argument(fold_threshold, "numeric")
  validate_argument(by_donor, "logical")
  if (!"FOLD" %in% names(x)) {
    stop("FOLD is required to assess down-regulation")
  }

  groups <- if (isTRUE(by_donor)) {
    c("OBJECT", "DONOR")
  } else {
    "OBJECT"
  }

  x |>
    filter(.data$SAMPLE == "test", .data$CONC > 0, !is.na(.data$FOLD)) |>
    reframe(
      n = n(),
      min_fold = min(.data$FOLD),
      max_fold = max(.data$FOLD),
      fold_maxc = .data$FOLD[which.max(.data$CONC)],
      rho = suppressWarnings(
        cor(.data$CONC, .data$FOLD, method = "spearman")
      ),
      downregulation = .data$fold_maxc < fold_threshold & .data$rho < 0,
      .by = all_of(groups)
    ) |>
    arrange(.data$OBJECT, .data$DONOR)
}


#' Analyze in vitro CYP induction data
#'
#' Fit Emax / EC50 curves to in vitro CYP mRNA induction
#'
#' @description
#' `r lifecycle::badge("experimental")`
#'
#' For each hepatocyte donor and CYP isoform, fit a hyperbolic Emax model
#' to vehicle-normalized mRNA fold-change (`FOLD`) from an
#' `induction_experiment`.
#'
#' @details
#' Only `SAMPLE == "test"` rows with `CONC > 0` are used. The model is
#'
#' \deqn{f(C) = 1 + \frac{E_{max} \cdot C}{EC_{50} + C}}
#'
#' which is an Emax relationship with baseline 1 (vehicle) and Hill slope
#' 1.
#' Emax is **fold-increase** (fold-change minus 1), not
#' fold-change. Observed Emax is `max(FOLD) - 1`.
#'
#' By default (`use_emax_obs = TRUE`) that observed value is held fixed
#' and only EC50 is estimated. If `use_emax_obs = FALSE`, both
#' Emax and EC50 are estimated.
#'
#' The fit assumes fold-change rises with concentration. Concentration-
#' dependent down-regulation (see [ddir::induction_downregulation()]) is not
#' modelled and will yield unusable parameters.
#'
#' `static_cyp_induction_risk()` treats `emax >= 2` as ICH ≥ 2-fold
#' induction. Values from this function are fold-increase, so a 2.4-fold
#' curve has `emax` 1.4 and would be called no risk if passed through
#' unchanged.
#'
#' @param x An `induction_experiment` object. A `FOLD` column is required.
#' @param use_emax_obs If `TRUE` (default), constrain \(E_{max}\) to the
#'   observed maximum fold-increase. If `FALSE`, \(E_{max}\) is fitted.
#' @param individual_donors If `TRUE` (default), `fold_plot` facets by
#'   donor and isoform. If `FALSE`, facets by isoform and colours donors.
#'
#' @returns A named list:
#' * `data` Nested tibble: raw curve, `nlsLM` fit, tidy coefficients, and
#'   `emax_obs` for each donor and isoform.
#' * `fold_plot` ggplot of observed `FOLD` and the fitted curves.
#' * `ind_param` Tidy coefficient table (`term`, `estimate`, …) by
#'   isoform and donor. With the default `use_emax_obs = TRUE` the only
#'   fitted term is `ec50`; `emax_obs` is a separate column.
#' * `downregulation` Tibble flagging donor/object combinations with suspected
#'   downregulation.
#'
#' @seealso [ddir::induction_experiment()], [ddir::induction_downregulation()]
#'
#' @import ggplot2
#' @export
#'
#' @examples
#' indmod(examplinib_in_vitro_ind)
#'
#' @examples
#' indmod(induction_experiment(examplinib_in_vitro_ind, "examplinib"))
#' indmod(induction_experiment(examplinib_in_vitro_ind1, "examplinib"))
indmod <- function(
    x,
    use_emax_obs = TRUE,
    individual_donors = TRUE
) {
  # input validation
  if (!inherits(x, "induction_experiment"))
    stop("input must be an in vitro induction experiment object!")

  validate_argument(use_emax_obs, "logical")
  validate_argument(individual_donors, "logical")

  # business logic
  data <- mutate(x, ID = paste0(.data$OBJECT, "_", .data$DONOR)) |>
    filter(.data$CONC > 0)

  precipitant <- attr(x, "precipitant")

  sigm <- function(c, emax, ec50) {
    1 + (emax / (1 + exp(log(ec50) - log(c))))
  }

  out <- list()

  # down-regulation check
  downreg <- induction_downregulation(x)
  out$downregulation <- downreg

  if (any(downreg$downregulation)) {
    temp <- downreg |>
      filter(.data$downregulation) |>
      distinct(.data$OBJECT)
    warning(paste(
      "Check for down-regulation of",
      nice_enumeration(temp$OBJECT),
      "by", precipitant
    ))
  }

  # prepare data set
  temp <- data |>
    filter(.data$SAMPLE == "test") |>
    nest_by(.data$DONOR, .data$OBJECT, .data$ID) |>
    mutate(emax_obs = max(data$FOLD, na.rm = TRUE) - 1)

  if (use_emax_obs == TRUE) {
    out$data <- temp |>
      mutate(mod = list(
        nlsLM(
          FOLD ~ sigm(CONC, emax_obs, ec50),
          data = data,
          start = list(ec50 = .1),
          lower = c(ec50 = 0),
          upper = c(ec50 = 100),
          control = nls.lm.control(maxiter = 1000)
        )
      )) |>
      mutate(modpar = list(broom::tidy(.data$mod)))
  } else {
    out$data <- temp |>
      mutate(mod = list(
        nlsLM(
          FOLD ~ sigm(CONC, emax, ec50),
          data = data,
          start = list(emax = 2, ec50 = .1),
          lower = c(emax = NA, ec50 = 0),
          upper = c(emax = NA, ec50 = 100),
          control = nls.lm.control(maxiter = 1000)
        )
      )) |>
      mutate(modpar = list(broom::tidy(.data$mod)))
  }

  # make curve plot data set
  pred <- data.frame(
    CONC = 10^seq(
      log10(min(data$CONC, na.rm = TRUE)),
      log10(max(data$CONC, na.rm = TRUE)),
      length.out = 100)
  )

  temp <- lapply(out$data$mod, function(x) stats::predict(x, newdata = pred))
  names(temp) <- out$data$ID
  pred <- bind_cols(pred, temp) |>
    pivot_longer(cols = -1, names_to = "ID", values_to = "FOLD") |>
    separate(.data$ID, c("OBJECT", "DONOR"), "_", remove = FALSE)

  # make output plot
  out$fold_plot <- ggplot(data = NULL, aes(
      x = .data$CONC,
      y = .data$FOLD,
      color = .data$OBJECT,
      group = .data$ID
    )) +
    geom_line(data = pred) +
    geom_point(data = filter(data, .data$SAMPLE == "test", !is.na(.data$FOLD)), size = 2) +
    scale_x_log10() +
    expand_limits(y = 0) +
    labs(title = paste("In vitro CYP induction by", precipitant)) +
    theme_bw()

  if (isTRUE(individual_donors)) {
    out$fold_plot <- out$fold_plot +
      facet_wrap(~ID, scales = "free") +
      theme(legend.position = "none")
  } else {
    out$fold_plot <- out$fold_plot +
      aes(color = .data$DONOR) +
      facet_wrap(~OBJECT, scales = "free") +
      theme(legend.position = "bottom")
  }

  out$ind_param <- out$data |>
    unnest(.data$modpar) |>
    ungroup() |>
    select(-c("ID", "data", "mod")) |>
    relocate("OBJECT", "DONOR", "term") |>
    arrange(.data$OBJECT, .data$DONOR, .data$term)

  # make inducer object from the donor that has the respective highers emax
  max_c <- data |>
    reframe(
      max_c = max(.data$CONC, na.rm = TRUE), .by = c("ID", "SOURCE"))

  temp <- out$data |>
    unnest(.data$modpar) |>
    select(-c("data", "mod")) |>
    group_by(.data$OBJECT) |>
    filter(.data$emax_obs == max(.data$emax_obs, na.rm = TRUE)) |>
    filter(.data$term == "ec50") |>
    left_join(max_c, by = "ID") |>
    select(
      object = .data$OBJECT,
      emax = .data$emax_obs,
      ec50 = .data$estimate,
      max_c,
      source = .data$SOURCE
    ) |>
    mutate(ec50 = round(.data$ec50, 2))

  out
}
