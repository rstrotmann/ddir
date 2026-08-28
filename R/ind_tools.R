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


#' Analyze in vitro CYP induction data
#'
#' Fit a 3-parameter hill function to the mRNA induction data from each
#' donor and DDI object.
#'
#' @details
#'
#' @description
#' `r lifecycle::badge("experimental")`
#'
#' \deqn{f = 1 + \frac{(E_{max}) * C^n}{(EC_{50}^n + C^n)}}
#'
#' @param x In vitro induction data as induction_experiment object.
#' @param use_emax_obs Constrain Emax to the observed Emax (default). If FALSE,
#' Emax will be fitted.
#' @param individual_donors Plot donors on individual panels.
#'
#' @returns A list with the elements:
#' * data The original data, model and model parameters for each fit.
#' * fold_plot A ggplot object showing the original data and the fits.
#' * ind_param The ec50, hill (n) and emax parameters by donor and DDI object.
#' * inducer An inducer object representing the model parameters from the donor
#' with the highest Emax per DDI object.
#'
#' @import ggplot2
#' @export
#'
#' @examples
#' indmod(induction_experiment(examplinib_in_vitro_ind, "examplinib"))
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
  data <- mutate(x@data, ID = paste0(.data$OBJECT, "_", .data$DONOR)) |>
    filter(CONC > 0)

  # sigm <- function(c, emax, ec50, h) {
  #   1 + (emax - 1) / (1 + exp(-(log(c) - log(ec50))/h))
  # }

  sigm <- function(c, emax, ec50) {
    1 + (emax / (1 + exp(log(ec50) - log(c))))
  }

  out <- list()

  # prepare data set
  temp <- data |>
    filter(.data$SAMPLE == "test") |>
    nest_by(.data$DONOR, .data$OBJECT, .data$ID) |>
    mutate(emax_obs = max(data$FOLD, na.rm = TRUE) - 1)

  # non-linear modeling
  # if (use_emax_obs == TRUE) {
  #   out$data <- temp |>
  #     mutate(mod = list(
  #       nlsLM(
  #         FOLD ~ sigm(CONC, emax_obs, ec50, n),
  #         data = data,
  #         start = list(ec50 = .1, n = 1),
  #         lower = c(ec50 = 0, n = 1),
  #         upper = c(ec50 = 100, n = 5),
  #         control = nls.lm.control(maxiter = 1000)
  #       )
  #     )) |>
  #     mutate(modpar = list(broom::tidy(.data$mod)))
  # } else {
  #   out$data <- temp |>
  #     mutate(mod = list(
  #       nlsLM(
  #         FOLD ~ sigm(CONC, emax, ec50, n),
  #         data = data,
  #         start = list(emax = 2, ec50 = .1, n = 1),
  #         lower = c(emax = NA, ec50 = 0, n = 1),
  #         upper = c(emax = NA, ec50 = 100, n = 5),
  #         control = nls.lm.control(maxiter = 1000)
  #       )
  #     )) |>
  #     mutate(modpar = list(broom::tidy(.data$mod)))
  # }

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
    labs(title = paste("In vitro CYP induction by", x@precipitant)) +
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

  # out$inducer <- inducer(temp, precipitant = x@precipitant)

  out

}
