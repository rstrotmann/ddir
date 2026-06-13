#' Derive TDI parameters from enzyme activity data
#'
#' @details
#' kobs is fitted to the enzyme activity data as a first-order process for each
#' inhibitor concentration.
#' kinact (the maximal inhibition rate) and KI (inhibitor concentration at the
#' half-maximal kobs) are fitted from kobs over inhibitor concentration using a
#' Emax model.
#'
#' @param x In vitro data as data frame. The following columns are expected:
#' * TIME Duration of Pre-incubation with the precipitant and NADPH in minutes
#' * CONC Precipitant concentration in uM
#' * ACT Enzyme activity in percent
#' @param object The DDI target object as character or NULL.
#'
#' @returns A list with the following items:
#' * data The input data with linear models of ACT over TIME by CONC.
#' * kobs The kobs parameters by precipitant concentration.
#' * kobs_plot A ggplot object showing the kobs fit to the ACT data.
#' * tdi_param The TDI parameters, Kinact and KI, derived from an Emax model of
#'   kobs over CONC.
#' * tdi_plot A ggplot object showing the Emnax modeling of kobs over CONC.
#' @import ggplot2
#' @import tidyr
#' @import minpack.lm
#' @import broom
#' @export
#'
#' @examples
#' tdimod(examplinib_in_vitro_tdi)
tdimod <- function(x, object = NULL) {
  # input validation
  if (!is.data.frame(x))
    stop("x must be a data frame")

  expected_cols <- c("TIME", "CONC", "ACT")
  missing_cols <- setdiff(expected_cols, names(x))
  if (length(missing_cols) > 0)
    stop(paste0("Missing column(s) in input: ", nice_enumeration(missing_cols)))

  validate_argument(object, "character", allow_null = TRUE)

  # business logic
  out <- list()

  # derive kobs
  out$data <- x |>
    nest_by(CONC) |>
    mutate(mod = list(lm(log(ACT) ~ TIME, data = data))) |>
    mutate(modpar = list(broom::tidy(mod)))

  temp <- out$data |>
    unnest(modpar)

  out$kobs <- temp |>
    filter(term == "TIME") |>
    mutate(kobs = -estimate) |>
    select(CONC, kobs, std.error, p.value)

  # kobs_plot
  pred <- temp |>
    select(CONC, term, estimate) |>
    pivot_wider(names_from = "term", values_from = "estimate")

  out$kobs_plot <- x |>
    arrange(CONC, TIME) |>
    mutate(ACT = log(ACT)) |>
    ggplot(aes(x = TIME, y = ACT, group = as.factor(CONC), color = as.factor(CONC))) +
    geom_point(size = 2) +
    geom_line(linetype = "dashed") +
    labs(x = "time (min)", y = "ln(activity)", color = "concentration (µM)") +
    geom_abline(
      aes(color = as.factor(CONC), slope = TIME, intercept = `(Intercept)`),
      data = pred,
      show.legend = FALSE) +
    theme_bw() +
    theme(legend.position = "bottom")

  # TDI parameters from Emax model
  emax <- function(x, kinact, kI) {
    return(kinact * x / (kI + x))
  }

  if (nrow(out$kobs) < 3) {
    warning(
      "At least 3 concentration levels are needed to fit Kinact and KI; ",
      "skipping Emax fit.",
      call. = FALSE
    )
    fit <- NULL
  } else {
    fit <- tryCatch(
      nlsLM(
        kobs ~ emax(CONC, kinact, kI),
        data = out$kobs,
        start = emax_nls_start(out$kobs),
        control = nls.lm.control(maxiter = 1000)
      ),
      error = function(e) {
        warning(
          "Emax fit of kobs over CONC failed: ",
          conditionMessage(e),
          call. = FALSE
        )
        NULL
      }
    )
  }

  if (is.null(fit)) {
    out$tdi_param <- NULL

    out$tdi_plot <- out$kobs |>
      ggplot(aes(x = CONC, y = kobs)) +
      geom_point(size = 2) +
      geom_errorbar(
        aes(ymin = kobs - std.error, ymax = kobs + std.error),
        width = 2) +
      labs(x = "concentration (µM)", y = "kobs (1/min)") +
      theme_bw()
  } else {
    out$tdi_param <- broom::tidy(fit)

    pred <- data.frame(CONC = seq(0, max(out$kobs$CONC), length.out = 100))
    pred$kobs <- predict(fit, newdata = pred)

    out$tdi_plot <- out$kobs |>
      ggplot(aes(x = CONC, y = kobs)) +
      geom_point(size = 2) +
      geom_errorbar(
        aes(ymin = kobs - std.error, ymax = kobs + std.error),
        width = 2) +
      geom_line(data = pred, aes(x = CONC, y = kobs)) +
      labs(x = "concentration (µM)", y = "kobs (1/min)") +
      theme_bw()
  }

  return(out)
}


#' Derive nls starting values for an Emax model of kobs over CONC
#'
#' Linearizes kobs = kinact * CONC / (kI + CONC) as CONC / kobs ~ CONC.
#'
#' @param kobs_df Data frame with CONC and kobs columns.
#' @returns A list with kinact and kI starting values.
#' @noRd
emax_nls_start <- function(kobs_df) {
  fallback <- list(kinact = 0.03, kI = 10)

  kobs_fit <- kobs_df[is.finite(kobs_df$kobs) & kobs_df$kobs > 0, ]
  if (nrow(kobs_fit) < 2)
    return(fallback)

  lin <- lm(I ~ CONC, data = transform(kobs_fit, I = CONC / kobs))
  slope <- coef(lin)[["CONC"]]
  intercept <- coef(lin)[["(Intercept)"]]

  if (!is.finite(slope) || slope <= 0)
    return(fallback)

  kinact_start <- 1 / slope
  kI_start <- intercept * kinact_start

  if (!is.finite(kinact_start) || kinact_start <= 0 ||
      !is.finite(kI_start) || kI_start <= 0)
    return(fallback)

  list(kinact = kinact_start, kI = kI_start)
}
