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
#'
#' @returns A ggplot object.
#' @export
#'
#' @examples
#' induction_plot(examplinib_in_vitro_ind)
#' induction_plot(examplinib_in_vitro_ind, type = "FOLD")
induction_plot <- function(
    x,
    type = "REL"
) {
  # input validation
  validate_argument(type, values = c("FOLD", "REL"))
  validate_df_argument(x, values = c("DONOR", "SAMPLE", "CONC", "OBJECT", type))
  #
  # if (!is.data.frame(x))
  #   stop("x must be a data frame")
  #
  # # validate_argument(precipitant, "character")
  # validate_argument(type, "character")
  #
  # expected_cols <- c("DONOR", "SAMPLE", "CONC", "OBJECT", type)
  # missing_cols <- setdiff(expected_cols, names(x))
  # if (length(missing_cols) > 0)
  #   stop(paste0("Missing column(s) in input: ", nice_enumeration(missing_cols)))

  # business logic
  p <- x |>
    filter(SAMPLE == "test") |>
    filter(!is.na(.data[[type]])) |>
    ggplot(aes(x = CONC, y = .data[[type]], group = DONOR)) +
    geom_point() +
    geom_line() +
    facet_wrap(~OBJECT) +
    theme_bw()

  p
}


#' Analyze in vitro CYP induction data
#'
#' Fit a 3-parameter hill function to the mRNA induction data from each
#' donor and DDI object.
#'
#' @details
#'
#' \deqn{f = 1 + \frac{(E_{max}) * C^n}{(EC_{50}^n + C^n)}}
#'
#' @param x In vitro induction data as data frame with the columns:
#' * DONOR The hepatocyte donor
#' * SAMPLE "test" or "positive_control"
#' * CONC The preipitant concentration in uM.
#' * OBJECT The DDI object.
#' * FOLD The mRNA or activity fold change.
#' * REL The fold change relative to positive control.
#' * SOURCE Source information
#' @param precipitant The precipitant as character.
#' @param use_emax_obs Constrain Emax to the observed Emax (default). If FALSE,
#' Emax will be fitted.
#'
#' @returns A list with the elements:
#' * data The original data, model and model parameters for each fit.
#' * fold_plot A ggplot object showing the original data and the fits.
#' * ind_param The ec50, hill (n) and emax parameters by donor and DDI object.
#' * inducer An inducer object representing the model parameters from the donor
#' with the highest Emax per DDI object.
#'
#' @export
#'
#' @examples
#' indmod(examplinib_in_vitro_ind, precipitant = "examplinib")
indmod <- function(x, precipitant = "", use_emax_obs = TRUE) {
  # input validation
  if (!is.data.frame(x))
    stop("x must be a data frame")

  validate_argument(precipitant, "character")

  expected_cols <- c("DONOR", "SAMPLE", "CONC", "OBJECT", "FOLD")
  missing_cols <- setdiff(expected_cols, names(x))
  if (length(missing_cols) > 0)
    stop(paste0("Missing column(s) in input: ", nice_enumeration(missing_cols)))

  if (any(!x$SAMPLE %in% c("test", "positive_control")))
    stop("SAMPLE must be 'test' or 'positive_control'")

  # business logic
  x <- mutate(x, ID = paste0(.data$OBJECT, "_", .data$DONOR))

  sigm3 <- function(c, emax, ec50, n) {
    1 + (emax - 1) / (1 + exp(-(c - ec50) / n))
  }

  hill3 <- function(c, emax, ec50, n) {
    1 + (emax - 1) * c^n / (ec50^n + c^n)
  }

  out <- list()

  # out$data <- x |>
  #   filter(.data$SAMPLE == "test") |>
  #   nest_by(DONOR, OBJECT, ID) |>
  #   mutate(emax_obs = max(data$FOLD, na.rm = TRUE) - 1) |>
    # mutate(mod = list(
    #   nlsLM(
    #     FOLD ~ hill3(CONC, emax, ec50, n),
    #     data = data,
    #     start = list(emax = 2, ec50 = .1, n = 1),
    #     lower = c(emax = NA, ec50 = 0, n = 1),
    #     upper = c(emax = NA, ec50 = 100, n = 5),
    #     control = nls.lm.control(maxiter = 1000)
    #   )
    # )) |>
    # mutate(modpar = list(broom::tidy(mod)))

  temp <- x |>
    filter(.data$SAMPLE == "test") |>
    nest_by(DONOR, OBJECT, ID) |>
    mutate(emax_obs = max(data$FOLD, na.rm = TRUE) - 1)

  if (use_emax_obs == TRUE) {
    out$data <- temp |>
      mutate(mod = list(
        nlsLM(
          FOLD ~ hill3(CONC, emax_obs, ec50, n),
          data = data,
          start = list(ec50 = .1, n = 1),
          lower = c(ec50 = 0, n = 1),
          upper = c(ec50 = 100, n = 5),
          control = nls.lm.control(maxiter = 1000)
        )
      )) |>
      mutate(modpar = list(broom::tidy(mod)))
  } else {
    out$data <- temp |>
      mutate(mod = list(
        nlsLM(
          FOLD ~ hill3(CONC, emax, ec50, n),
          data = data,
          start = list(emax = 2, ec50 = .1, n = 1),
          lower = c(emax = NA, ec50 = 0, n = 1),
          upper = c(emax = NA, ec50 = 100, n = 5),
          control = nls.lm.control(maxiter = 1000)
        )
      )) |>
      mutate(modpar = list(broom::tidy(mod)))
  }



  # curve plot data set
  pred <- data.frame(
    CONC = 10^seq(
      log10(min(x$CONC, na.rm = TRUE)),
      log10(max(x$CONC, na.rm = TRUE)),
      length.out = 100)
  )

  temp <- lapply(out$data$mod, function(x) predict(x, newdata = pred))
  names(temp) <- out$data$ID
  pred <- bind_cols(pred, temp) |>
    pivot_longer(cols = -1, names_to = "ID", values_to = "FOLD") |>
    separate(ID, c("OBJECT", "DONOR"), "_", remove = FALSE)

  out$fold_plot <- ggplot(data = NULL, aes(x = CONC, y = FOLD, color = OBJECT)) +
    geom_line(data = pred) +
    geom_point(data = filter(x, SAMPLE == "test", !is.na(FOLD)), size = 2) +
    facet_wrap(~ID, scales = "free") +
    scale_x_log10() +
    expand_limits(y = 0) +
    theme_bw() +
    theme(legend.position = "none")

  out$ind_param <- out$data |>
    unnest(modpar) |>
    ungroup() |>
    select(-c("ID", "data", "mod"))

  # make inducer object from the donor that has the respective highers emax
  max_c <- x |>
    reframe(
      max_c = max(CONC, na.rm = TRUE), .by = c("ID", "SOURCE"))

  temp <- out$data |>
    unnest(modpar) |>
    select(-c("data", "mod")) |>
    group_by(OBJECT) |>
    filter(emax_obs == max(emax_obs, na.rm = TRUE)) |>
    filter(term == "ec50") |>
    left_join(max_c, by = "ID") |>
    select(object = OBJECT, emax = emax_obs, ec50 = estimate, max_c, source = SOURCE) |>
    mutate(ec50 = round(ec50, 2))

  out$inducer <- inducer(temp, precipitant = precipitant)

  out

}
