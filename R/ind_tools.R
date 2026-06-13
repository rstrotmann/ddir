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
#' @param precipitant The DDI precipitant, defaults to 'test'.
#' @param type The DDI metric, e.g. FOLD or REL, as character. A different type
#' can be given but a column of that name must be present in the data set.
#'
#' @returns A ggplot object.
#' @export
#'
#' @examples
#' induction_plot(examplinib_in_vitro_ind)
#' induction_plot(examplinib_in_vitro_ind, type = "FOLD")
induction_plot <- function(x, precipitant = "test", type = "REL") {
  # input validation
  if (!is.data.frame(x))
    stop("x must be a data frame")

  validate_argument(precipitant, "character")
  validate_argument(type, "character")

  expected_cols <- c("DONOR", "PRECIPITANT", "CONC", "OBJECT", type)
  missing_cols <- setdiff(expected_cols, names(x))
  if (length(missing_cols) > 0)
    stop(paste0("Missing column(s) in input: ", nice_enumeration(missing_cols)))

  if (!precipitant %in% unique(x$PRECIPITANT))
    stop(paste0("Precipitant ", precipitant, " not in data set."))

  # business logic
  x |>
    filter(PRECIPITANT == precipitant) |>
    filter(!is.na(.data[[type]])) |>
    ggplot(aes(x = CONC, y = .data[[type]], group = DONOR)) +
    geom_point() +
    geom_line() +
    facet_wrap(~OBJECT) +
    theme_bw()
}


#' Analyze in vitro CYP induction data
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
#'
#' @returns A list.
#' @export
#'
#' @examples
#' indmod(examplinib_in_vitro_ind, precipitant = "examplinib")
    indmod <- function(x, precipitant = "") {
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
  x <- mutate(x, ID = paste0(OBJECT, "_", DONOR))

  sigm3 <- function(c, emax, ec50, n) {
    1 + (emax - 1) / (1 + exp(-(c - ec50) / n))
  }

  hill3 <- function(c, emax, ec50, n) {
    1 + (emax) * c^n / (ec50^n +c^n)
  }

  out <- list()

  out$data <- x |>
    filter(SAMPLE == "test") |>
    nest_by(DONOR, OBJECT, ID) |>
    mutate(emax_obs = max(data$FOLD, na.rm = TRUE) - 1) |>
    mutate(mod = list(
      nlsLM(
        FOLD ~ hill3(CONC, emax_obs, ec50, n),
        data = data,
        start = list(ec50 = 1, n = 1),
        control = nls.lm.control(maxiter = 1000)
      )
    )) |>
    mutate(modpar = list(broom::tidy(mod)))

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
    facet_wrap(~ID) +
    scale_x_log10() +
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
