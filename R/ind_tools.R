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



# temp <- tibble::tribble(
#   ~DONOR, ~CONC, ~OBJECT, ~FOLD, ~PRECIPITANT,
#   "A", 3, "CYP1A2", 1.38, "test",
#   "A", 10, "CYP1A2", 2.07, "test",
#   "A", 30, "CYP1A2", 3.08, "test",
#   "A", 100, "CYP1A2", 7.37, "test",
#   "A", NA, "CYP1A2", 23.8, "positive_control",
#
#   "B", 3, "CYP1A2", 1.24, "test",
#   "B", 10, "CYP1A2", 1.94, "test",
#   "B", 30, "CYP1A2", 2.86, "test",
#   "B", 100, "CYP1A2", 5.09, "test",
#   "B", NA, "CYP1A2", 29, "positive_control",
#
#   "C", 3, "CYP1A2", 1.56, "test",
#   "C", 10, "CYP1A2", 2.25, "test",
#   "C", 30, "CYP1A2", 3.21, "test",
#   "C", 100, "CYP1A2", 4.15, "test",
#   "C", NA, "CYP1A2", 34.9, "positive_control"
# ) |>
#   mutate(SOURCE = "Study examplinib_ind") |>
#   mutate(REL = FOLD / FOLD[PRECIPITANT == "positive_control"] * 100, .by = "DONOR")



indmod <- function(x, precipitant = "test") {
  # input validation
  if (!is.data.frame(x))
    stop("x must be a data frame")

  validate_argument(precipitant, "character")

  expected_cols <- c("DONOR", "PRECIPITANT", "CONC", "OBJECT", "FOLD")
  missing_cols <- setdiff(expected_cols, names(x))
  if (length(missing_cols) > 0)
    stop(paste0("Missing column(s) in input: ", nice_enumeration(missing_cols)))

  if (!precipitant %in% unique(x$PRECIPITANT))
    stop(paste0("Precipitant ", precipitant, " not in data set."))

  # business logic
  x <- mutate(x, ID = paste0(OBJECT, "_", DONOR))

  sigm3 <- function(c, emax, ec50, n) {
    # emax / (1 + exp(-(c - ec50)/n))
    1 + (emax - 1) / (1 + exp(-(c - ec50) / n))
  }

  # hill3 <- function(c, emax, ec50, n) {
  #   1 + (emax - 1) * c^n / (ec50^n +c^n)
  # }

  hill3 <- function(c, emax, ec50, n) {
    1 + (emax) * c^n / (ec50^n +c^n)
  }

  out <- list()

  out$data <- x |>
    filter(PRECIPITANT == precipitant) |>
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
    geom_point(data = filter(x, PRECIPITANT == precipitant, !is.na(FOLD)), size = 2) +
    facet_wrap(~ID, scales = "free") +
    scale_x_log10() +
    theme_bw() +
    theme(legend.position = "none")

  out$ind_param <- out$data |>
    unnest(modpar) |>
    ungroup() |>
    select(-c("ID", "data", "mod"))

    ###############

  temp <- out$data |>
    unnest(modpar) |>
    select(-c("data", "mod"))

  data <- data.frame(
    target = character(0),
    emax = numeric(0),
    ec50 = numeric(0),
    max_c = numeric(0),
    source = character(0)
  )

}



d <- data.frame(
  CONC = 10^seq(1, 3, length.out = 10)
) |>
  mutate(REL = sigmoid3par(CONC, 100, 100)) |>
  mutate(REL = REL * (1 + rnorm(10, 0, .1)))

d |>
  ggplot(aes(x = CONC, y = REL)) +
  geom_point() +
  geom_line() +
  scale_x_log10() +
  theme_bw()

nlsLM(
  REL ~ hill3(CONC, emax, ec50, n),
  data = d,
  start = list(emax = 10, ec50 = 1, n = 1),
  control = nls.lm.control(maxiter = 1000)
)



temp <- x |>
  filter(PRECIPITANT == "test") |>
  # filter(DONOR == "A") |>
  nest_by(DONOR, OBJECT) |>
  mutate(mod = list(
    # lm(log(VALUE) ~ TIME, data = data)
    nlsLM(
      REL ~ hill3(CONC, emax, ec50, n),
      data = data,
      start = list(emax = 10, ec50 = 1, n = 1),
      control = nls.lm.control(maxiter = 1000)
    )
  )) |>
  mutate(modpar = list(broom::tidy(mod))) |>
  unnest(modpar)



pred <- data.frame(CONC = seq(0, 100, 1))
pred$kobs <- predict(fit, newdata = pred)
